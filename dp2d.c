#include <math.h>
#include <search.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "dp2d.h"
#include "fp.h"
#include "globals.h"
#include "histogram.h"
#include "itstree.h"
#include "rs.h"

#define MICROSECONDS 1000000L

/* scale factor for noise */
#ifndef SCALE_FACTOR
#define SCALE_FACTOR 2.3 /* log(10) = 2.3, 90% of noise */
#endif
/* print the noisy counts for each item */
#ifndef PRINT_ITEM_TABLE
#define PRINT_ITEM_TABLE 0
#endif
/* print the returned rules */
#ifndef PRINT_FINAL_RULES
#define PRINT_FINAL_RULES 0
#endif
/* asymmetric quality function */
#ifndef ASYMMETRIC_Q
#define ASYMMETRIC_Q 0
#endif
/* use EM to select first item too, instead of noisy count */
#ifndef EM_1ST_ITEM
#define EM_1ST_ITEM 0
#endif
/* force last selection quality */
#ifndef EM_FORCED_LAST
#define EM_FORCED_LAST 0
#endif
/* use only last item added to itemset */
#ifndef EM_LAST_ITEM
#define EM_LAST_ITEM 1
#endif
/* Ensemble reduce function (min/max), valid only if EM_LAST_ITEM is 0 */
#ifndef EM_REDFUN
#define EM_REDFUN max
#endif

enum quality_fun {
	EM_QD = 0,
	EM_QDELTA,
	EM_QSIGMA
};
static const enum quality_fun QMETHOD = EM_QSIGMA;

struct item_count {
	int value;
	int real_count;
	double noisy_count;
};

static int ic_noisy_cmp(const void *a, const void *b)
{
	const struct item_count *ia = a, *ib = b;
	return double_cmp_r(&ia->noisy_count, &ib->noisy_count);
}

static size_t build_items_table(const struct fptree *fp, struct item_count *ic,
		double eps, struct drand48_data *buffer)
{
	size_t i;

	printf("Compute noisy counts for items with eps = %lf\n", eps);
	for (i = 0; i < fp->n; i++) {
		ic[i].value = i + 1;
		ic[i].real_count = fpt_item_count(fp, i);
		ic[i].noisy_count = laplace_mechanism(ic[i].real_count, eps,
				1, buffer);
		if (ic[i].noisy_count < 0)
			ic[i].noisy_count = 0;
	}

	qsort(ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);

#if PRINT_ITEM_TABLE
	print_item_table(ic, fp->n);
#endif

	printf("Noise scale: %5.2f\n", SCALE_FACTOR/eps);
	for (i = 0; i < fp->n; i++)
		if (ic[i].noisy_count < SCALE_FACTOR / eps)
			return i;

	return fp->n;
}

#if PRINT_FINAL_RULES
static void print_this_rule(const int *A, const int* AB,
		size_t a_length, size_t ab_length, double c)
{
	size_t i, j;

	for (i = 0; i < a_length; i++)
		printf("%d ", A[i]);
	printf("-> ");
	for (i = 0; i < ab_length; i++) {
		for (j = 0; j < a_length; j++)
			if (AB[i] == A[j])
				j = 2 * a_length;
		if (j == a_length)
			printf("%d ", AB[i]);
	}
	printf("| c=%7.6f\n", c);
}
#endif

/**
 * Checks whether the current itemset has been generated previously
 */
static int its_already_seen(int *its, size_t itslen,
		const struct itstree_node *itst)
{
	int *cf = calloc(itslen, sizeof(cf[0]));
	size_t i, ret = 0;

	for (i = 0; i < itslen; i++)
		cf[i] = its[i];

	qsort(cf, itslen, sizeof(cf[0]), int_cmp);
	ret = search_its_private(itst, cf, itslen);

	free(cf);
	return ret;
}

/**
 * Updates the list of itemsets that were generated.
 */
static void update_seen_its(const int *its, size_t itslen,
		struct itstree_node *itst)
{
	int *cf = calloc(itslen, sizeof(cf[0]));
	size_t i;

	for (i = 0; i < itslen; i++)
		cf[i] = its[i];

	qsort(cf, itslen, sizeof(cf[0]), int_cmp);
	record_its_private(itst, cf, itslen);
	free(cf);
}

static void generate_rules_from_itemset(const int *AB, size_t ab_length,
		const struct fptree *fp, double *minc, double *maxc,
		struct histogram *h)
{
	int *A = calloc(ab_length, sizeof(A[0]));
	size_t i, j, max, a_length;
	int sup_ab, sup_a;
	double c;

	max = (1 << ab_length) - 1;
	sup_ab = fpt_itemset_count(fp, AB, ab_length);
	for (i = 1; i < max; i++) {
		a_length = 0;
		for (j = 0; j < ab_length; j++)
			if (i & (1 << j))
				A[a_length++] = AB[j];

		sup_a = fpt_itemset_count(fp, A, a_length);
		c = div_or_zero(sup_ab, sup_a);
		if (c < *minc) *minc = c;
		if (c > *maxc) *maxc = c;
		histogram_register(h, c);

#if PRINT_FINAL_RULES
		print_this_rule(A, AB, a_length, ab_length, c);
#endif
	}

	free(A);
}

static void generate_rules(const int *items, size_t lmax,
		const struct fptree *fp,
		double *minc, double *maxc, struct histogram *h,
		struct itstree_node *itst)
{
	int *AB = calloc(lmax, sizeof(AB[0]));
	size_t i, j, max=1<<lmax, ab_length;

	for (i = 0; i < max; i++) {
		ab_length = 0;
		for (j = 0; j < lmax; j++)
			if (i & (1 << j))
				AB[ab_length++] = items[j];
		if (ab_length < 2)
			continue;
		if (its_already_seen(AB, ab_length, itst))
			continue;
		update_seen_its(AB, ab_length, itst);
		generate_rules_from_itemset(AB, ab_length, fp, minc, maxc, h);
	}

	free(AB);
}

struct reservoir_item {
	int *items;
	size_t sz;
	int support;
	double q;
};

static void print_reservoir_item(const void *it)
{
	const struct reservoir_item *ri = it;
	size_t i;

	printf("[%3d", ri->items[0]);
	for (i = 1; i < ri->sz; i++)
		printf(", %3d", ri->items[i]);
	printf("], s=%5d, q=%7.2lf", ri->support, ri->q);
}

static void *clone_reservoir_item(const void *it)
{
	struct reservoir_item *ret = calloc(1, sizeof(*ret));
	const struct reservoir_item *ri = it;
	size_t i;

	ret->items = calloc(ri->sz, sizeof(ret->items[0]));
	for (i = 0; i < ri->sz; i++)
		ret->items[i] = ri->items[i];
	ret->sz = ri->sz;
	ret->support = ri->support;
	ret->q = ri->q;

	return ret;
}

static void free_reservoir_item(void *it)
{
	struct reservoir_item *ri = it;
	free(ri->items);
	free(ri);
}

static inline double quality_d(int x, int y, double c0)
{
	double q = -x + y / c0;
#if ASYMMETRIC_Q
	if (q > 0) q = 0;
#endif
	return -fabs(q);
}

static inline double compute_d_quality(const struct fptree *fp,
		double c0, int sup_ab, struct reservoir_item *rit)
{
	double bq = quality_d(fpt_item_count(fp, rit->items[rit->sz - 1]),
			sup_ab, c0);
#if !EM_LAST_ITEM
	size_t i, ep = rit->sz - 1;

	for (i = 0; i < ep; i++)
		bq = EM_REDFUN(bq,
			quality_d(fpt_item_count(fp, rit->items[i]),
				sup_ab, c0));
#endif

	return bq;
}

static inline double compute_delta_quality(const struct fptree *fp,
		int sup_ab, struct reservoir_item *rit)
{
	double bq = sup_ab - fpt_itemset_count(fp, rit->items, rit->sz - 1);

#if !EM_LAST_ITEM
	size_t i, j, ep = rit->sz - 1;
	int t;

	for (i = 0; i < rit->sz; i++) {
		t = rit->items[ep];
		rit->items[ep] = rit->items[i];
		rit->items[i] = t;
		for (j = 1; j < ep; j++)
			bq = EM_REDFUN(bq, sup_ab - fpt_itemset_count(fp,
						rit->items, j));
		t = rit->items[ep];
		rit->items[ep] = rit->items[i];
		rit->items[i] = t;
	}
#endif

	return bq;
}

static inline double compute_quality(const struct fptree *fp, double c0,
		const struct item_count *ic, size_t ix_item,
		struct reservoir_item *rit, size_t lmax)
{
	int sup_ab = fpt_itemset_count(fp, rit->items, rit->sz);
	(void)lmax; /* used only if EM_FORCED_LAST */

	/* select first item: use either real or noisy count */
	if (rit->sz == 1)
#if EM_1ST_ITEM
		return ic[ix_item].real_count;
#else
		return ic[ix_item].noisy_count;
#endif

#if EM_FORCED_LAST
	if (rit->sz == lmax)
		return compute_d_quality(fp, c0, sup_ab, rit);
#endif

	switch(QMETHOD) {
	case EM_QD: return compute_d_quality(fp, c0, sup_ab, rit);
	case EM_QDELTA: return compute_delta_quality(fp, sup_ab, rit);
	default: return sup_ab;
	}
}

static inline int generated_above(const int *celms, size_t level)
{
	size_t i;

	for (i = 0; i < level; i++)
		if (celms[i] == celms[level])
			return 1;
	return 0;
}

static void mine_level(const struct fptree *fp, const struct item_count *ic,
		size_t numits, size_t lmax, const int *celms, size_t level,
		double c0, double *epss, size_t *spls, struct histogram *h,
		double *minc, double *maxc, struct itstree_node *itst,
		struct drand48_data *randbuffer)
{
	struct reservoir_item *rit = calloc(1, sizeof(*rit));
	const struct reservoir_item *crit;
	struct reservoir_iterator *ri;
	struct reservoir *r;
	double eps_round;
	size_t i;

	r = init_reservoir(spls[level], print_reservoir_item,
			clone_reservoir_item, free_reservoir_item);
	eps_round = epss[level] / spls[level];

	/* init common part of rit */
	rit->sz = level + 1;
	rit->items = calloc(rit->sz, sizeof(rit[0]));
	for (i = 0; i < level; i++)
		rit->items[i] = celms[i];

	/* generate last element */
	for (i = 0; i < numits; i++) {
		rit->items[level] = ic[i].value;
		if (generated_above(rit->items, level))
			continue;
		if (level == lmax - 1 &&
				its_already_seen(rit->items, lmax, itst))
			continue;

		rit->support = fpt_itemset_count(fp, rit->items, rit->sz);
		rit->q = compute_quality(fp, c0, ic, i, rit, lmax);
		add_to_reservoir_log(r, rit, eps_round * rit->q/2, randbuffer);
	}
	free_reservoir_item(rit);

	ri = init_reservoir_iterator(r);
	/* TODO: generate all subtrees after a level? */
	if (level == lmax - 1)
		while ((crit = next_item(ri)))
			generate_rules(crit->items, lmax, fp, minc, maxc, h,
					itst);
	else while ((crit = next_item(ri)))
		mine_level(fp, ic, numits, lmax, crit->items, level + 1, c0,
				epss, spls, h, minc, maxc, itst, randbuffer);
	free_reservoir_iterator(ri);
	free_reservoir(r);
}

static void print_mining_scenario()
{
	size_t i;
	enum quality_fun qf = QMETHOD;

	printf("Methods used: ");
#if EM_1ST_ITEM
	printf("em ");
#else
	printf("noisy ");
#endif

	for (i = 0; i < 2; i++) {
#if EM_FORCED_LAST
		if (i) qf = EM_QD;
#endif
#if !EM_LAST_ITEM
		printf("m%c%c(", -EM_REDFUN(-'i',-'a'), EM_REDFUN('n', 'x'));
#endif
		switch(qf) {
		case EM_QD: printf("qd"); break;
		case EM_QDELTA: printf("qdelta"); break;
		default: printf("qsigma");
		}
#if !EM_LAST_ITEM
		printf(")");
#endif
		if (i) printf("\n");
		else printf(" ");
	}
}

/**
 * Step 2 of mining, private.
 */
static void mine_rules(const struct fptree *fp, const struct item_count *ic,
		double eps, double c0, size_t numits, size_t lmax,
		struct histogram *h, double *minc, double *maxc, size_t cspl,
		struct drand48_data *randbuffer)
{
	double *epsilons = calloc(lmax, sizeof(epsilons[0]));
	struct itstree_node *itst = init_empty_itstree();
	size_t *spl = calloc(lmax, sizeof(spl[0]));
	size_t i, f = 1;
	double cf = 0;

	printf("Mining with eps %lf, numitems=%lu\n", eps, numits);
	print_mining_scenario();

#if !EM_1ST_ITEM
	cf = 1;
#endif
	for (i = 0; i < lmax; i++) {
		/* TODO: better formulas here */
		epsilons[i] = eps/(lmax - cf);
		spl[i] = cspl;

		epsilons[i] /= f; /* branching factor */
		f *= spl[i];
	}
#if !EM_1ST_ITEM
	epsilons[0] = spl[0] * 2; /* use noisy count */
#endif
	printf("Total leaves %lu\n", f);

	mine_level(fp, ic, numits, lmax, NULL, 0, c0, epsilons, spl, h,
			minc, maxc, itst, randbuffer);

	free(epsilons);
	free_itstree(itst);
	free(spl);
}

#if PRINT_ITEM_TABLE
static inline void print_item_table(const struct item_count *ic, size_t n)
{
	size_t i;

	printf("\n");
	for (i = 0; i < n; i++)
		printf("%5lu[%5.2lf] %5d %7d %9.2lf\n", i, (i + 1.0)/n,
				ic[i].value, ic[i].real_count,
				ic[i].noisy_count);
}
#endif

void dp2d(const struct fptree *fp, double eps, double eps_ratio1,
		double c0, size_t lmax, size_t ni, long int seed, size_t cspl)
{
	struct item_count *ic = calloc(fp->n, sizeof(ic[0]));
	double epsilon_step1 = eps * eps_ratio1;
	struct histogram *h = init_histogram();
	struct timeval starttime, endtime;
	struct drand48_data randbuffer;
	double minc, maxc, t1, t2;
	size_t numits;

	printf("eps=%lf, eps_step1=%lf, c0=%5.2lf, rmax=%lu\n",
			eps, epsilon_step1, c0, lmax);

	init_rng(seed, &randbuffer);
	build_items_table(fp, ic, epsilon_step1, &randbuffer);
	minc = 1;
	maxc = 0;
	numits = min(ni, fp->n);
	eps = eps - epsilon_step1;

	gettimeofday(&starttime, NULL);
	mine_rules(fp, ic, eps, c0, numits, lmax, h, &minc, &maxc, cspl,
			&randbuffer);
	gettimeofday(&endtime, NULL);
	t1 = starttime.tv_sec + (0.0 + starttime.tv_usec) / MICROSECONDS;
	t2 = endtime.tv_sec + (0.0 + endtime.tv_usec) / MICROSECONDS;

	printf("Rules saved: %lu, minconf: %3.2lf, maxconf: %3.2lf\n",
			histogram_get_all(h), minc, maxc);
	printf("Total time: %5.2lf\n", t2 - t1);
	printf("%ld %ld %ld %ld\n", starttime.tv_sec, starttime.tv_usec,
			endtime.tv_sec, endtime.tv_usec);

	printf("Final histogram:\n");
	histogram_dump(stdout, h, 1, "\t");

	free_histogram(h);
	free(ic);
}
