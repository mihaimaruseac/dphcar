#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "dp2d.h"
#include "fp.h"
#include "globals.h"
#include "histogram.h"
#include "rule.h"
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
/* print the rule lattice generation step debug info */
#ifndef PRINT_RULE_LATTICE
#define PRINT_RULE_LATTICE 0
#endif
/* print changes to the selected rule */
#ifndef PRINT_RULE_LATTICE_TRACE
#define PRINT_RULE_LATTICE_TRACE 0
#endif
/* print the returned rules */
#ifndef PRINT_FINAL_RULES
#define PRINT_FINAL_RULES 0
#endif
/* asymmetric quality function */
#ifndef ASYMMETRIC_Q
#define ASYMMETRIC_Q 0
#endif
/* clique cutoff */
#ifndef CLIQUE_CUTTOF
#define CLIQUE_CUTTOF 0
#endif
/* uniform ("proper") reduction method */
#ifndef UNIFORM_REDUCTION
#define UNIFORM_REDUCTION 1
#endif

static double quality(int x, int y, double c0, struct drand48_data *buffer)
{
	(void)buffer;
	double q = -x + y / c0;
#if ASYMMETRIC_Q
	if (q > 0) q = 0;
#endif
	return -fabs(q);
}

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
		double eps, struct drand48_data *buffer, int private)
{
	size_t i;

	for (i = 0; i < fp->n; i++) {
		ic[i].value = i + 1;
		ic[i].real_count = fpt_item_count(fp, i);
		if (private) {
			ic[i].noisy_count = laplace_mechanism(
					ic[i].real_count, eps, 1, buffer);
			if (ic[i].noisy_count < 0)
				ic[i].noisy_count = 0;
		} else
			ic[i].noisy_count = ic[i].real_count;
	}

	qsort(ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);

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
static inline int its_already_seen(int *its, size_t itslen,
		const int *seen, size_t seenlen)
{
	size_t ix = 0, i, j;

	qsort(its, itslen, sizeof(its[0]), int_cmp);
	while (ix < seenlen) {
		/* check length first */
		if (seen[ix] != (int)itslen)
			goto nothere;

		/* same length, check by items */
		for (i = 0, j = ix + 1; i < itslen; i++, j++)
			if (its[i] != seen[j])
				goto nothere;

		return 1;
nothere:
		ix += seen[ix] + 1;
	}

	return 0;
}

/**
 * Updates the list of itemsets that were generated.
 */
static inline void update_seen_its(const int *its, size_t itslen,
		int *seen, size_t *seenlen)
{
	size_t ix = *seenlen, i;

	seen[ix++] = itslen;
	for (i = 0; i < itslen; i++)
		seen[ix++] = its[i];

	*seenlen = ix;
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
	for (i = 1; i < max; i++) {
		a_length = 0;
		for (j = 0; j < ab_length; j++)
			if (i & (1 << j))
				A[a_length++] = AB[j];

		sup_ab = fpt_itemset_count(fp, AB, ab_length);
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

static void generate_rules(const size_t *items, size_t lmax,
		const struct fptree *fp, const struct item_count *ic,
		double *minc, double *maxc, struct histogram *h,
		int *seen, size_t *seenix)
{
	int *AB = calloc(lmax, sizeof(AB[0]));
	size_t i, j, max=1<<lmax, ab_length;

	for (i = 0; i < max; i++) {
		ab_length = 0;
		for (j = 0; j < lmax; j++)
			if (i & (1 << j))
				AB[ab_length++] = ic[items[j]].value;
		if (ab_length < 2)
			continue;
		if (its_already_seen(AB, ab_length, seen, *seenix))
			continue;
		update_seen_its(AB, ab_length, seen, seenix);
		generate_rules_from_itemset(AB, ab_length, fp, minc, maxc, h);
	}

	free(AB);
}

/**
 * Analyze the current items to see if we can select a good rule lattice.
 */
static void analyze_items(const size_t *items, size_t lmax,
		struct reservoir *reservoir,
		const struct fptree *fp, const struct item_count *ic,
		double c0, double eps, struct drand48_data *randbuffer)
{
	size_t *citems = calloc(lmax, sizeof(citems[0]));
	int *AB = calloc(lmax, sizeof(AB[0]));
	int sup_ab, sup_a;
	size_t i, j, k;
	double q;

	for (i = 0; i < lmax; i++)
		AB[i] = ic[items[i]].value;

	sup_ab = fpt_itemset_count(fp, AB, lmax);

#if PRINT_RULE_LATTICE
	printf("Analyzing new set of items: ");
	for (i = 0; i < lmax; i++)
		printf("%3d ", AB[i]);
	printf(" | support: %d\n", sup_ab);
#endif

	/* try for each corner item */
	for (i = 0; i < lmax; i++) {
		sup_a = ic[items[i]].real_count;
		q = quality(sup_a, sup_ab, c0, randbuffer);
		/* need a copy to insert in reservoir */
		for (j = 0; j < lmax; j++) citems[j] = items[j];
		k = citems[0]; citems[0] = citems[i]; citems[i] = k;
		add_to_reservoir_log(reservoir, citems, lmax,
				sizeof(citems[0]), eps * q / 2, randbuffer);
	}

	free(citems);
	free(AB);
}

/**
 * Checks whether the current items vector is forbidden (already generated).
 */
static inline int already_seen(const size_t *items, size_t lmax,
		const size_t *seen, size_t seenlen)
{
	size_t i, j, ix = 0;

	while (ix < seenlen) {
		for (i = 0, j = ix; i < lmax; i++, j++)
			if (items[i] != seen[j]) {
				ix += lmax;
				break;
			}

		if (i == lmax)
			return 1;
	}

	return 0;
}

/**
 * Constructs the next items vector, the next set of rules to be analyzed.
 */
static int update_items(const struct item_count *ic, double c0, size_t *items,
		size_t lmax, size_t n, const size_t *seen, size_t seenlen,
		size_t rf)
{
	size_t ix = lmax - 1, ok;
#if CLIQUE_CUTTOF
	double t = c0 * ic[items[0]].noisy_count;
#else
	(void)ic;
	(void)c0;
#endif

	do {
		do {
			/* try next */
			ok = 1;
#if UNIFORM_REDUCTION
			if (ix == lmax - 1) { /* here we can skip */
#endif
				items[ix] += rf;
#if UNIFORM_REDUCTION
				rf = items[ix] - n + 1;
			} else
				items[ix]++;
#endif

			/* if impossible or undesirable, move to one below */
			if (items[ix] >= n
#if CLIQUE_CUTTOF
					|| ic[items[ix]].noisy_count < t
#endif
					) {
				if (!ix) return 1; /* end of generation */
				ix--;
				ok = 0; /* this level was not finished */
				break; /* restart process for ix - 1 */
			}

			/* check to not generate seen set */
			if (ix == lmax - 1 &&
				already_seen(items, lmax, seen, seenlen)) {
#if UNIFORM_REDUCTION
				rf = 1;
#endif
				ok = 0; /* get next set */
			}
		} while (!ok);

		if (ok) {
			ix++; /* move to next */
			if (ix == lmax)
				return 0; /* all generated */
			items[ix] = items[ix - 1]; /* before first option */
		}
	} while (ix < lmax);

	return 0; /* notreached */
}

/**
 * Initialize the items vector, the first set of rules to be analyzed.
 */
static inline int init_items(const struct item_count *ic, double c0,
		size_t *items, size_t lmax, size_t n,
		const size_t *seen, size_t seenlen)
{
	size_t i;

	for (i = 0; i < lmax; i++)
		items[i] = i;

	if (already_seen(items, lmax, seen, seenlen))
		return update_items(ic, c0, items, lmax, n, seen, seenlen, 1);

	return 0;
}

/**
 * Step 2 of mining, private.
 */
static void mine_rules(const struct fptree *fp, const struct item_count *ic,
		struct histogram *h, size_t numits, size_t lmax,
		double *minc, double *maxc, double c0, size_t k, double eps,
		struct drand48_data *randbuffer, size_t rf)
{
	struct reservoir_iterator *reservoir_iterator;
	size_t i, j, end, seenix, seenitsix;
	struct timeval starttime, endtime;
	size_t *items, *seen;
	const size_t *bitems;
	struct reservoir *reservoir;
	double t1, t2;
	int *seenits;

	printf("Mining %lu steps each with eps %lf, numitems=%lu\n",
			k, eps, numits);

	items = calloc(lmax, sizeof(items[0]));
	seenits = calloc(k * (lmax+ 1) * (1 << lmax), sizeof(seenits[0]));
	seen = calloc(k * lmax, sizeof(seen[0]));
	seenix = 0;
	seenitsix = 0;

	for (i = 0; i < k; i++) {
		gettimeofday(&starttime, NULL);
		t1 = starttime.tv_sec + (0.0 + starttime.tv_usec) / MICROSECONDS;
		end = init_items(ic, c0, items, lmax, numits, seen, seenix);
		if (end) break;
		bitems = calloc(lmax, sizeof(bitems[0]));

		reservoir = init_reservoir(1, print_size_t_array,
				shallow_clone, free);
		while (!end) {
			analyze_items(items, lmax, reservoir,
					fp, ic, c0, eps, randbuffer);
			end = update_items(ic, c0, items, lmax, numits,
					seen, seenix, rf);
		}

		reservoir_iterator = init_reservoir_iterator(reservoir);
		bitems = next_item(reservoir_iterator, NULL, NULL);

#if PRINT_RULE_LATTICE_TRACE || PRINT_FINAL_RULES
		printf("Selected items: %lu(%d) -> ",
				bitems[0], ic[bitems[0]].value);
		for (j = 1; j < lmax; j++)
			printf("%lu(%d) ", bitems[j],
					ic[bitems[j]].value);
		printf("\n");
#endif

		generate_rules(bitems, lmax, fp, ic, minc, maxc, h,
				seenits, &seenitsix);

		/* remove bitems from future items */
		qsort((void*)bitems, lmax, sizeof(bitems[0]), int_cmp);
		for (j = 0; j < lmax; j++)
			seen[seenix++] = bitems[j];

		free_reservoir_iterator(reservoir_iterator);
		free_reservoir(reservoir);

		gettimeofday(&endtime, NULL);
		t2 = endtime.tv_sec + (0.0 + endtime.tv_usec) / MICROSECONDS;
		printf("Round %lu time: %5.2lf\n", i, t2 - t1);
	}

	free(seenits);
	free(items);
	free(seen);
}

/**
 * Step 2 of mining, not private.
 */
static void mine_rules_np(const struct fptree *fp, const struct item_count *ic,
		struct histogram *h, size_t numits, size_t lmax, double c0,
		double *minc, double *maxc)
{
	size_t *items = calloc(lmax, sizeof(items[0]));
	int *AB = calloc(lmax, sizeof(AB[0]));
	size_t clen, i, end;

	printf("Mining all rules of top %lu items\n", numits);

	for (clen = 2; clen <= lmax; clen++) {
		end = init_items(ic, c0, items, clen, numits, NULL, 0);
		while (!end) {
			for (i = 0; i < clen; i++)
				AB[i] = ic[items[i]].value;
			generate_rules_from_itemset(AB, clen, fp,
					minc, maxc, h);
			end = update_items(ic, c0, items, clen, numits,
					NULL, 0, 1);
		}
	}

	free(AB);
	free(items);
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

#if CLIQUE_CUTTOF
static void count_cliques(const struct fptree *fp, const struct item_count *ic,
		double c0, size_t lmax)
{

	size_t i, j;
	double cliques = 0;
	double t;

	for (i = lmax-1; i < fp->n; i++) {
		t = c0 * ic[i-lmax+1].noisy_count;
		if (t == 0) break;
		for (j = i; j < fp->n; j++)
			if (ic[j].noisy_count < t)
				break;
		cliques += pow((j - i), lmax);
	}
	printf("Graph would have %lf cliques.\n", cliques);
}
#endif

void dp2d(const struct fptree *fp, double eps, double eps_ratio1,
		double c0, size_t lmax, size_t k, long int seed, int private,
		size_t ni, size_t rf)
{
	struct item_count *ic = calloc(fp->n, sizeof(ic[0]));
	double epsilon_step1 = eps * eps_ratio1;
	struct histogram *h = init_histogram();
	struct timeval starttime, endtime;
	struct drand48_data randbuffer;
	double minc, maxc, t1, t2;
	size_t numits;

	init_rng(seed, &randbuffer);

	if (private) {
		printf("Running private method with eps=%lf, eps_step1=%lf, "
				"k=%lu, c0=%5.2lf, rmax=%lu\n", eps,
				epsilon_step1, k, c0, lmax);
		printf("Compute noisy counts for items with "
				"eps_1 = %lf\n", epsilon_step1);
	} else {
		printf("Running non-private method with k=%lu, c0=%5.2lf, "
				"rmax=%lu\n", k, c0, lmax);
	}

	numits = build_items_table(fp, ic, epsilon_step1, &randbuffer, private);
#if PRINT_ITEM_TABLE
	print_item_table(ic, fp->n);
#endif
#if CLIQUE_CUTTOF
	count_cliques(fp, ic, c0, lmax);
#endif

	minc = 1;
	maxc = 0;
	numits = ni;
	if (numits > fp->n) numits = fp->n;
	eps = eps - epsilon_step1;
	gettimeofday(&starttime, NULL);
	if (private)
		mine_rules(fp, ic, h, numits, lmax, &minc, &maxc, c0, k,
				eps / k, &randbuffer, rf);
	else
		mine_rules_np(fp, ic, h, numits, lmax, c0, &minc, &maxc);
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
