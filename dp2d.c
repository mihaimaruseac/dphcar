#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include <gmp.h>
#include <mpfr.h>

#include "dp2d.h"
#include "fp.h"
#include "globals.h"
#include "histogram.h"
#include "rule.h"

#define MICROSECONDS 1000000L

/* print the noisy counts for each item */
#ifndef PRINT_ITEM_TABLE
#define PRINT_ITEM_TABLE 0
#endif
/* print the rules generated at each step and their quality */
#ifndef PRINT_RULE_DOMAIN
#define PRINT_RULE_DOMAIN 0
#endif
/* print actions to the reservoir */
#ifndef PRINT_RS_TRACE
#define PRINT_RS_TRACE 0
#endif
/* print the returned rules */
#ifndef PRINT_FINAL_RULES
#define PRINT_FINAL_RULES 0
#endif
/* print probabilities */
#ifndef PRINT_PROBS
#define PRINT_PROBS 0
#endif

static double quality(size_t x, size_t y, size_t m, size_t sfactor, double c0)
{
	double c;

#if 0 /* original quality fun*/
	if (x < m) x = m;
	c = (y + 0.0) / (x + 0.0);
	c -= c0;
	if (c < 0) c = 0;

	return c * m / sfactor;
#else
	(void)m;
#if 0 /* cosine quality */
	c = x + c0 * y;
	return -c / (1 + c0) / sfactor;
#else /* sine quality */
	c = y - c0 * x;
	if (c < 0) c = 0;

	return c / sfactor;
#endif
#endif
}

struct item_count {
	size_t value;
	size_t real_count;
	double noisy_count;
	double smin;
	double smax;
};

struct reservoir {
	struct rule *r;
	double c;
	double v;
	int sxy;
	int sx;
};

static void free_reservoir_array(struct reservoir *reservoir, int size)
{
	int i;

	for (i = 0 ; i < size; i++)
		free_rule(reservoir[i].r);
	free(reservoir);
}

static int reservoir_cmp(const void *a, const void *b)
{
	const struct reservoir *ra = a, *rb = b;
	return double_cmp(&ra->v, &rb->v);
}

static void build_items_table(const struct fptree *fp, struct item_count *ic,
		double eps, struct drand48_data *buffer)
{
	double scale = log(10) * fp->l_max_t / eps;
	size_t i;

#if PRINT_PROBS
	printf("Scale parameter: %lf\n", scale);
#endif

	for (i = 0; i < fp->n; i++) {
		ic[i].value = i + 1;
		ic[i].real_count = fpt_item_count(fp, i + 1);
		ic[i].noisy_count = laplace_mechanism(ic[i].real_count,
				eps, fp->l_max_t, buffer);
		/* TODO: filter some noise */
		if (ic[i].noisy_count < 0)
			ic[i].noisy_count = 0;
		ic[i].smin = ic[i].noisy_count - scale;
		ic[i].smax = ic[i].noisy_count + scale;
		if (ic[i].smin < 0) ic[i].smin = 0;
		if (ic[i].smax < 0) ic[i].smax = 0;
		if (ic[i].smin > fp->t) ic[i].smin = fp->t;
		if (ic[i].smax > fp->t) ic[i].smax = fp->t;
	}
}

#if PRINT_RS_TRACE || PRINT_FINAL_RULES
static void print_reservoir(struct reservoir *reservoir, size_t rs)
{
	size_t i;

	for (i = 0; i < rs; i++) {
		printf("\t%lu\t%3.2lf\t% 11.5lf\t%7d\t%7d | ", i,
				reservoir[i].c, reservoir[i].v,
				reservoir[i].sxy, reservoir[i].sx);
		print_rule(reservoir[i].r);
		printf("\n");
	}
}
#endif

static void process_rule(const struct fptree *fp,
		const size_t *AB, int ab_length, const size_t *A, int a_length,
		double eps, size_t *rs, struct reservoir *reservoir, size_t k,
		double u, size_t m, size_t sfactor, double c0)
{
	struct itemset *iA, *iAB;
	struct rule *r = NULL;
	int sup_a, sup_ab;
	double q, v, c;

	sup_a = fpt_itemset_count(fp, A, a_length);
	sup_ab = fpt_itemset_count(fp, AB, ab_length);
	iA = build_itemset(A, a_length);
	iAB = build_itemset(AB, ab_length);
	r = build_rule_A_AB(iA, iAB);
	q = quality(sup_a, sup_ab, m, sfactor, c0);

	// TODO: rule probability & cut early

	v = log(log(1/u)) - eps * q / 2;
	c = (sup_ab + 0.0) / (sup_a + 0.0);

#if PRINT_RULE_DOMAIN || PRINT_RS_TRACE
	printf("\tRule: ");
	print_rule(r);
	printf(" %d/%d (%3.2lf) %5.4lf %5.4lf %5.4lf\n",
			sup_a, sup_ab, c, q, u, v);
#endif

	if (*rs < k) {
		reservoir[*rs].r = r;
		reservoir[*rs].v = v;
		reservoir[*rs].c = c;
		reservoir[*rs].sxy = sup_ab;
		reservoir[*rs].sx = sup_a;
		*rs = *rs + 1;

		if (*rs == k) {
			qsort(reservoir, k, sizeof(reservoir[0]), reservoir_cmp);
#if PRINT_RS_TRACE
			printf("Initial reservoir:\n");
			print_reservoir(reservoir, *rs);
#endif
		}
	} else if (v < reservoir[k-1].v) {
		free_rule(reservoir[k-1].r);
		reservoir[k-1].r = r;
		reservoir[k-1].v = v;
		reservoir[k-1].c = c;
		reservoir[k-1].sxy = sup_ab;
		reservoir[k-1].sx = sup_a;
		qsort(reservoir, k, sizeof(reservoir[0]), reservoir_cmp);
#if PRINT_RS_TRACE
		printf("Inserted into reservoir, now:\n");
		print_reservoir(reservoir, *rs);
#endif
	} else
		free_rule(r);

	free_itemset(iA);
	free_itemset(iAB);
}

static void generate_and_add_all_rules(const struct fptree *fp,
		const size_t *items, size_t num_items, double eps,
		size_t *rs, struct reservoir *reservoir,
		size_t k, struct drand48_data *randbuffer,
		size_t m, size_t a_length, size_t sfactor, double c0)
{
	size_t *A = calloc(a_length, sizeof(*A));
	size_t i;
	double u;

	for (i = 0; i < a_length; i++)
		A[i] = items[i];

	drand48_r(randbuffer, &u);
	process_rule(fp, items, num_items, A, a_length, eps,
			rs, reservoir, k, u, m, sfactor, c0);

	free(A);
}

static void mine_rules_path(const struct fptree *fp,
		const struct item_count *ic,
		struct reservoir *reservoir,
		size_t *rs, size_t rlen, size_t k, double eps,
		size_t minalpha, double c0,
		size_t *items, size_t cn, size_t pos,
		struct drand48_data *randbuffer,
		double smin, double smax,
		size_t *pch, size_t pchsz)
{
	size_t *ch, i, chsz, sf, li;
	double snmin, snmax, t;

	if (!pos) {
		smin = ic[cn-1].smin;
		smax = ic[cn-1].smax;
	} else {
		li = items[pos - 1];
		snmin = snmax = 0;
		for (i = 0; i < pchsz; i++) {
			if (pch[i] == cn)
				continue;
			snmin += ic[pch[i]-1].smin;
			snmax += ic[pch[i]-1].smax;
		}

#if PRINT_PROBS
		printf("For path of length %lu:", pos);
		for (i = 0; i < pos; i++)
			printf(" %lu", items[i]);
		printf("\n\tbounds:\t[%lf, %lf]\n", smin, smax);
		printf("\tlast %lu [%lf, %lf]\n", li,
				ic[li-1].smin, ic[li-1].smax);
		printf("\titem %lu [%lf, %lf]\n", cn,
				ic[cn-1].smin, ic[cn-1].smax);
		printf("\tnone -- [%lf, %lf]\n", snmin, snmax);
#endif

		/* smin = max{0, limin - snmax, smin - snmax} */
		smin = smin - snmax;
		t = ic[li-1].smin - snmax;
		if (smin < t) smin = t;
		if (smin < 0) smin = 0;

		/* smax = min{limax, cnmax} */
		smax = ic[li-1].smax;
		if (smax > ic[cn-1].smax) smax = ic[cn-1].smax;

#if PRINT_PROBS
		printf("\tnew:\t[%lf, %lf]\n", smin, smax);
#endif

		// TODO: cut early
	}

	items[pos++] = cn;

	if (pos > rlen) {
#if PRINT_RULE_DOMAIN || PRINT_RS_TRACE
		for (i = 0; i < pos; i++)
			printf("%lu ", items[i]);
		printf("\n");
#endif

		sf = fp->has_returns ? fp->l_max_t / rlen : 1;
		generate_and_add_all_rules(fp, items, pos, eps,
				rs, reservoir, k, randbuffer, minalpha,
				rlen, sf, c0);
	}

	/* stop recursion */
	if (pos == fp->l_max_r)
		return;

	ch = fp_grph_children(fp, cn, &chsz);
	for (i = 0; i < chsz; i++)
		mine_rules_path(fp, ic, reservoir, rs, rlen, k, eps,
				minalpha, c0, items, ch[i], pos, randbuffer,
				smin, smax, ch, chsz);
	free(ch);
}

static void mine_rules_length(const struct fptree *fp,
		const struct item_count *ic,
		struct histogram *h,
		size_t rlen, size_t k, double eps,
		size_t minalpha, double c0,
		struct drand48_data *randbuffer)
{
	struct reservoir *reservoir = calloc(k, sizeof(reservoir[0]));
	size_t *items = calloc(rlen, sizeof(items[0]));
	double minc, maxc;
	size_t i, rs = 0;

	printf("\tlength %s%lu=>*: %lu rules with budget %lf each\n",
			fp->has_returns?"==":">=", rlen, k, eps);

	for (i = 1; i <= fp->n; i++) {
		mine_rules_path(fp, ic, reservoir, &rs, rlen, k, eps,
				minalpha, c0, items, i, 0, randbuffer,
				0, 0, NULL, 0);
	}

#if PRINT_FINAL_RULES
	print_reservoir(reservoir, rs);
#endif

	/* move rules from reservoir to histogram */
	minc = 1; maxc = 0;
	for (i = 0; i < rs; i++) {
		if (reservoir[i].c < minc)
			minc = reservoir[i].c;
		if (reservoir[i].c > maxc)
			maxc = reservoir[i].c;
		histogram_register(h, reservoir[i].c);
	}

	printf("Rules saved: %lu, minconf: %3.2lf, maxconf: %3.2lf\n",
			rs, minc, maxc);

	free_reservoir_array(reservoir, rs);
	free(items);
}

static struct histogram *non_private_mining(const struct fptree *fp)
{
	struct histogram *h = init_histogram();
	fpt_mine(fp, h);
	return h;
}

static void display_histograms(size_t k,
		const struct histogram *h,
		const struct histogram *nph)
{
	size_t i, n = histogram_get_count_bins();

	printf("Final private histogram:\n");
	histogram_dump(stdout, h, 1, "\t");
	printf("Final non-private histogram:\n");
	histogram_dump(stdout, nph, 1, "\t");

	printf("\n\t\t");
	for (i = 0; i < n; i++)
		printf(" %3.2lf", histogram_bin_bound(i));
	printf("\nPrecision:\t");
	for (i = 0; i < n; i++)
		printf(" %3.2lf", histogram_get_bin_full(h, i)/(k + 0.0));
	printf("\nRecall:\t\t");
	for (i = 0; i < n; i++) {
		size_t p = histogram_get_bin_full(h, i);
		size_t q = histogram_get_bin_full(nph, i);
		printf(" %3.2lf", p / (q + 0.0));
	}
	printf("\nF1:\t\t");
	for (i = 0; i < n; i++) {
		size_t p = histogram_get_bin_full(h, i);
		size_t q = histogram_get_bin_full(nph, i);
		printf(" %3.2lf", 2 * (p + 0.0) / (q + k));
	}
	printf("\n");
}

void dp2d(const struct fptree *fp, double eps, double eps_share,
		size_t k, size_t minalpha, double c0, long int seed)
{
	struct item_count *ic = calloc(fp->n, sizeof(ic[0]));
	size_t *ks = calloc(fp->l_max_r - 1, sizeof(ks[0])); /* number of rules */
	size_t *ls = calloc(fp->l_max_r - 1, sizeof(ls[0])); /* rule lengths */
	double *es = calloc(fp->l_max_r - 1, sizeof(es[0])); /* epsilons */
	struct histogram *h = init_histogram();
	double epsilon_step1 = eps * eps_share;
	struct timeval starttime, endtime;
	struct drand48_data randbuffer;
	struct histogram *nph;
	size_t i, lens;
	double t1, t2;

	printf("Running dp2D with eps=%lf, eps_share=%lf, "
			"k=%lu, minalpha=%lu c0=%lf\n",
			eps, eps_share, k, minalpha, c0);
	printf("Step 1: compute noisy counts for items with eps_1 = %lf\n",
			epsilon_step1);
	init_rng(seed, &randbuffer);
	build_items_table(fp, ic, epsilon_step1, &randbuffer);

#if PRINT_ITEM_TABLE
	printf("\n");
	for (i = 0; i < fp->n; i++)
		printf("%lu %lu %lf [%lf, %lf]\n",
				ic[i].value, ic[i].real_count,
				ic[i].noisy_count,
				ic[i].smin, ic[i].smax);
#endif

	eps = eps - epsilon_step1;
	printf("Step 2: mining %lu rules with remaining eps: %lf\n", k, eps);

	/* TODO: better split into sets, round robin for now */
	lens = fp->has_returns ? fp->l_max_r - 1 : 1;
	eps /= lens;
	for (i = 0; i < lens; i++) {
		ls[i] = fp->has_returns ? lens - i : 1;
		ks[i] = (k / lens) + ((k % lens) > i);
		es[i] = eps / ks[i];
	}

	gettimeofday(&starttime, NULL);
	for (i = 0; i < lens; i++)
		mine_rules_length(fp, ic, h, ls[i], ks[i], es[i],
				minalpha, c0, &randbuffer);
	gettimeofday(&endtime, NULL);
	t1 = starttime.tv_sec + (0.0 + starttime.tv_usec) / MICROSECONDS;
	t2 = endtime.tv_sec + (0.0 + endtime.tv_usec) / MICROSECONDS;
	printf("Total time: %5.2lf\n", t2 - t1);
	printf("%ld %ld %ld %ld\n", starttime.tv_sec, starttime.tv_usec, endtime.tv_sec, endtime.tv_usec);

	printf("Step 3: non private mining\n");
	gettimeofday(&starttime, NULL);
	nph = non_private_mining(fp);
	gettimeofday(&endtime, NULL);
	t1 = starttime.tv_sec + (0.0 + starttime.tv_usec) / MICROSECONDS;
	t2 = endtime.tv_sec + (0.0 + endtime.tv_usec) / MICROSECONDS;
	printf("Total time: %5.2lf\n", t2 - t1);
	printf("%ld %ld %ld %ld\n", starttime.tv_sec, starttime.tv_usec, endtime.tv_sec, endtime.tv_usec);

	display_histograms(k, h, nph);

	free_histogram(nph);
	free_histogram(h);
	free(ls);
	free(ks);
	free(es);
	free(ic);
}
