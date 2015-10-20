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
/* print every step in computing probabilities */
#ifndef PRINT_PROBS_TRACE
#define PRINT_PROBS_TRACE 0
#endif

static double quality(size_t x, size_t y, size_t sfactor, double c0)
{
	double c;

	c = y - c0 * x;
	if (c < 0) c = 0;

	return c / sfactor;
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

#if PRINT_PROBS_TRACE
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

static void compute_seq_bound_at(const struct fptree *fp,
		const struct item_count *ic,
		const size_t *items, size_t num_items, size_t cur_item,
		double *pmin, double *pmax)
{
	size_t cn, li, i, *ch, chsz;
	double snmin, snmax, t;

	if (cur_item == num_items)
		return;

	if (!cur_item) {
		*pmin = ic[items[0]-1].smin;
		*pmax = ic[items[0]-1].smax;
		goto next;
	}

	cn = items[cur_item]; /* current item */
	li = items[cur_item - 1];
	ch = fp_grph_children(fp, li, &chsz);
	snmin = snmax = 0;
	for (i = 0; i < chsz; i++) {
		if (ch[i] == cn)
			continue;
		snmin += ic[ch[i]-1].smin;
		snmax += ic[ch[i]-1].smax;
	}
	free(ch);

#if PRINT_PROBS_TRACE
	printf("At item %lu:\n", cn);
	printf("\tbounds:\t[%lf, %lf]\n", *pmin, *pmax);
	printf("\tlast %lu [%lf, %lf]\n", li,
			ic[li-1].smin, ic[li-1].smax);
	printf("\titem %lu [%lf, %lf]\n", cn,
			ic[cn-1].smin, ic[cn-1].smax);
	printf("\trest -- [%lf, %lf]\n", snmin, snmax);
#endif

	/* smin = max{0, limin - snmax, smin - snmax} */
	*pmin = *pmin - snmax;
	t = ic[li-1].smin - snmax;
	if (*pmin < t) *pmin = t;
	if (*pmin < 0) *pmin = 0;

	/* smax = min{limax, cnmax} */
	*pmax = ic[li-1].smax;
	if (*pmax > ic[cn-1].smax) *pmax = ic[cn-1].smax;

#if PRINT_PROBS_TRACE
	printf("\tnew:\t[%lf, %lf]\n", *pmin, *pmax);
#endif

next:
	compute_seq_bound_at(fp, ic, items, num_items, cur_item + 1, pmin, pmax);
}

static void compute_seq_bounds(const struct fptree *fp,
		const struct item_count *ic,
		const size_t *items, size_t num_items,
		double *pmin, double *pmax)
{
	size_t i;

#if PRINT_PROBS || PRINT_PROBS_TRACE
	printf("Computing bounds on seq of length %lu:", num_items);
	for (i = 0; i < num_items; i++)
		printf(" %lu", items[i]);
	printf("\n");
#endif

	*pmin = *pmax = 0;
	compute_seq_bound_at(fp, ic, items, num_items, 0, pmin, pmax);

	/* pmax cannot be higher than min{max{sup}} */
	for (i = 0; i < num_items; i++)
		if (*pmax > ic[items[i]-1].smax)
			*pmax = ic[items[i]-1].smax;

	/* ensure valid interval */
	if (*pmin > *pmax) *pmin = *pmax;

#if PRINT_PROBS || PRINT_PROBS_TRACE
	printf("\t\t seq bounds [%lf, %lf]\n", *pmin, *pmax);
#endif
}

static void compute_rule_bounds(const struct fptree *fp,
		const struct item_count *ic,
		const size_t *AB, size_t ab_length, size_t a_length,
		double c0, double *pmin, double *pmax)
{
	double b_length = ab_length - a_length;
	double minsupAB, maxsupAB;
	double minsupA, maxsupA;
	double minsupB, maxsupB;
	double t;

#if PRINT_PROBS || PRINT_PROBS_TRACE
	size_t i;

	printf("Computing bounds on rule: ");
	for (i = 0; i < a_length; i++)
		printf(" %lu", AB[i]);
	printf(" => ");
	for (; i < ab_length; i++)
		printf(" %lu", AB[i]);
	printf("\n");
#endif

	compute_seq_bounds(fp, ic, AB,           ab_length, &minsupAB, &maxsupAB);
	compute_seq_bounds(fp, ic, AB,            a_length,  &minsupA,  &maxsupA);
	compute_seq_bounds(fp, ic, AB + a_length, b_length,  &minsupB,  &maxsupB);

	/* smin = max{0, minsupAB - c0 * maxsupA, minsupB - c0 * maxsupA} */
	*pmin = minsupAB;
	if (*pmin < minsupB) *pmin = minsupB;
	*pmin -= c0 * maxsupA;
	if (*pmin < 0) *pmin = 0;

	/* smax = min{(1 - c0) * maxsupAB, (1- c0) * maxsupA, maxsupAB - c0 * minsupA} */
	*pmax = maxsupAB;
	if (*pmax > maxsupA) *pmax = maxsupA;
	*pmax *= 1 - c0;
	t = maxsupAB - c0 * minsupA;
	if (*pmax > t) *pmax = t;

	/* ensure valid interval */
	if (*pmin < 0) *pmin = 0;
	if (*pmax < 0) *pmax = 0;
	if (*pmin > *pmax) *pmin = *pmax;

#if PRINT_PROBS || PRINT_PROBS_TRACE
	printf("\t\trule bounds [%lf, %lf]\n", *pmin, *pmax);
#endif
}

static void process_rule(const struct fptree *fp,
		const struct item_count *ic,
		const size_t *AB, int ab_length, const size_t *A, int a_length,
		double eps, size_t *rs, struct reservoir *reservoir, size_t k,
		double u, size_t sfactor, double c0, double cmin)
{
	struct itemset *iA, *iAB;
	struct rule *r = NULL;
	int sup_a, sup_ab;
	double pmin, pmax;
	double q, v, c;

	sup_a = fpt_itemset_count(fp, A, a_length);
	sup_ab = fpt_itemset_count(fp, AB, ab_length);
	iA = build_itemset(A, a_length);
	iAB = build_itemset(AB, ab_length);
	r = build_rule_A_AB(iA, iAB);
	q = quality(sup_a, sup_ab, sfactor, c0);

	compute_rule_bounds(fp, ic, AB, ab_length, a_length, c0, &pmin, &pmax);
	if (pmax < cmin) {
#if PRINT_PROBS
		printf("Rule cutoff\n");
#endif
		goto end;
	}

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
	} else {
end:
		free_rule(r);
	}

	free_itemset(iA);
	free_itemset(iAB);
}

static void generate_and_add_all_rules(const struct fptree *fp,
		const struct item_count *ic,
		const size_t *items, size_t num_items, double eps,
		size_t *rs, struct reservoir *reservoir,
		size_t k, struct drand48_data *randbuffer,
		size_t a_length, size_t sfactor, double c0, double cmin)
{
	size_t *A = calloc(a_length, sizeof(*A));
	size_t i;
	double u;

	for (i = 0; i < a_length; i++)
		A[i] = items[i];

	drand48_r(randbuffer, &u);
	process_rule(fp, ic, items, num_items, A, a_length, eps,
			rs, reservoir, k, u, sfactor, c0, cmin);

	free(A);
}

static void mine_rules_path(const struct fptree *fp,
		const struct item_count *ic,
		struct reservoir *reservoir,
		size_t *rs, size_t rlen, size_t k, double eps,
		double c0, double sigma_min, double cmin, double fr,
		size_t *items, size_t cn, size_t pos,
		struct drand48_data *randbuffer)
{
	size_t *ch, i, chsz, sf;
	double smin, smax;

	items[pos++] = cn;
	compute_seq_bounds(fp, ic, items, pos, &smin, &smax);
	if (smax < sigma_min) {
#if PRINT_PROBS
		printf("Cutoff sequence\n");
#endif
		return;
	}

	if (pos > rlen) {
#if PRINT_RULE_DOMAIN || PRINT_RS_TRACE
		for (i = 0; i < pos; i++)
			printf("%lu ", items[i]);
		printf("\n");
#endif

		sf = fp->has_returns ? fp->l_max_t / rlen : 1;
		generate_and_add_all_rules(fp, ic, items, pos, eps,
				rs, reservoir, k, randbuffer,
				rlen, sf, c0, cmin);
	}

	/* stop recursion */
	if (pos == fp->l_max_r)
		return;

	ch = fp_grph_children(fp, cn, &chsz);
	// TODO: get only fr children, sorted/permuted
	for (i = 0; i < chsz; i++)
		mine_rules_path(fp, ic, reservoir, rs, rlen, k, eps,
				c0, sigma_min, cmin, fr,
				items, ch[i], pos, randbuffer);
	free(ch);
}

static void mine_rules_length(const struct fptree *fp,
		const struct item_count *ic,
		struct histogram *h,
		size_t rlen, size_t k, double eps,
		double c0, double sigma_min, double cmin, double fr,
		struct drand48_data *randbuffer)
{
	struct reservoir *reservoir = calloc(k, sizeof(reservoir[0]));
	size_t *items = calloc(fp->l_max_r, sizeof(items[0]));
	double minc, maxc;
	size_t i, rs = 0;

	printf("\tlength %s%lu=>*: %lu rules with budget %lf each and "
			"fillratio %lf\n", fp->has_returns?"==":">=",
			rlen, k, eps, fr);

	// TODO: get only fr children, sorted/permuted
	for (i = 1; i <= fp->n; i++) {
		mine_rules_path(fp, ic, reservoir, &rs, rlen, k, eps,
				c0, sigma_min, cmin, fr, items, i, 0,
				randbuffer);
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
		size_t k, double c0, double sigma_min, double cmin,
		long int seed)
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
	double t1, t2, fr;
	size_t i, lens;

	printf("Running dp2D with eps=%lf, eps_share=%lf, "
			"k=%lu, c0=%lf\n",
			eps, eps_share, k, c0);
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
	fr = 1; // TODO: formula
	printf("Step 2: mining %lu rules with remaining eps: %lf\n", k, eps);

	/* TODO: better split into sets, round robin for now */
	lens = fp->has_returns ? fp->l_max_r - 1 : 1;
	// no need to eps /= lens; because parallel comp
	for (i = 0; i < lens; i++) {
		ls[i] = fp->has_returns ? lens - i : 1;
		ks[i] = (k / lens) + ((k % lens) > i);
		es[i] = eps / ks[i];
	}

	gettimeofday(&starttime, NULL);
	for (i = 0; i < lens; i++)
		if (ks[i])
			mine_rules_length(fp, ic, h, ls[i], ks[i], es[i],
					c0, sigma_min, cmin, fr,
					&randbuffer);
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
