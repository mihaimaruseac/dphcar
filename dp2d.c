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
#define PRINT_RS_TRACE 1
#endif
/* print the returned rules */
#ifndef PRINT_FINAL_RULES
#define PRINT_FINAL_RULES 0
#endif

// TODO: m -> size_t
static double quality(int x, int y, double m, size_t sfactor)
{
	double c;

	if (x < m) x = m;
	c = (y + 0.0) / (x + 0.0);

	return c * m / sfactor;
}

struct item_count {
	size_t value;
	size_t real_count;
	double noisy_count;
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
	size_t i;

	for (i = 0; i < fp->n; i++) {
		ic[i].value = i + 1;
		ic[i].real_count = fpt_item_count(fp, i + 1);
		ic[i].noisy_count = laplace_mechanism(ic[i].real_count,
				eps, fp->l_max_t, buffer);
		/* TODO: filter some noise */
		if (ic[i].noisy_count < 0)
			ic[i].noisy_count = 0;
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
		double u, double m, size_t sfactor)
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
	q = quality(sup_a, sup_ab, m, sfactor);
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
		double m, size_t sfactor)
{
	size_t *A = calloc(num_items, sizeof(*A));
	size_t a_length;
	double u;

	for (a_length = 1; a_length < num_items; a_length++) {
		A[a_length-1] = items[a_length-1];
		drand48_r(randbuffer, &u);
		process_rule(fp, items, num_items, A, a_length, eps,
				rs, reservoir, k, u, m, sfactor);
	}

	free(A);
}

static void mine_rules_path(const struct fptree *fp,
		const struct item_count *ic,
		struct reservoir *reservoir,
		size_t *rs, size_t rlen, size_t k, double eps, double minalpha,
		size_t *items, size_t cn, size_t pos,
		struct drand48_data *randbuffer)
{
	size_t *ch, i, chsz, sf;

	items[pos++] = cn;

	if (pos == rlen || (!fp->has_returns && pos > 1)) {
#if PRINT_RULE_DOMAIN || PRINT_RS_TRACE
		for (i = 0; i < pos; i++)
			printf("%lu ", items[i]);
		printf("\n");
#endif

		sf = fp->has_returns ? fp->l_max_t / rlen : 1;
		generate_and_add_all_rules(fp, items, pos, eps,
				rs, reservoir, k, randbuffer, minalpha, sf);
	}

	/* stop recursion */
	if (pos == rlen)
		return;

	// TODO: check probability and cut early

	ch = fp_grph_children(fp, cn, &chsz);
	for (i = 0; i < chsz; i++)
		mine_rules_path(fp, ic, reservoir, rs, rlen, k, eps, minalpha,
				items, ch[i], pos, randbuffer);
	free(ch);
}

static void mine_rules_length(const struct fptree *fp,
		const struct item_count *ic,
		struct histogram *h,
		size_t rlen, size_t k, double eps, double minalpha,
		struct drand48_data *randbuffer)
{
	struct reservoir *reservoir = calloc(k, sizeof(reservoir[0]));
	size_t *items = calloc(rlen, sizeof(items[0]));
	double minc, maxc;
	size_t i, rs = 0;

	printf("\tlength %s%lu: %lu rules with budget %lf each\n",
			fp->has_returns?"==":"<=", rlen, k, eps);

	for (i = 1; i <= fp->n; i++) {
		mine_rules_path(fp, ic, reservoir, &rs, rlen, k, eps, minalpha,
				items, i, 0, randbuffer);
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

void dp2d(const struct fptree *fp, double eps, double eps_share,
		size_t k, double minalpha, long int seed)
{
	struct item_count *ic = calloc(fp->n, sizeof(ic[0]));
	size_t *ks = calloc(fp->l_max_r - 1, sizeof(ks[0])); /* number of rules */
	size_t *ls = calloc(fp->l_max_r - 1, sizeof(ls[0])); /* rule lengths */
	double *es = calloc(fp->l_max_r - 1, sizeof(es[0])); /* epsilons */
	double epsilon_step1 = eps * eps_share;
	struct histogram *h = init_histogram();
	struct timeval starttime, endtime;
	struct drand48_data randbuffer;
	size_t i, lens;
	double t1, t2;

	printf("Running dp2D with eps=%lf, eps_share=%lf, "
			"k=%lu, minalpha=%lf\n",
			eps, eps_share, k, minalpha);
	printf("Step 1: compute noisy counts for items with eps_1 = %lf\n",
			epsilon_step1);
	init_rng(seed, &randbuffer);
	build_items_table(fp, ic, epsilon_step1, &randbuffer);

#if PRINT_ITEM_TABLE
	printf("\n");
	for (i = 0; i < fp->n; i++)
		printf("%lu %lu %lf\n", ic[i].value, ic[i].real_count, ic[i].noisy_count);
#endif

	eps = eps - epsilon_step1;
	printf("Step 2: mining %lu rules with remaining eps: %lf\n", k, eps);

	// TODO: better split into sets, round robin for now
	lens = fp->has_returns ? fp->l_max_r - 1 : 1;
	eps /= lens;
	for (i = 0; i < lens; i++) {
		ls[i] = fp->has_returns ? i + 2 : fp->l_max_r;
		ks[i] = (k / lens) + ((k % lens) > i);
		es[i] = eps / ks[i];
	}

	gettimeofday(&starttime, NULL);
	for (i = 0; i < lens; i++)
		mine_rules_length(fp, ic, h, ls[i], ks[i], es[i], minalpha, &randbuffer);

	gettimeofday(&endtime, NULL);
	t1 = starttime.tv_sec + (0.0 + starttime.tv_usec) / MICROSECONDS;
	t2 = endtime.tv_sec + (0.0 + endtime.tv_usec) / MICROSECONDS;

	printf("Total time: %5.2lf\n", t2 - t1);
	printf("%ld %ld %ld %ld\n", starttime.tv_sec, starttime.tv_usec, endtime.tv_sec, endtime.tv_usec);

	printf("Final histogram:\n");
	histogram_dump(stdout, h, 1, "\t");

	free_histogram(h);
	free(ls);
	free(ks);
	free(es);
	free(ic);
}
