#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>
#include <mpfr.h>

#include "dp2d.h"
#include "fp.h"
#include "globals.h"
#include "histogram.h"
#include "rule.h"

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

static double quality(int x, int y)
{
	double c;

	if (x < 1000) x = 1000;
	c = (y + 0.0) / (x + 0.0);

	return 1000 * pow(c, 1);
}

struct item_count {
	int value;
	int real_count;
	double noisy_count;
};

struct reservoir {
	struct rule *r;
	double c;
	double v;
};

static void free_reservoir_array(struct reservoir *reservoir, int size)
{
	int i;

	for (i = 0 ; i < size; i++)
		free_rule(reservoir[i].r);
	free(reservoir);
}

static int ic_noisy_cmp(const void *a, const void *b)
{
	const struct item_count *ia = a, *ib = b;
	return double_cmp_r(&ia->noisy_count, &ib->noisy_count);
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
		ic[i].real_count = fpt_item_count(fp, i);
		ic[i].noisy_count = laplace_mechanism(ic[i].real_count,
				eps, 1, buffer);
		/* TODO: filter some noise */
		if (ic[i].noisy_count < 0)
			ic[i].noisy_count = 0;
	}

	qsort(ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);
}

#if PRINT_RS_TRACE || PRINT_FINAL_RULES
static void print_reservoir(struct reservoir *reservoir, size_t rs)
{
	size_t i;

	for (i = 0; i < rs; i++) {
		printf("\t%lu\t%3.2lf\t% 11.5lf | ", i,
				reservoir[i].c, reservoir[i].v);
		print_rule(reservoir[i].r);
		printf("\n");
	}
}
#endif

static void process_rule(const struct fptree *fp,
		const int *AB, int ab_length, const int *A, int a_length,
		double eps, size_t *rs, struct reservoir *reservoir, size_t k,
		double u)
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
	q = quality(sup_a, sup_ab);
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
		const int *items, size_t num_items, size_t st, double eps,
		size_t *rs, size_t *rs2,
		struct reservoir *reservoir, struct reservoir *reservoir2,
		size_t k, struct drand48_data *randbuffer)
{
	size_t i, j, l, max=1<<num_items, a_length, ab_length, max2;
	int *A, *AB;
	double u;

	k /= 2;
	A = calloc(num_items, sizeof(A[0]));
	AB = calloc(num_items, sizeof(AB[0]));

	/* generate rule's AB */
	for (i = st; i < max; i++) {
		ab_length = 0;
		for (j = 0; j < num_items; j++)
			if (i & (1 << j))
				AB[ab_length++] = items[j];
		if (ab_length < 2) continue;

		max2 = (1 << ab_length) - 1;
		for (j = 1; j < max2; j++) {
			a_length = 0;
			for (l = 0; l < ab_length; l++)
				if (j & (1 << l))
					A[a_length++] = AB[l];
			drand48_r(randbuffer, &u);
			process_rule(fp, AB, ab_length, A, a_length,
					eps, rs, reservoir, k, u);
			process_rule(fp, AB, ab_length, A, a_length,
					eps, rs2, reservoir2, k, 1 - u);
		}
	}

	free(A);
	free(AB);
}

void dp2d(const struct fptree *fp, double eps, double eps_share, int minth,
		size_t mis, size_t k, long int seed)
{
	struct reservoir *reservoir = calloc(k/2, sizeof(reservoir[0]));
	struct reservoir *reservoir2 = calloc(k/2, sizeof(reservoir[0]));
	struct item_count *ic = calloc(fp->n, sizeof(ic[0]));
	int *items = calloc(mis + 1, sizeof(items[0]));
	struct rule_table *rt = init_rule_table();
	struct histogram *h = init_histogram();
	double epsilon_step1 = eps * eps_share;
	struct drand48_data randbuffer;
	size_t i, fm, rs, rs2, st;
	double maxc, minc;

	init_rng(seed, &randbuffer);

	printf("Running dp2D with minth=%d, eps=%lf, eps_share=%lf, "
			"mis=%lu, k=%lu\n", minth, eps, eps_share, mis, k);

	printf("Step 1: compute noisy counts for items with eps_1 = %lf\n",
			epsilon_step1);
	build_items_table(fp, ic, epsilon_step1, &randbuffer);

#if PRINT_ITEM_TABLE
	printf("\n");
	for (i = 0; i < fp->n; i++)
		printf("%d %d %lf\n", ic[i].value, ic[i].real_count, ic[i].noisy_count);
#endif

	eps = (eps - epsilon_step1) / k;
	printf("Step 2: mining with remaining eps (per rule): %lf\n", eps);

	/* select mining domains */
	rs = 0; /* empty reservoir */
	rs2 = 0; /* empty reservoir */
	st = 3;

	/* initial items */
	for (fm = 0; fm < mis; fm++)
		items[fm] = ic[fm].value;

	while (fm < fp->n) {
#if PRINT_RULE_DOMAIN || PRINT_RS_TRACE
		printf("Domain: %lu: ", fm);
		for (i = 0; i < mis; i++)
			printf("%d ", items[i]);
		printf("\n");
#endif

		generate_and_add_all_rules(fp, items, mis, st, eps,
				&rs, &rs2, reservoir, reservoir2,
				k, &randbuffer);
		st = (1 << (mis - 1)) + 1;

		for (i = 0; i < mis - 1; i++)
			items[i] = items[i+1];
		items[mis - 1] = ic[fm++].value;
		if (ic[fm-1].noisy_count < minth)
			break;
	}

#if PRINT_FINAL_RULES
	print_reservoir(reservoir, rs);
#endif

	minc = 1; maxc = 0;
	for (i = 0; i < rs; i++) {
		if (reservoir[i].c < minc)
			minc = reservoir[i].c;
		if (reservoir[i].c > maxc)
			maxc = reservoir[i].c;
		save_rule2(rt, reservoir[i].r, reservoir[i].c);
	}
	for (i = 0; i < rs2; i++) {
		if (reservoir2[i].c < minc)
			minc = reservoir2[i].c;
		if (reservoir2[i].c > maxc)
			maxc = reservoir2[i].c;
		save_rule2(rt, reservoir2[i].r, reservoir2[i].c);
		//histogram_register(h, reservoir[i].c);
	}
	for (i = 0; i < rt->sz; i++)
		histogram_register(h, rt->c[i]);
	printf("Rules saved: %lu/%lu, minconf: %3.2lf, maxconf: %3.2lf\n",
			rt->sz, 2*rs, minc, maxc);

	printf("Final histogram:\n");
	histogram_dump(h, 1, "\t");

	free_reservoir_array(reservoir, k/2);
	//free_reservoir_array(reservoir2, k); don't free this since we might have shared rules :(
	free_rule_table2(rt);
	free_histogram(h);
	free(items);
	free(ic);
}
