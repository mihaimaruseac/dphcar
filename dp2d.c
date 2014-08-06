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
/* print the items considered for each rule generation step */
#ifndef PRINT_ITEM_DOMAIN
#define PRINT_ITEM_DOMAIN 0
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
#define PRINT_FINAL_RULES 1
#endif

static double quality(int x, int y)
{
#if 0
	return y - x;
#else
	if (x < 1000) x = 1000;
	return (y + 0.0) / (0.0 +  x);
#endif
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

static void generate_and_add_all_rules(const struct fptree *fp,
		const int *items, size_t num_items, double eps,
		size_t *rs, struct reservoir *reservoir, size_t k,
		struct drand48_data *randbuffer)
{
#define RULE_A 1
#define RULE_B 2
#define RULE_END 3

	size_t i, a_length, b_length, ab_length;
	struct itemset *iA, *iAB;
	struct rule *r = NULL;
	unsigned char *ABi;
	double q, u, v, c;
	int sup_a, sup_ab;
	int *A, *B, *AB;

	ABi = calloc(num_items, sizeof(ABi[0]));
	AB = calloc(num_items, sizeof(AB[0]));
	A = calloc(num_items, sizeof(A[0]));
	B = calloc(num_items, sizeof(B[0]));

	ABi[num_items - 1] = 1;
	ABi[num_items - 2] = 1;

	/* O(3^num_items) !!!! */
	while (1) {
		ABi[num_items - 1]++;
		a_length = b_length = ab_length = 0;

		for (i = num_items - 1; i > 0; i--)
			if (ABi[i] == RULE_END) {
				ABi[i] = 0;
				ABi[i-1]++;
		}

		if (ABi[0] == RULE_END) break;

		for (i = 0; i < num_items; i++)
			if (ABi[i]) {
				AB[ab_length++] = items[i];
				switch (ABi[i]) {
				case RULE_B: B[b_length++] = items[i]; break;
				case RULE_A: A[a_length++] = items[i]; break;
				}
			}

		if (a_length == 0) continue;
		if (b_length == 0) continue;

		/* build new rule and compute stats */
		sup_a = fpt_itemset_count(fp, A, a_length);
		sup_ab = fpt_itemset_count(fp, AB, ab_length);
		iA = build_itemset(A, a_length);
		iAB = build_itemset(AB, ab_length);
		r = build_rule_A_AB(iA, iAB);
		q = quality(sup_a, sup_ab);
		drand48_r(randbuffer, &u);
		v = eps * q / 2 + log(log(1/u));
		c = (sup_ab + 0.0) / (sup_a + 0.0);

#if PRINT_RULE_DOMAIN || PRINT_RS_TRACE
		printf("\tRule: ");
		print_rule(r);
		printf(" %d/%d (%3.2lf) %5.2lf %5.2lf %5.2lf\n", 
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

	free(ABi);
	free(AB);
	free(A);
	free(B);

#undef RULE_A
#undef RULE_B
#undef RULE_END
}

static size_t update_fm(size_t fm, size_t fM, int mis, double mu,
		size_t n, int minth, struct item_count *ic)
{
	while (fm < n /* don't exit the buffer */
		&& ic[fm].noisy_count >= minth /* don't pass below threshold */
		&& fm < fM + mis /* number of items condition */
		&& ic[fm].noisy_count >= mu * ic[fM].noisy_count /* spread */)
		fm++;
	return fm;
}

void dp2d(const struct fptree *fp, double eps, double eps_share, int minth,
		double mu, size_t mis, size_t k)
{
	struct reservoir *reservoir = calloc(k, sizeof(reservoir[0]));
	struct item_count *ic = calloc(fp->n, sizeof(ic[0]));
	int *items = calloc(mis + 1, sizeof(items[0]));
	double epsilon_step1 = eps * eps_share;
	struct drand48_data randbuffer;
	size_t i, fm, fM, rs;
	double maxc, minc;

	init_rng(&randbuffer);

	printf("Running dp2D with minth=%d, eps=%lf, eps_share=%lf, mu=%lf, "
			"mis=%lu, k=%lu\n", minth, eps, eps_share, mu, mis, k);

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
	fm = fM = 0;
	fm = update_fm(fm, fM, mis, mu, fp->n, minth, ic);
	while (fm != fM) {
#if PRINT_ITEM_DOMAIN || PRINT_RULE_DOMAIN || PRINT_RS_TRACE
		printf("Domain: %lu-%lu:\n", fm, fM);
#endif

		/* extract items from which to generate rules */
		for (i = fM; i <= fm; i++) {
#if PRINT_ITEM_DOMAIN
			printf("\t%d %d %lf\n", ic[i].value,
					ic[i].real_count, ic[i].noisy_count);
#endif
			items[i - fM] = ic[i].value;
		}

#if PRINT_ITEM_DOMAIN || PRINT_RULE_DOMAIN || PRINT_RS_TRACE
		printf("\n");
#endif

		generate_and_add_all_rules(fp, items, 1 + fm - fM, eps,
				&rs, reservoir, k, &randbuffer);

		fM++;
		fm = update_fm(fm, fM, mis, mu, fp->n, minth, ic);
	}

	/* TODO: compute min c inside the array of results (simulate outputing rules) */
#if PRINT_FINAL_RULES
	print_reservoir(reservoir, k);
#endif

	minc = 1; maxc = 0;
	for (i = 0; i < k; i++) {
		if (reservoir[i].c < minc)
			minc = reservoir[i].c;
		if (reservoir[i].c > maxc)
			maxc = reservoir[i].c;
	}
	printf("minconf: %3.2lf, maxconf: %3.2lf\n", minc, maxc);

	free_reservoir_array(reservoir, k);
	free(items);
	free(ic);
}
