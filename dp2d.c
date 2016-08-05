#include <float.h>
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

static double quality(int x, int y, double c0)
{
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

static void generate_rules(const size_t *items, size_t lmax,
		const struct fptree *fp, const struct item_count *ic,
		double *minc, double *maxc, struct histogram *h,
		int *seen, size_t *seenix)
{
	size_t i, j, l, max=1<<lmax, max2, a_length, ab_length;
	int sup_ab, sup_a;
	int *A, *AB;
	double c;

	A = calloc(lmax, sizeof(A[0]));
	AB = calloc(lmax, sizeof(AB[0]));

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

		max2 = (1 << ab_length) - 1;
		for (j = 1; j < max2; j++) {
			a_length = 0;
			for (l = 0; l < ab_length; l++)
				if (j & (1 << l))
					A[a_length++] = AB[l];

			sup_ab = fpt_itemset_count(fp, AB, ab_length);
			sup_a = fpt_itemset_count(fp, A, a_length);
			c = (sup_ab + 0.0) / sup_a;
			if (c < *minc) *minc = c;
			if (c > *maxc) *maxc = c;
			histogram_register(h, c);

#if PRINT_FINAL_RULES
			print_this_rule(A, AB, a_length, ab_length, c);
#endif
		}
	}

	free(A);
	free(AB);
}

/**
 * Analyze the current items to see if we can select a good rule lattice.
 */
static void analyze_items(const size_t *items, size_t lmax,
		double *bv, size_t *bitems,
		const struct fptree *fp, const struct item_count *ic,
		double c0, double eps, struct drand48_data *randbuffer)
{
	int *AB = calloc(lmax, sizeof(AB[0]));
	int sup_ab, sup_a;
	double q, u, v;
	size_t i, j, k;

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
		q = quality(sup_a, sup_ab, c0);
		drand48_r(randbuffer, &u);
		v = log(log(1/u)) - eps * q / 2;

#if PRINT_RULE_LATTICE
		printf("\t%3d -> {}: c=%7.6f q=%5.2f u=%5.2f v=%5.2f\n",
				AB[i], (sup_ab + 0.0)/sup_a, q, u, v);
#endif

		if (v < *bv) {
			*bv = v;
			for (j = 0; j < lmax; j++)
				bitems[j] = items[j];
			k = bitems[0];
			bitems[0] = bitems[i];
			bitems[i] = k;

#if PRINT_RULE_LATTICE_TRACE
			printf("Current best items: %lu(%d) -> ",
					bitems[0], ic[bitems[0]].value);
			for (j = 1; j < lmax; j++)
				printf("%lu(%d) ", bitems[j],
						ic[bitems[j]].value);
			printf(": c=%7.6f q=%5.2f u=%5.2f v=%5.2f\n",
					(sup_ab + 0.0)/sup_a, q, u, v);
#endif
		}
	}

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
static int update_items(size_t *items, size_t lmax, size_t n,
		const size_t *seen, size_t seenlen)
{
	size_t ix = lmax - 1, ok;

	do {
		do {
			/* try next */
			ok = 1;
			items[ix]++;

			/* if impossible, move to one below */
			if (items[ix] == n) {
				if (!ix) return 1; /* end of generation */
				ix--;
				ok = 0; /* this level was not finished */
				break; /* restart process for ix - 1 */
			}

			/* check to not generate seen set */
			if (already_seen(items, lmax, seen, seenlen))
				ok = 0; /* get next set */
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
static inline void init_items(size_t *items, size_t lmax, size_t n,
		const size_t *seen, size_t seenlen)
{
	size_t i;

	for (i = 0; i < lmax; i++)
		items[i] = i;

	if (already_seen(items, lmax, seen, seenlen))
		update_items(items, lmax, n, seen, seenlen);
}

void dp2d(const struct fptree *fp, double eps, double eps_ratio1,
		double c0, size_t lmax, size_t k, long int seed)
{
	struct item_count *ic = calloc(fp->n, sizeof(ic[0]));
	size_t i, j, end, seenix, seenitsix, numits;
	double epsilon_step1 = eps * eps_ratio1;
	struct histogram *h = init_histogram();
	struct drand48_data randbuffer;
	size_t *items, *bitems, *seen;
	double bv, minc, maxc;
	int *seenits;

	init_rng(seed, &randbuffer);
	items = calloc(lmax, sizeof(items[0]));

	printf("Running dp2D with eps=%lf, eps_step1=%lf, k=%lu, c0=%5.2lf, "
			"rmax=%lu\n", eps, epsilon_step1, k, c0, lmax);

	printf("Step 1: compute noisy counts for items with eps_1 = %lf\n",
			epsilon_step1);
	build_items_table(fp, ic, epsilon_step1, &randbuffer);

#if PRINT_ITEM_TABLE
	printf("\n");
	for (i = 0; i < fp->n; i++)
		printf("%d %d %lf\n", ic[i].value, ic[i].real_count, ic[i].noisy_count);
#endif

	eps = eps - epsilon_step1;
	eps /= k;
	seen = calloc(k * lmax, sizeof(seen[0]));
	seenits = calloc(k * (lmax+ 1) * (1 << lmax), sizeof(seenits[0]));
	printf("Step 2: mining %lu steps each with eps %lf\n", k, eps);

	struct timeval starttime;
	gettimeofday(&starttime, NULL);

	seenix = 0;
	seenitsix = 0;
	/* TODO: should be fp->n or a constant that is determined by code */
	numits = 20;
	minc = 1;
	maxc = 0;
	for (i = 0; i < k; i++) {
		init_items(items, lmax, numits, seen, seenix);
		bitems = calloc(lmax, sizeof(bitems[0]));
		bv = DBL_MAX;
		do {
			analyze_items(items, lmax, &bv, bitems, fp, ic, c0,
					eps, &randbuffer);
			end = update_items(items, lmax, numits, seen, seenix);
		} while (!end);

#if PRINT_RULE_LATTICE_TRACE || PRINT_FINAL_RULES
		printf("Selected items: %lu(%d) -> ",
				bitems[0], ic[bitems[0]].value);
		for (j = 1; j < lmax; j++)
			printf("%lu(%d) ", bitems[j],
					ic[bitems[j]].value);
		printf("\n");
#endif

		generate_rules(bitems, lmax, fp, ic, &minc, &maxc, h,
				seenits, &seenitsix);

		/* remove bitems from future items */
		qsort(bitems, lmax, sizeof(bitems[0]), int_cmp);
		for (j = 0; j < lmax; j++)
			seen[seenix++] = bitems[j];

		free(bitems);
	}

	printf("Rules saved: %lu, minconf: %3.2lf, maxconf: %3.2lf\n",
			histogram_get_all(h), minc, maxc);

	struct timeval endtime;
	gettimeofday(&endtime, NULL);
	double t1 = starttime.tv_sec + (0.0 + starttime.tv_usec) / MICROSECONDS;
	double t2 = endtime.tv_sec + (0.0 + endtime.tv_usec) / MICROSECONDS;

	printf("Total time: %5.2lf\n", t2 - t1);
	printf("%ld %ld %ld %ld\n", starttime.tv_sec, starttime.tv_usec, endtime.tv_sec, endtime.tv_usec);

	printf("Final histogram:\n");
	histogram_dump(stdout, h, 1, "\t");

	free_histogram(h);
	free(seenits);
	free(items);
	free(seen);
	free(ic);
}

void dp2d_np(const struct fptree *fp, double c0, size_t lmax, long int seed)
{
}
