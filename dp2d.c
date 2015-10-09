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

static double quality(int x, int y, double m)
{
	double c;

	if (x < m) x = m;
	c = (y + 0.0) / (x + 0.0);

	return m * pow(c, 1);
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
		ic[i].real_count = fpt_item_count(fp, i + 1);
		ic[i].noisy_count = laplace_mechanism(ic[i].real_count,
				eps, fp->l_max_t, buffer);
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
		double u, double m)
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
	q = quality(sup_a, sup_ab, m);
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
		size_t k, struct drand48_data *randbuffer, double m)
{
	size_t i, j, a_length, ab_length;
	size_t *A, *AB, **next, *nsz;
	int kp, *ptr;
	double u;

	A = calloc(num_items, sizeof(A[0]));
	AB = calloc(num_items, sizeof(AB[0]));

	next = calloc(num_items, sizeof(next[0]));
	ptr = calloc(num_items, sizeof(ptr[0]));
	nsz = calloc(num_items, sizeof(nsz[0]));

	next[0] = calloc(num_items, sizeof(next[0][0]));
	for (i = 0; i < num_items; i++)
		next[0][i] = items[i];

	kp = 0; ptr[0] = -1; nsz[0] = num_items;
	while (kp >= 0) {
		ptr[kp]++;

		/* backtrack */
		if (ptr[kp] == (int)nsz[kp]) {
			free(next[kp]);
			next[kp] = 0;
			kp--;
			continue;
		}

		/* check that current item is in current set of items */
		for (j = 0; j < num_items; j++)
			if (items[j] == next[kp][ptr[kp]])
				break;
		if (j == num_items)
			continue;

		/* combination generated */
		if (kp) { /* need at least 2 items */
			/* fill AB */
			for (j = 0; j <= (size_t)kp; j++)
				AB[j] = next[j][ptr[j]];
			ab_length = kp + 1;

			/* fill A, generate & process rule */
			for (a_length = 1; a_length < ab_length; a_length++) {
				A[a_length-1] = AB[a_length-1];
				drand48_r(randbuffer, &u);
				process_rule(fp, AB, ab_length, A, a_length,
						eps, rs, reservoir, k, u, m);
			}
		}

		if (kp == (int)(num_items-1)) /* max length */
			continue;

		/* next */
		kp++;
		ptr[kp] = -1;
		next[kp] = fp_grph_children(fp, next[kp-1][ptr[kp-1]], &nsz[kp]);
	}

	free(next);
	free(ptr);
	free(nsz);
	free(AB);
	free(A);
}

void dp2d(const struct fptree *fp, double eps, double eps_share,
		size_t mis, size_t nt, size_t k,
		double minalpha, long int seed)
{
	struct item_count *ic = calloc(fp->n, sizeof(ic[0]));
	size_t *ls = calloc(fp->l_max_r - 1, sizeof(ls[0])); /* rule lengths */
	size_t *ks = calloc(fp->l_max_r - 1, sizeof(ks[0])); /* number of rules */
	double *es = calloc(fp->l_max_r - 1, sizeof(es[0])); /* epsilons */
	double epsilon_step1 = eps * eps_share;
	size_t i, lens;
#if 0 /* moving to graphs */
	size_t i, j, fm = 0, rs, rsi, ct, csz, tmp, tmp2;
#else
	(void)nt; (void)k;
#endif
	struct histogram *h = init_histogram();
	struct drand48_data randbuffer;
#if 0 /* moving to graphs */
	size_t *items, *ch;
	double maxc, minc;
#endif

	init_rng(seed, &randbuffer);
#if 0 /* moving to graphs */
	items = calloc(mis + 1, sizeof(items[0]));
#endif

	printf("Running dp2D with eps=%lf, eps_share=%lf, "
			"mis=%lu, k=%lu, minalpha=%lf\n",
			eps, eps_share, mis, k, minalpha);

	printf("Step 1: compute noisy counts for items with eps_1 = %lf\n",
			epsilon_step1);
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
		printf("\tlength %s%lu: %lu rules with budget %lf each\n",
				fp->has_returns?"==":"<=",
				ls[i], ks[i], es[i]);
	}

	struct timeval starttime;
	gettimeofday(&starttime, NULL);

#if 0 /* moving to graphs */
	/* select mining domains */
	struct reservoir *reservoir = calloc(k, sizeof(reservoir[0]));
	rs = 0; /* empty reservoir */

	/* initial items */
	if (nt) {
		fm = 0;
		items[0] = ic[fm++].value;

		for (j = 1, i = 0; j < mis;) {
			ch = fp_grph_children(fp, items[i++], &csz);

			tmp2 = 0;
			for (rsi = 0; rsi < fp->n; rsi++)
				for (ct = tmp2; ct < csz; ct++)
					if (ic[rsi].value == ch[ct]) {
						tmp = ch[tmp2];
						ch[tmp2++] = ch[ct];
						ch[ct] = tmp;
					}

			for (ct = 0; ct < csz && j < mis; ct++) {
				for (rsi = 0; rsi < j; rsi++)
					if (items[rsi] == ch[ct])
						break;
				if (rsi < j)
					continue;
				items[j++] = ch[ct];
			}

			free(ch);
		}
	} else {
		for (fm = 0, j = 0; j < mis && fm < fp->n; fm++)
			items[j++] = ic[fm].value;

		if (j < mis)
			mis = j;
	}

	while (1) {
#if PRINT_RULE_DOMAIN || PRINT_RS_TRACE
		printf("Domain: %lu: ", fm);
		for (i = 0; i < mis; i++)
			printf("%lu ", items[i]);
		printf("\n");
#endif

		generate_and_add_all_rules(fp, items, mis, eps/k,
				&rs, reservoir, k, &randbuffer, minalpha);

		/* move to next set of items */
		if (nt) {
			if (fm == nt)
				break;

			items[0] = ic[fm++].value;

			for (j = 1, i = 0; j < mis;) {
				ch = fp_grph_children(fp, items[i++], &csz);

				tmp2 = 0;
				for (rsi = 0; rsi < fp->n; rsi++)
					for (ct = tmp2; ct < csz; ct++)
						if (ic[rsi].value == ch[ct]) {
							tmp = ch[tmp2];
							ch[tmp2++] = ch[ct];
							ch[ct] = tmp;
						}

				for (ct = 0; ct < csz && j < mis; ct++) {
					for (rsi = 0; rsi < j; rsi++)
						if (items[rsi] == ch[ct])
							break;
					if (rsi < j)
						continue;
					items[j++] = ch[ct];
				}

				free(ch);
			}
		} else {
			if (fm == fp->n)
				break;

			for (i = 0; i < mis - 1; i++)
				items[i] = items[i+1];
			items[mis - 1] = ic[fm++].value;
		}
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

	struct timeval endtime;
	gettimeofday(&endtime, NULL);
	double t1 = starttime.tv_sec + (0.0 + starttime.tv_usec) / MICROSECONDS;
	double t2 = endtime.tv_sec + (0.0 + endtime.tv_usec) / MICROSECONDS;

	printf("Total time: %5.2lf\n", t2 - t1);
	printf("%ld %ld %ld %ld\n", starttime.tv_sec, starttime.tv_usec, endtime.tv_sec, endtime.tv_usec);

	printf("Final histogram:\n");
	histogram_dump(stdout, h, 1, "\t");
#endif

	free_histogram(h);
#if 0 /* moving to graphs */
	free(items);
#endif
	free(ls);
	free(ks);
	free(es);
	free(ic);
}
