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
/* do rule expansion */
#ifndef RULE_EXPAND
#define RULE_EXPAND 0
#endif

static double quality(int x, int y, double m)
{
	double c;

	if (x < m) x = m;
	c = (y + 0.0) / (x + 0.0);

	return m * pow(c, 1);
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
		printf("\t%lu\t%3.2lf\t% 11.5lf\t%7d\t%7d | ", i,
				reservoir[i].c, reservoir[i].v,
				reservoir[i].sxy, reservoir[i].sx);
		print_rule(reservoir[i].r);
		printf("\n");
	}
}
#endif

static void process_rule(const struct fptree *fp,
		const int *AB, int ab_length, const int *A, int a_length,
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
		const int *items, size_t num_items, size_t st, double eps,
		size_t *rs, struct reservoir *reservoir,
		size_t k, struct drand48_data *randbuffer, double m)
{
	size_t i, j, l, max=1<<num_items, a_length, ab_length, max2;
	int *A, *AB;
	double u;

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
					eps, rs, reservoir, k, u, m);
		}
	}

	free(A);
	free(AB);
}

static void split_in_partitions(const struct fptree *fp,
		const struct item_count *ic, int minth,
		size_t ***parts, size_t **parlens, size_t *parcts,
		struct drand48_data *randbuffer)
{
	size_t i, j, k, c;

	/* compute how many items there are */
	for (c = 0; c < fp->n; c++)
		if (ic[c].noisy_count < minth)
			break;
	printf("Items above threshold: %lu\n", c);

	/* TODO: for now is a single partition */
	*parcts = 1;
	*parlens = calloc(*parcts, sizeof(**parlens));
	*parts = calloc(*parcts, sizeof(**parts));

	j = 0;
	for (i = 0; i < *parcts; i++) {
		(*parlens)[i] = c; // TODO: different updates
		(*parts)[i] = calloc((*parlens)[i], sizeof(***parts));
		/* TODO: for now we assume all items go to the single bin */
		for (k = 0; k < *parlens[i]; k++)
			(*parts)[i][k] = j++;
	}
	/* TODO: random ordering of items */
}

void dp2d(const struct fptree *fp, const char *npfile,
		double eps, double eps_share, int minth, size_t mis, size_t k,
		double minalpha, long int seed)
{
	struct item_count *ic = calloc(fp->n, sizeof(ic[0]));
	size_t **partitions = NULL, *parlens = NULL, parcts;
#if 0
	struct histogram *nph = init_histogram();
#endif
	struct histogram *h = init_histogram();
	double epsilon_step1 = eps * eps_share;
	size_t i, j, ip, fm, rs, st, nci, npr;
	int *control_items = NULL, *items;
	struct drand48_data randbuffer;
	double maxc, minc;

	if (strncmp(npfile, "-", strlen("-"))) {
		die("Don't use control items!!!. Pass - for the file!");
#if 0
		FILE *f = fopen(npfile, "r");
		if (!f)
			die("Invalid/not found npfile!");
		if (fscanf(f, "%lu", &nci) != 1)
			die("Invalid npfile: cannot read number of items!");
		control_items = calloc(nci, sizeof(control_items[0]));
		for (i = 0; i < nci; i++)
			if (fscanf(f, "%d", &control_items[i]) != 1)
				die("Invalid npfile: cannot read control items!");
		if (fscanf(f, "%lu", &npr) != 1)
			die("Invalid npfile: cannot read number of non private rules!");
		histogram_load(f, nph, 1, "\t");
		fclose(f);
#endif
	} else {
		nci = 0;
		control_items = calloc(1, sizeof(control_items[0]));
		npr = 0;
	}

	init_rng(seed, &randbuffer);
	items = calloc(mis + nci + 1, sizeof(items[0]));

	printf("Running dp2D with minth=%d, eps=%lf, eps_share=%lf, "
			"mis=%lu, k=%lu\n", minth, eps, eps_share, mis, k);
	printf("Using %lu control items:", nci);
	for (i = 0; i < nci; i++)
		printf(" %d", control_items[i]);
	printf("\nMatching %lu non-private rules\n", npr);

	printf("Step 1: compute noisy counts for items with eps_1 = %lf\n",
			epsilon_step1);
	build_items_table(fp, ic, epsilon_step1, &randbuffer);

#if PRINT_ITEM_TABLE
	printf("\n");
	for (i = 0; i < fp->n; i++)
		printf("%d %d %lf\n", ic[i].value, ic[i].real_count, ic[i].noisy_count);
#endif

	split_in_partitions(fp, ic, minth,
			&partitions, &parlens, &parcts,
			&randbuffer);
	printf("Partitions: %lu |", parcts);
	for (i = 0; i < parcts; i++)
		printf(" %lu", parlens[i]);
	printf("\n");

	eps = (eps - epsilon_step1) / k;
	printf("Step 2: mining with remaining eps (per rule): %lf\n", eps);

	struct timeval starttime;
	gettimeofday(&starttime, NULL);

	for (ip = 0; ip < parcts; ip++) {
		/* select mining domains */
		struct reservoir *reservoir = calloc(k, sizeof(reservoir[0]));
		rs = 0; /* empty reservoir */
		st = 0;

		/* initial items */
		for (i = 0; i  < nci; i++) {
			items[mis + i] = control_items[i];
			st |= 1 << (mis + i);
		}
		for (fm = 0, j = 0; j < mis; fm++) {
#if 0
			for (i = 0; i < nci; i++) {
				if (ic[fm].value == control_items[i])
					break;
			}

			if (i != nci)
				continue;
#endif
			//TODO: remove MM: if (ic[fm].noisy_count > 810) continue;
			if (ic[fm].noisy_count < minth)
				break;
			items[j++] = ic[fm].value;
		}

		if (j < mis)
			mis = j;

		if (mis == 0)
			goto end;

		while (fm < fp->n) {
#if PRINT_RULE_DOMAIN || PRINT_RS_TRACE
			printf("Domain: %lu: ", fm);
			for (i = 0; i < mis + nci; i++)
				printf("%d ", items[i]);
			printf("\n");
#endif

			generate_and_add_all_rules(fp, items, mis + nci, st, eps,
					&rs, reservoir, k, &randbuffer, minalpha);
			st |= (1 << (mis - 1));

			for (i = 0; i < mis - 1; i++)
				items[i] = items[i+1];
			for (; ic[fm].noisy_count >= minth; fm++) {
				for (i = 0; i < nci; i++)
					if (ic[fm].value == control_items[i])
						break;
				if (i != nci)
					continue;
				items[mis - 1] = ic[fm].value;
				break;
			}
			if (ic[fm++].noisy_count < minth)
				break;
		}
end:
		printf("Stopped at fm=%ld item=%d nc=%lf minth=%d\n", fm, ic[fm].value, ic[fm].noisy_count, minth);

#if PRINT_FINAL_RULES
		print_reservoir(reservoir, rs);
#endif

		minc = 1; maxc = 0;
#if RULE_EXPAND
		struct rule_table *rt = init_rule_table();
		for (i = 0; i < rs; i++) {
			struct rule *r = reservoir[i].r, *nr1, *nr2;
			int *items = calloc(r->B->length, sizeof(items[0]));
			struct itemset *A, *AB, *ABprime;
			int supA, supAB;
			size_t k, nis;

			save_rule2(rt, r, reservoir[i].c);

			st = 1 << r->B->length; st--;
			AB = build_itemset_add_items(r->A, r->B->items, r->B->length);
			for (j = 1; j < st; j++) {
				nis = 0;
				for (k = 0; k < r->B->length; k++)
					if ((1 << k) & j)
						items[nis++] = r->B->items[k];

				ABprime = build_itemset_del_items(AB, items, nis);
				nr1 = build_rule_A_AB(r->A, ABprime);
				supAB = fpt_itemset_count(fp, ABprime->items, ABprime->length);
				supA = fpt_itemset_count(fp, r->A->items, r->A->length);
				save_rule2(rt, nr1, supAB / (supA + 0.0));

				A = build_itemset_add_items(r->A, items, nis);
				nr2 = build_rule_A_AB(A, AB);
				supAB = fpt_itemset_count(fp, AB->items, AB->length);
				supA = fpt_itemset_count(fp, A->items, A->length);
				save_rule2(rt, nr2, supAB / (supA + 0.0));
			}
		}

		for (i = 0; i < rt->sz; i++) {
			if (rt->c[i] < minc)
				minc = rt->c[i];
			if (rt->c[i] > maxc)
				maxc = rt->c[i];
			histogram_register(h, rt->c[i]);
		}

		rs = rt->sz;
#else
		for (i = 0; i < rs; i++) {
			if (reservoir[i].c < minc)
				minc = reservoir[i].c;
			if (reservoir[i].c > maxc)
				maxc = reservoir[i].c;
			histogram_register(h, reservoir[i].c);
		}
#endif
		printf("Rules saved: %lu, minconf: %3.2lf, maxconf: %3.2lf\n",
				rs, minc, maxc);

		free_reservoir_array(reservoir, rs);
	}

	struct timeval endtime;
	gettimeofday(&endtime, NULL);
	double t1 = starttime.tv_sec + (0.0 + starttime.tv_usec) / MICROSECONDS;
	double t2 = endtime.tv_sec + (0.0 + endtime.tv_usec) / MICROSECONDS;

	printf("Total time: %5.2lf\n", t2 - t1);
	printf("%ld %ld %ld %ld\n", starttime.tv_sec, starttime.tv_usec, endtime.tv_sec, endtime.tv_usec);

	printf("Final histogram:\n");
	histogram_dump(stdout, h, 1, "\t");

#if 0
	printf("Non-private histogram:\n");
	histogram_dump(stdout, nph, 1, "\t");
#endif

	for (i = 0; i < parcts; i++)
		free(partitions[i]);
	free(parlens);
	free(partitions);
	free_histogram(h);
	free(control_items);
	free(items);
	free(ic);
}
