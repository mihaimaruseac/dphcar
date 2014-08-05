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
#define PRINT_ITEM_DOMAIN 1
#endif
/* print the rules generated at each step and their quality */
#ifndef PRINT_RULE_DOMAIN
#define PRINT_RULE_DOMAIN 1
#endif
#if 0
/* debug defines */
#ifndef PRINT_TRIANGLES_ITEMS
#define PRINT_TRIANGLES_ITEMS 0
#endif
#ifndef PRINT_TRIANGLE_CONTENT
#define PRINT_TRIANGLE_CONTENT 0
#endif
#ifndef PRINT_RULE_TABLE
#define PRINT_RULE_TABLE 0
#endif
#ifndef PRINT_RULE_RESULTS
#define PRINT_RULE_RESULTS 0
#endif
#ifndef RULE_LEN_LIMIT
#define RULE_LEN_LIMIT 1
#endif
#ifndef RULE_LEN_LIMITED_STEP
#define RULE_LEN_LIMITED_STEP 100
#endif
#ifndef RULE_ITEMSET
#define RULE_ITEMSET 0
#endif
#ifndef UNIFORM_SAMPLE
#define UNIFORM_SAMPLE 1
#endif

#define PRECISION 2048
#define ROUND_MODE MPFR_RNDN

static double quality(int x, int y, double X, double Y, double c);

/**
 * Reads `bytes` bytes from /dev/urandom to initialize random number generator
 * state.
 */
static char *read_from_devurandom(int bytes)
{
	FILE *f = fopen("/dev/urandom", "r");
	char *buff = malloc(bytes);
#if 0
	// TODO: enable seeding
	char *p = buff;
	size_t read;

	while (bytes) {
		read = fread(p, 1, bytes, f);
		p += read;
		bytes -= read;
	}
#else
	memset(buff, 'a', bytes);
#endif

	fclose(f);
	return buff;
}

static void initialize_random(gmp_randstate_t state, int bytes)
{
	void *buff = read_from_devurandom(bytes);
	mpz_t s;

	gmp_randinit_default(state);
	mpz_init(s);
	mpz_import(s, bytes, 1, 1, 0, 0, buff);
	gmp_randseed(state, s);
	mpz_clear(s);
	free(buff);
}
#endif

struct item_count {
	int value;
	int real_count;
	double noisy_count;
};

struct reservoir {
	struct rule *r;
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

#if 0
static int get_first_theta_value(const struct fptree *fp,
		const struct item_count *ic, int ni, double c, double try)
{
	struct item_count x;
	int ix;

	x.noisy_count = try;
	ix = bsearch_i(&x, ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);

	while (ix < ni) {
		x.noisy_count *= c;
		ix = bsearch_i(&x, ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);
	}

	return floor(x.noisy_count);
}

static void get_next_triangle(const struct fptree *fp,
		const struct item_count *ic, double c, int *m, int *M)
{
#if RULE_LEN_LIMIT
	*m -= RULE_LEN_LIMITED_STEP;
	fp = fp;
	ic = ic;
#else
	struct item_count x;
	int nm;

	/* get one more item */
	x.noisy_count = *m;
	nm = bsearch_i(&x, ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);
	*m = floor(ic[nm].noisy_count - 1);
#endif
	*M = *m / c;
}

static int get_triangles(const struct fptree *fp, const struct item_count *ic,
		double c, int m, int M, int minth)
{
#if RULE_LEN_LIMIT
	fp = fp; ic = ic; c = c; M = M;
	return (m - minth) / RULE_LEN_LIMITED_STEP;
#else
	int num_triangles = 0;

	if (m <= minth)
		die("Minimum threshold set too high");

#if PRINT_TRIANGLES_ITEMS
	struct item_count x;
	int nm, nM;

	printf("\n");

	x.noisy_count = M;
	nM = bsearch_i(&x, ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);
	x.noisy_count = m;
	nm = bsearch_i(&x, ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);
#endif

	while (m > minth) {
#if PRINT_TRIANGLES_ITEMS
		printf("triangle %d (%d %d) (items %d %d) (noise %lf %lf) "
				"inside: %d\n",
				num_triangles, M, m, nM, nm,
				ic[nM].noisy_count, ic[nm].noisy_count,
				nm - nM);

		x.noisy_count = M;
		nM = bsearch_i(&x, ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);
		x.noisy_count = m;
		nm = bsearch_i(&x, ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);
#endif

		get_next_triangle(fp, ic, c, &m, &M);
		num_triangles++;
	}

	return num_triangles;
#endif
}

#if RULE_LEN_LIMIT
#else
static int num_allowed_items(const struct fptree *fp,
		const struct item_count *ic, int m, int M)
{
	struct item_count x;
	int nm;

	x.noisy_count = m;
	nm = bsearch_i(&x, ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);
	x.noisy_count = M;
	return nm - bsearch_i(&x, ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);
}
#endif

static void all_allowed_items(const struct fptree *fp,
		const struct item_count *ic, int m, int M, int *items)
{
	int i = 0;
	struct item_count x;
	int nm, nM;

	x.noisy_count = m;
	nm = bsearch_i(&x, ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);
#if RULE_LEN_LIMIT
	nM = nm - M;
	if (nM < 0) nM = 0;
#else
	x.noisy_count = M;
	nM = bsearch_i(&x, ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);
#endif

	for (i = nM;  i < nm; i++)
		items[i - nM] = ic[i].value;
}

static void do_output_rule(const struct fptree *fp,
		struct rule_table *rt, struct histogram *h,
		const int *A, const int *AB,
		int a_length, int ab_length,
		int m, int M)
{
	struct itemset *a, *ab;
	int sup_a, sup_ab;
	struct rule *r;

	sup_a = fpt_itemset_count(fp, A, a_length);
	sup_ab = fpt_itemset_count(fp, AB, ab_length);

	a = build_itemset(A, a_length);
	ab = build_itemset(AB, ab_length);
	r = build_rule_A_AB(a, ab);

#if PRINT_RULE_RESULTS
	printf("\n");
	if (sup_ab > M)
		printf("a ");
	else if (sup_ab >= m)
		if (sup_a > M)
			printf("r ");
		else
			printf("  ");
	else if (sup_a > M)
		printf("+ ");
	else if (sup_a >= m)
		printf("b ");
	else
		printf("* ");

	print_rule(r);

	printf(" (%d, %d) %lf", sup_a, sup_ab,
			(sup_ab + 0.0) / (0.0 + sup_a));
#else
	m = M = m;
#endif

	if (rt)
		save_rule(rt, r, sup_a, sup_ab);
	else  {/* if (h) */
		histogram_register(h, (0.0 + sup_ab) / (0.0 + sup_a));
		free_rule(r);
	}

	free_itemset(a);
	free_itemset(ab);
}

static void output_rule_and_expansions(const struct fptree *fp,
		struct rule_table *rt, struct histogram *h,
		const int *A, const int *B,
		int a_length, int b_length, int m, int M)
{
	int *ab = calloc(a_length + b_length, sizeof(*ab));
	int *a = calloc(a_length + b_length, sizeof(*a));
	int *b = calloc(a_length + b_length, sizeof(*b));
	char *x = calloc(b_length, sizeof(*b));
	int i, j = 0;

	if (rt == NULL && h == NULL)
		die("Either rule table or histogram must be non NULL");
	if (rt != NULL && h != NULL)
		die("Either rule table or histogram must be NULL");

	for (i = 0; i < a_length; i++) {
		a[i] = A[i];
		ab[i] = A[i];
	}

	for (i = 0; i < b_length; i++) {
		b[i] = B[i];
		ab[i + a_length] = B[i];
	}

	do_output_rule(fp, rt, h, a, ab, a_length, a_length + b_length, m, M);
	x[b_length - 1]++;

	while (x[0] < 2) {
		/* remove from b */
		for (i = 0; i < b_length; i++)
			if (x[i]) {
				b[i] = 0;
				ab[i + a_length] = 0;
			}

		j = 0;
		for (i = 0; i < b_length; i++)
			if (b[i])
				j++;

		if (j == 0)
			goto next;

		do_output_rule(fp, rt, h, a, ab, a_length, a_length + b_length, m, M);

		j = a_length;

		/* add item to a */
		for (i = 0; i < b_length; i++)
			if (x[i]) {
				a[a_length++] = B[i];
				ab[i + j] = B[i];
			}

		do_output_rule(fp, rt, h, a, ab, a_length, j + b_length, m, M);

		/* remove item from a */
		for (i = j; i < a_length; i++)
			a[i] = 0;
		a_length = j;

next:
		/* add item back to b */
		for (i = 0; i < b_length; i++)
			if (x[i]) {
				b[i] = B[i];
				ab[i + a_length] = B[i];
			}

		x[b_length - 1]++;
		if (x[0] == 2) break;
		for (i = b_length - 1; i > 0; i--)
			if (x[i] == 2) {
				x[i] = 0;
				x[i-1]++;
			}
	}

	free(a);
	free(ab);
	free(b);
	free(x);
}

static void compute_pdf(const struct fptree *fp, int m, int M,
		double c, double epsilon,
		const int *A, const int *AB,
		int a_length, int ab_length,
		int *sup_a, int *sup_ab, double *q, mpfr_t pdf)
{
	*sup_a = fpt_itemset_count(fp, A, a_length);
	*sup_ab = fpt_itemset_count(fp, AB, ab_length);
	*q = quality(*sup_a, *sup_ab, M, m, c);
#if RULE_ITEMSET
	int i;
	for (i = 0; i < ab_length; i++)
		*q *= fpt_item_score(fp, AB[i]);
#endif
	mpfr_set_d(pdf, *q, ROUND_MODE);
	mpfr_mul_d(pdf, pdf, epsilon / 2, ROUND_MODE);
	mpfr_exp(pdf, pdf, ROUND_MODE);
}

static void compute_cdf(const struct fptree *fp, struct rule_table *rt,
		int m, int M,
		double c, double epsilon,
		const int *A, const int *B, const int *AB,
		int a_length, int b_length, int ab_length,
		mpfr_t pdf, mpfr_t cdf)
{
	int sup_a, sup_ab;
	double q;

	compute_pdf(fp, m, M, c, epsilon, A, AB, a_length, ab_length,
			&sup_a, &sup_ab, &q, pdf);
#if UNIFORM_SAMPLE
	mpfr_add_ui(cdf, cdf, 1, ROUND_MODE);
#else
	mpfr_add(cdf, cdf, pdf, ROUND_MODE);
#endif

#if PRINT_RULE_TABLE == 1
	printf("\n");
	print_rule(A, B, a_length, ab_length);
	printf(" (%d, %d) %lf ", sup_a, sup_ab, q);
	mpfr_out_str(stdout, 10, 5, pdf, ROUND_MODE);
	printf(" ");
	mpfr_out_str(stdout, 10, 5, cdf, ROUND_MODE);
#else
	B = B;
#endif
	b_length = b_length;
	rt = rt; /* unused */
}

static void sample_rule(const struct fptree *fp, struct rule_table *rt,
		int m, int M,
		double c, double epsilon,
		const int *A, const int *B, const int *AB,
		int a_length, int b_length, int ab_length,
		mpfr_t pdf, mpfr_t rnd)
{
	int sup_a, sup_ab;
	double q;

	if (mpfr_sgn(rnd) < 0)
		return; /* one rule was already sampled */

	compute_pdf(fp, m, M, c, epsilon, A, AB, a_length, ab_length,
			&sup_a, &sup_ab, &q, pdf);
#if UNIFORM_SAMPLE
	mpfr_sub_ui(rnd, rnd, 1, ROUND_MODE);
#else
	mpfr_sub(rnd, rnd, pdf, ROUND_MODE);
#endif
	if (mpfr_sgn(rnd) < 0)
		output_rule_and_expansions(fp, rt, NULL, A, B, a_length, b_length, m, M);
}

#define SHORT 1 /* item -> anything, O(n^2) */
#define LONG 2 /* anything -> anything, O(n^3) */
#define METHOD LONG
static void for_all_rules(const struct fptree *fp, struct rule_table *rt,
		const int *items, int num_items, int m, int M,
		double c, double epsilon, mpfr_t arg1, mpfr_t arg2,
		void (*f)(const struct fptree*, struct rule_table*, int, int,
			double, double, const int*, const int*, const int*,
			int, int, int, mpfr_t, mpfr_t))
{
	int i, a_length, b_length, ab_length;
	unsigned char *ABi;
	int *A, *B, *AB;

	ABi = calloc(num_items, sizeof(ABi[0]));
	AB = calloc(num_items, sizeof(AB[0]));
	A = calloc(num_items, sizeof(A[0]));
	B = calloc(num_items, sizeof(B[0]));

#if METHOD == SHORT
#define RULE_END 2
#elif METHOD == LONG
#define RULE_A 1
#define RULE_B 2
#define RULE_END 3
	ABi[num_items - 1] = 1;
	ABi[num_items - 2] = 1;
#endif

	while (1) {
		ABi[num_items - 1]++;
		a_length = b_length = ab_length = 0;

		for (i = num_items - 1; i > 0; i--)
			if (ABi[i] == RULE_END) {
				ABi[i] = 0;
				ABi[i-1]++;
		}

		if (ABi[0] == RULE_END) break;

#if METHOD == SHORT
		a_length = a_length;
		for (i = 0; i < num_items; i++)
			if (ABi[i])
				AB[ab_length++] = items[i];

		if (ab_length < 2) continue;

		for (i = 0; i < ab_length; i++) {
			int j;

			A[0] = AB[i];
			b_length = 0;
			for (j = 0; j < ab_length; j++)
				if (i != j)
					B[b_length++] = AB[j];

			f(fp, rt, m, M, c, epsilon, A, B, AB,
				1, b_length, ab_length, arg1, arg2);
		}

#elif METHOD == LONG
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

		f(fp, rt, m, M, c, epsilon, A, B, AB,
			a_length, b_length, ab_length, arg1, arg2);
#endif
	}

	free(ABi);
	free(AB);
	free(A);
	free(B);
#undef RULE_A
#undef RULE_B
#undef RULE_END
#undef METHOD
#undef SHORT
#undef LONG
}

static struct rule_table *mine(const struct fptree *fp,
		const struct item_count *ic,
		double c, double epsilon, int num_triangles,
		int m, int M, int ni)
{
	struct rule_table *rt = init_rule_table();
	int c_triangle, num_items, *items;
	gmp_randstate_t rnd_state;
	mpfr_t cdf, pdf, rnd;
	double c_epsilon;

	ni = ni;
	mpfr_init2(pdf, PRECISION);
	mpfr_init2(cdf, PRECISION);
	mpfr_init2(rnd, PRECISION);
	initialize_random(rnd_state, PRECISION * num_triangles);

	for (c_triangle = 0; c_triangle < num_triangles; c_triangle++) {
		/* TODO: split based on triangles */
		c_epsilon = epsilon / num_triangles;

#if RULE_LEN_LIMIT
		num_items = ni;
#else
		num_items = num_allowed_items(fp, ic, m, M);
#endif
#if PRINT_TRIANGLE_CONTENT == 1
		printf("\nTriangle %d %d %d: %d", c_triangle, m, M, num_items);
#endif
		if (num_items < 2)
			goto end_loop;

		items = calloc(num_items, sizeof(items[0]));
#if RULE_LEN_LIMIT
		all_allowed_items(fp, ic, m, ni, items);
#else
		all_allowed_items(fp, ic, m, M, items);
#endif
#if PRINT_TRIANGLE_CONTENT == 1
		int i;

		printf("| ");
		for (i = 0; i < num_items; i++)
			printf("%d ", items[i]);
#endif

		mpfr_set_zero(cdf, 1);
		for_all_rules(fp, rt, items, num_items, m, M, c, c_epsilon,
				pdf, cdf, compute_cdf);

		mpfr_urandomb(rnd, rnd_state);
		mpfr_mul(rnd, rnd, cdf, ROUND_MODE); /* rnd in [0, cdf) */
		for_all_rules(fp, rt, items, num_items, m, M, c, c_epsilon,
				pdf, rnd, sample_rule);

		free(items);

end_loop:
#if PRINT_TRIANGLE_CONTENT == 1
		printf("\n");
#endif
		get_next_triangle(fp, ic, c, &m, &M);
	}

	mpfr_clear(pdf);
	mpfr_clear(cdf);
	mpfr_clear(rnd);
	mpfr_free_cache();
	gmp_randclear(rnd_state);

	return rt;
}

static void get_rules_np(const struct fptree *fp,
		int *itemset, int nits, struct histogram *hnp)
{
	int *A = calloc(1, sizeof(A[0]));
	int i;

	for (i = 0; i < nits; i++) {
		A[0] = itemset[i];
		itemset[i] = 0;
		output_rule_and_expansions(fp, NULL, hnp, A, itemset,
				1, nits - 1, 0, 0);
		itemset[i] = A[0];
	}

	free(A);
}

static void mine_np(const struct fptree *fp, int ni,
		const char *ifname, int *top_items, int hic,
		struct histogram *hnp)
{
	int *itemset = calloc(ni, sizeof(itemset[0])), nits, tmp;
	FILE *f = fopen(ifname, "r");
	int i, j;

	if (!f || fscanf(f, "(%d)", &tmp) != 1 || tmp != fp->t)
		die("Invalid/non-matching transaction support file!");


	while (1) {
		for (nits = 0; nits < ni; nits++)
			if(fscanf(f, "%d", &itemset[nits]) != 1)
				break;

		if (fscanf(f, "%d", &tmp) == 1) {
			/* extra items, read till end and continue to next
			 * while-loop
			 */
			while(fscanf(f, "%d", &tmp) == 1);
			goto end;
		}

		if (nits == 0) /* eof */
			break;

		tmp = 0;
		for (i = 0; i < hic; i++)
			for (j = 0; j < nits; j++)
				if (top_items[i] == itemset[j]) {
					tmp++;
					break;
				}

		if (tmp != hic)
			goto end;

		if (nits == 1)
			goto end;

		get_rules_np(fp, itemset, nits, hnp);

end:
		if (fscanf(f, "(%d)", &tmp) != 1)
			die("Invalid support file :: Extra line data");
	}

	fclose(f);
	free(itemset);
}
#endif

static void generate_and_add_all_rules(const struct fptree *fp,
		const int *items, size_t num_items,
		size_t *rs, struct reservoir *reservoir,
		struct drand48_data *randbuffer)
{
#define RULE_A 1
#define RULE_B 2
#define RULE_END 3

	size_t i, a_length, b_length, ab_length;
	struct itemset *iA, *iAB;
	struct rule *r = NULL;
	unsigned char *ABi;
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

		// TODO: print rule first
		free_rule(r);
		iA = build_itemset(A, a_length);
		iAB = build_itemset(AB, ab_length);
		r = build_rule_A_AB(iA, iAB);
		free_itemset(iA);
		free_itemset(iAB);

#if PRINT_RULE_DOMAIN
		print_rule(r);
		// TODO: supports and quality
		printf("\n");
#endif

#if 0
		f(fp, rt, m, M, c, epsilon, A, B, AB,
			a_length, b_length, ab_length, arg1, arg2);
#endif
	}

	free_rule(r);
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

	init_rng(&randbuffer);

	printf("Running dp2D with minth=%d, eps=%lf, eps_share=%lf, mu=%lf, "
			"mis=%lu, k=%lu\n", minth, eps, eps_share, mu, mis, k);

	printf("Step 1: compute noisy counts for items with eps_1 = %lf: ",
			epsilon_step1);
	build_items_table(fp, ic, epsilon_step1, &randbuffer);
	eps = (eps - epsilon_step1) / k;
	printf("remaining eps (per rule): %lf\n", eps);

#if PRINT_ITEM_TABLE
	printf("\n");
	for (i = 0; i < fp->n; i++)
		printf("%d %d %lf\n", ic[i].value, ic[i].real_count, ic[i].noisy_count);
#endif

	/* select mining domains */
	rs = 0; /* empty reservoir */
	fm = fM = 0;
	fm = update_fm(fm, fM, mis, mu, fp->n, minth, ic);
	while (fm != fM) {
#if PRINT_ITEM_DOMAIN || PRINT_RULE_DOMAIN
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

#if PRINT_ITEM_DOMAIN || PRINT_RULE_DOMAIN
		printf("\n");
#endif

		/* TODO: generate and keep the rules in this domain */
		generate_and_add_all_rules(fp, items, 1 + fm - fM,
				&rs, reservoir, &randbuffer);

		fM++;
		fm = update_fm(fm, fM, mis, mu, fp->n, minth, ic);
	}

	/* TODO: compute min c inside the array of results (simulate outputing rules) */

	free_reservoir_array(reservoir, k);
	free(items);
	free(ic);
}

#if 0
#define DISPLACEMENT 1
#define DISTANCE 2
#define NEW 3
#define DELTA 4
#define CDIST 5
#define METHOD DISPLACEMENT
#if METHOD == DISPLACEMENT
static double displacement(int x, double m, double M, double v)
{
	double z = 1 - fabs(x - v) / (M - m);

	if (signbit(z))
		return 0;

	return v * z;
}

static double quality(int x, int y, double X, double Y, double c)
{
	c = c;
	return 0.5 * (1/Y - 1/X) *
		displacement(x, Y, X, X) * displacement(y, Y, X, Y);
}
#elif METHOD == DISTANCE
static double quality(int x, int y, double X, double Y, double c)
{
	double dx = X - x, dy = Y - y;
	c = c;
	return sqrt(dx * dx + dy * dy);
}
#elif METHOD == NEW
static double quality(int x, int y, double X, double Y, double c)
{
	c = c;
#define P1 1
#define P2 1
#define P3 1
	double s1, s2, s, x1, y1;

	s  = pow(1 + 1 / X, P3);
	s1 = pow(1 + 1 / Y, P2);
	s2 = 2 - pow(1 - 1 / (X - Y), P1);

	if (s1 > s2)
		s *= s1;
	else
		s *= s2;

	x1 = x / X;
	if (y < Y)
		y1 = pow(y / Y, P2);
	else
		y1 = pow((X - y) / (X - Y), P1);

	return x1 * y1 / s;
#undef P1
#undef P2
#undef P3
}
#elif METHOD == DELTA
static double quality(int x, int y, double X, double Y, double c)
{
	X = X; Y = Y; c = c;
	return y - x;
}
#elif METHOD == CDIST
static double quality(int x, int y, double X, double Y, double c)
{
	double s = c;
	X = X; Y = Y;

	if (1 - c > s)
		s = 1 - c;
	s /= sqrt(1 + c*c);

	return fabs(c * x - y) / s;
}
#else
#error ("Undefined quality method")
#endif
#undef DISPLACEMENT
#undef DISTANCE
#undef METHOD
#endif
