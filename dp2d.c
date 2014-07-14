#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>
#include <mpfr.h>

#include "dp2d.h"
#include "fp.h"
#include "globals.h"

/* debug defines */
#ifndef PRINT_TRIANGLES_ITEMS
#define PRINT_TRIANGLES_ITEMS 0
#endif
#ifndef PRINT_TRIANGLE_CONTENT
#define PRINT_TRIANGLE_CONTENT 1
#endif
#ifndef PRINT_ITEM_TABLE
#define PRINT_ITEM_TABLE 0
#endif
#ifndef PRINT_RULE_TABLE
#define PRINT_RULE_TABLE 1
#endif

#define PRECISION 2048
#define ROUND_MODE MPFR_RNDN

struct item_count {
	int value;
	int real_count;
	double noisy_count;
};

static double quality(int x, int y, double X, double Y);

static struct item_count *alloc_items(int sz)
{
	return calloc(sz, sizeof(struct item_count));
}

static void free_items(const struct item_count *ic)
{
	free((void*)ic);
}

static int ic_noisy_cmp(const void *a, const void *b)
{
	const struct item_count *ia = a, *ib = b;
	return double_cmp_r(&ia->noisy_count, &ib->noisy_count);
}

static void build_items_table(const struct fptree *fp, struct item_count *ic,
		double eps, struct drand48_data *buffer)
{
	int i;

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
	struct item_count x;
	int nm;

	/* get one more item */
	x.noisy_count = *m;
	nm = bsearch_i(&x, ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);
	*m = floor(ic[nm].noisy_count - 1);
	*M = *m / c;
}

static int get_triangles(const struct fptree *fp, const struct item_count *ic,
		double c, int m, int M, int minth)
{
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
}

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

static void all_allowed_items(const struct fptree *fp,
		const struct item_count *ic, int m, int M, int *items)
{
	int i = 0;
	struct item_count x;
	int nm, nM;

	x.noisy_count = m;
	nm = bsearch_i(&x, ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);
	x.noisy_count = M;
	nM = bsearch_i(&x, ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);

	for (i = nM;  i < nm; i++)
		items[i - nM] = ic[i].value;
}

static void compute_pdf(const struct fptree *fp, int m, int M, double epsilon,
		const int *A, const int *AB,
		int a_length, int ab_length,
		int *sup_a, int *sup_ab, double *q, mpfr_t pdf)
{
	*sup_a = fpt_itemset_count(fp, A, a_length);
	*sup_ab = fpt_itemset_count(fp, AB, ab_length);
	*q = quality(*sup_a, *sup_ab, M, m);
	mpfr_set_d(pdf, *q, ROUND_MODE);
	mpfr_mul_d(pdf, pdf, epsilon / 2, ROUND_MODE);
	mpfr_exp(pdf, pdf, ROUND_MODE);
}

static void print_rule(const int *A, const int *B, int a_length, int b_length)
{
	int i;

	for (i = 0; i < a_length; i++)
		printf("%d ", A[i]);
	printf("-> ");
	for (i = 0; i < b_length; i++)
		if (B[i])
			printf("%d ", B[i]);
}

static void do_print_rule(const struct fptree *fp,
		const int *A, const int *B, const int *AB,
		int a_length, int b_length, int ab_length)
{
	int sup_a, sup_ab;

	sup_a = fpt_itemset_count(fp, A, a_length);
	sup_ab = fpt_itemset_count(fp, AB, ab_length);

	printf("\n");
	print_rule(A, B, a_length, b_length);
	printf(" (%d, %d) %lf", sup_a, sup_ab,
			(sup_ab + 0.0) / (0.0 + sup_a));

}

static void print_rule_and_expansions(const struct fptree *fp,
		const int *A, const int *B,
		int a_length, int b_length)
{
	int *ab = calloc(a_length + b_length, sizeof(*ab));
	int *a = calloc(a_length + b_length, sizeof(*a));
	int *b = calloc(a_length + b_length, sizeof(*b));
	char *x = calloc(b_length, sizeof(*b));
	int i, j = 0;

	for (i = 0; i < a_length; i++) {
		a[i] = A[i];
		ab[i] = A[i];
	}

	for (i = 0; i < b_length; i++) {
		b[i] = B[i];
		ab[i + a_length] = B[i];
	}

	do_print_rule(fp, a, b, ab, a_length, b_length, a_length + b_length);
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

		do_print_rule(fp, a, b, ab, a_length, b_length, a_length + b_length);

		j = a_length;

		/* add item to a */
		for (i = 0; i < b_length; i++)
			if (x[i]) {
				a[a_length++] = B[i];
				ab[i + j] = B[i];
			}

		do_print_rule(fp, a, b, ab, a_length, b_length, j + b_length);

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

static void compute_cdf(const struct fptree *fp, int m, int M, double epsilon,
		const int *A, const int *B, const int *AB,
		int a_length, int b_length, int ab_length,
		mpfr_t pdf, mpfr_t cdf)
{
	int sup_a, sup_ab;
	double q;

	compute_pdf(fp, m, M, epsilon, A, AB, a_length, ab_length,
			&sup_a, &sup_ab, &q, pdf);
	mpfr_add(cdf, cdf, pdf, ROUND_MODE);

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
}

static void sample_rule(const struct fptree *fp, int m, int M, double epsilon,
		const int *A, const int *B, const int *AB,
		int a_length, int b_length, int ab_length,
		mpfr_t pdf, mpfr_t rnd)
{
	int sup_a, sup_ab;
	double q;

	if (mpfr_sgn(rnd) < 0)
		return; /* one rule was already sampled */

	compute_pdf(fp, m, M, epsilon, A, AB, a_length, ab_length,
			&sup_a, &sup_ab, &q, pdf);
	mpfr_sub(rnd, rnd, pdf, ROUND_MODE);
	if (mpfr_sgn(rnd) < 0)
		print_rule_and_expansions(fp, A, B, a_length, b_length);

	B = B;
	b_length = b_length;
}

static void for_all_rules(const struct fptree *fp, const int *items,
		int num_items, int m, int M, double epsilon,
		mpfr_t arg1, mpfr_t arg2,
		void (*f)(const struct fptree*, int, int, double,
			const int*, const int*, const int*, int, int, int,
			mpfr_t, mpfr_t))
{
#define SHORT 1 /* item -. anything, O(n^2) */
#define LONG 2 /* anything -> anything, O(n^3) */
#define METHOD SHORT

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

			f(fp, m, M, epsilon, A, B, AB, 1, b_length, ab_length,
					arg1, arg2);
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

		f(fp, m, M, epsilon, A, B, AB, a_length, b_length, ab_length,
				arg1, arg2);
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

static void mine(const struct fptree *fp, const struct item_count *ic,
		double c, double epsilon, int num_triangles,
		int m, int M)
{
	int c_triangle, num_items, *items;
	gmp_randstate_t rnd_state;
	mpfr_t cdf, pdf, rnd;
	double c_epsilon;

	mpfr_init2(pdf, PRECISION);
	mpfr_init2(cdf, PRECISION);
	mpfr_init2(rnd, PRECISION);
	initialize_random(rnd_state, PRECISION * 200 * num_triangles);

	for (c_triangle = 0; c_triangle < num_triangles; c_triangle++) {
		/* TODO: split based on triangles */
		c_epsilon = epsilon / num_triangles;

		num_items = num_allowed_items(fp, ic, m, M);
#if PRINT_TRIANGLE_CONTENT == 1
		printf("\nTriangle %d %d %d: %d", c_triangle, m, M, num_items);
#endif
		if (num_items < 2)
			goto end_loop;

		items = calloc(num_items, sizeof(items[0]));
		all_allowed_items(fp, ic, m, M, items);
#if PRINT_TRIANGLE_CONTENT == 1
		int i;

		printf("| ");
		for (i = 0; i < num_items; i++)
			printf("%d ", items[i]);
#endif

		mpfr_set_zero(cdf, 1);
		for_all_rules(fp, items, num_items, m, M, c_epsilon,
				pdf, cdf, compute_cdf);

		mpfr_urandomb(rnd, rnd_state);
		mpfr_mul(rnd, rnd, cdf, ROUND_MODE); /* rnd in [0, cdf) */
		for_all_rules(fp, items, num_items, m, M, c_epsilon,
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
}

void dp2d(const struct fptree *fp, double c, double eps, double eps_share,
		int ni, int minth)
{
	struct item_count *ic = alloc_items(fp->n);
	double epsilon_step1 = eps * eps_share;
	struct drand48_data randbuffer;
	int theta, nt;

	init_rng(&randbuffer);

	printf("Running dp2D with ni=%d, minth=%d, c=%lf, eps=%lf, "
			"eps_share=%lf\n", ni, minth, c, eps, eps_share);

	printf("Step 1: compute noisy counts for items: eps_1 = %lf\n",
			epsilon_step1);
	build_items_table(fp, ic, epsilon_step1, &randbuffer);

#if PRINT_ITEM_TABLE == 1
	printf("\n");
	int i;
	for (i = 0; i < fp->n; i++)
		printf("%d %d %lf\n", ic[i].value, ic[i].real_count, ic[i].noisy_count);
#endif

	printf("Step 2: get min support for %d items: ", ni);
	theta = get_first_theta_value(fp, ic, ni, c, c * fp->t);
	printf("%d\n", theta);

	printf("Step 3: get triangles: ");
	nt = get_triangles(fp, ic, c, theta, fp->t, minth);
	printf("%d triangles needed\n", nt);

	printf("Step 4: mining rules:\n");
	mine(fp, ic, c, eps - epsilon_step1, nt, theta, fp->t);

	free_items(ic);
}

#define DISPLACEMENT 1
#define DISTANCE 2
#define METHOD DISPLACEMENT
#if METHOD == DISPLACEMENT
static double displacement(int x, double m, double M, double v)
{
	double z = 1 - fabs(x - v) / (M - m);

	if (signbit(z))
		return 0;

	return v * z;
}

static double quality(int x, int y, double X, double Y)
{
	return 0.5 * (1/Y - 1/X) *
		displacement(x, Y, X, X) * displacement(y, Y, X, Y);
}
#elif METHOD == DISTANCE
static double quality(int x, int y, double X, double Y)
{
	double dx = X - x, dy = Y - y;
	return sqrt(dx * dx + dy * dy);
}
#else
#error ("Undefined quality method")
#endif
#undef DISPLACEMENT
#undef DISTANCE
#undef METHOD
