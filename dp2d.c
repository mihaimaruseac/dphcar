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
#define PRINT_TRIANGLE_CONTENT 0
#endif
#ifndef PRINT_ITEM_TABLE
#define PRINT_ITEM_TABLE 0
#endif
#ifndef PRINT_RULE_TABLE
#define PRINT_RULE_TABLE 0
#endif

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

static void free_items(struct item_count *ic)
{
	free(ic);
}

static int ic_noisy_cmp(const void *a, const void *b)
{
	const struct item_count *ia = a, *ib = b;
	return double_cmp_r(&ia->noisy_count, &ib->noisy_count);
}

static void build_items_table(struct fptree *fp, struct item_count *ic,
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

static int get_first_theta_value(struct fptree *fp, struct item_count *ic,
		int ni, double c, double try)
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

static void get_next_triangle(struct fptree *fp, struct item_count *ic,
		double c, int *m, int *M)
{
	struct item_count x;
	int nm;

	/* get one more item */
	x.noisy_count = *m;
	nm = bsearch_i(&x, ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);
	*m = floor(ic[nm].noisy_count - 1);
	*M = *m / c;
}

static int get_triangles(struct fptree *fp, struct item_count *ic, double c,
		int m, int M, int minth)
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

static int num_allowed_items(struct fptree *fp, struct item_count *ic,
		int m, int M)
{
	struct item_count x;
	int nm;

	x.noisy_count = m;
	nm = bsearch_i(&x, ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);
	x.noisy_count = M;
	return nm - bsearch_i(&x, ic, fp->n, sizeof(ic[0]), ic_noisy_cmp);
}

static void all_allowed_items(struct fptree *fp, struct item_count *ic,
		int m, int M, int *items)
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

static void compute_cdf(struct fptree *fp, int m, int M, double epsilon,
		int *A, int *B, int a_length, int b_length)
{
	int i;

#if PRINT_RULE_TABLE == 1
	for (i = 0; i < a_length; i++)
		printf("%d ", A[i]);
	printf("-> ");
	for (i = 0; i < b_length; i++)
		printf("%d ", B[i]);
	printf("\n");
#endif
}

static void for_all_rules(struct fptree *fp, int *items, int num_items,
		int m, int M, double epsilon,
		void (*f)(struct fptree *, int, int, double,
			int*, int *, int, int))
{
#define RULE_A 1
#define RULE_B 2
#define RULE_END 3
	int i, a_length, b_length;
	unsigned char *AB;
	int *A, *B;

	AB = calloc(num_items, sizeof(AB[0]));
	A = calloc(num_items, sizeof(A[0]));
	B = calloc(num_items, sizeof(B[0]));
	AB[num_items - 1] = 1;
	AB[num_items - 2] = 1;

	while (1) {
		AB[num_items - 1]++;
		a_length = b_length = 0;

		for (i = num_items - 1; i > 0; i--)
			if (AB[i] == RULE_END) {
				AB[i] = 0;
				AB[i-1]++;
			}

		for (i = 0; i < num_items; i++)
			switch (AB[i]) {
			case RULE_B: B[b_length++] = items[i]; break;
			case RULE_A: A[a_length++] = items[i]; break;
			}

		if (AB[0] == RULE_END) break;
		if (a_length == 0) continue;
		if (b_length == 0) continue;

		f(fp, m, M, epsilon, A, B, a_length, b_length);
	}

	free(AB);
#undef RULE_A
#undef RULE_B
#undef RULE_END
}

static void mine(struct fptree *fp, struct item_count *ic,
		double c, double epsilon, int num_triangles,
		int m, int M, struct drand48_data *buffer)
{
	int c_triangle, num_items, *items;
	double c_epsilon;

	for (c_triangle = 0; c_triangle < num_triangles; c_triangle++) {
		/* TODO: split based on triangles */
		c_epsilon = epsilon / num_triangles;

		num_items = num_allowed_items(fp, ic, m, M);
#if PRINT_TRIANGLE_CONTENT == 1
		printf("Triangle %d %d %d: %d", c_triangle, m, M, num_items);
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

		/* TODO: extract all rules with current items */
		for_all_rules(fp, items, num_items, m, M, c_epsilon,
				compute_cdf);

		/* TODO: sample */

		free(items);

end_loop:
#if PRINT_TRIANGLE_CONTENT == 1
		printf("\n");
#endif
		get_next_triangle(fp, ic, c, &m, &M);
	}
}

void dp2d(struct fptree *fp, double c, double eps, double eps_share, int ni,
		int minth)
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
	mine(fp, ic, c, eps - eps_share, nt, theta, fp->t, &randbuffer);

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
	return displacement(x, Y, X, X) * displacement(y, Y, X, Y);
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
