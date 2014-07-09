#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include <mpfr.h>

#include "dp2d.h"
#include "fp.h"
#include "globals.h"

struct item_count {
	int value;
	int real_count;
	double noisy_count;
};

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
	if (ib->noisy_count < ia->noisy_count)
		return -1;
	if (ib->noisy_count > ia->noisy_count)
		return 1;
	return 0;

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

static double get_first_theta_value(struct fptree *fp, struct item_count *ic,
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

	return x.noisy_count;
}

static double quality(int x, int y, double X, double Y);
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

static void mine(struct fptree *fp, struct item_count *ic,
		double c, double epsilon, double m, double M,
		struct drand48_data *buffer)
{
}

void dp2d(struct fptree *fp, double c, double eps, double eps_share, int ni)
{
	struct item_count *ic = alloc_items(fp->n);
	double epsilon_step1 = eps * eps_share;
	struct drand48_data randbuffer;
	double theta;

	init_rng(&randbuffer);

	printf("Running dp2D with ni=%d, c=%lf, eps=%lf, eps_share=%lf\n",
			ni, c, eps, eps_share);

	printf("Step 1: compute noisy counts for items: eps_1 = %lf\n", epsilon_step1);
	build_items_table(fp, ic, epsilon_step1, &randbuffer);

	printf("Step 2: get min support for %d items: ", ni);
	theta = get_first_theta_value(fp, ic, ni, c, c * fp->t);
	printf("%lf\n", theta);

	printf("Step 3: mining rules\n");
	mine(fp, ic, c, eps - eps_share, theta, fp->t, &randbuffer);


#if 0
	int i;
	for (i = 0; i < fp->n; i++)
		printf("%d %d %lf\n", ic[i].value, ic[i].real_count, ic[i].noisy_count);
#endif

	free_items(ic);
}
