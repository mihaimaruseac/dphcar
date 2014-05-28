#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include <mpfr.h>

#include "dp2d.h"
#include "fp.h"
#include "globals.h"

/** settings for an experiment */
struct experiment {
	/* center of exp-mechanism invocation */
	int x, y;
	/* min-conf used */
	float c;
	/* exponent */
	float exponent;
	/* epsilon */
	float epsilon;
	/* Z values */
	mpfr_t Z0, Z;
};

/**
 * Constant values used to generate experiments and to run them.
 */
float exponents[] = {1, 1.2, 1.4, 1.6, 1.8, 2, 2.0, 2.2, 2.4, 2.6, 2.8, 3};
int exponents_sz = sizeof(exponents) / sizeof(exponents[0]);
float cs[] = {.1, .2, .25, .3, .4, .5, .6, .7, .75, .8, .9};
int cs_sz = sizeof(cs) / sizeof(cs[0]);
float epsilons[] = {.1, .2, .3, .4, .5, .6, .7, .8, .9, 1};
int epsilons_sz = sizeof(epsilons) / sizeof(epsilons[0]);

#define ROUND_MODE MPFR_RNDN

static void setup_experiment(struct experiment *xp,
		int x, int y, float c, float exponent, float epsilon)
{
	xp->x = x;
	xp->y = y;
	xp->c = c;
	xp->exponent = exponent;
	xp->epsilon = epsilon;
	mpfr_init_set_ui(xp->Z0, 0, ROUND_MODE);
	mpfr_init_set_ui(xp->Z, 0, ROUND_MODE);
}

#define EXTRA 5
#define PPC 5
static struct experiment *prepare_experiments(int t, int minsup, int k, int *exps_sz)
{
	int sz = exponents_sz * epsilons_sz * (cs_sz * PPC + EXTRA), eix = 0, i, j, l;
	struct experiment *exps = calloc(sz, sizeof(exps[0]));

	/* extra experiments first */
	for (i = 0; i < exponents_sz; i++)
		for (l = 0; l < epsilons_sz; l++) {
			/* upper corner */
			setup_experiment(&exps[eix++], t, t, 0, exponents[i], epsilons[l] / k);
			/* lower corner */
			setup_experiment(&exps[eix++], t, 0, 0, exponents[i], epsilons[l] / k);
			/* line bisector meeting point */
			setup_experiment(&exps[eix++], t/2, t/2, 0, exponents[i], epsilons[l] / k);
			/* median meeting point */
			setup_experiment(&exps[eix++], 2*t/3, t/3, 0, exponents[i], epsilons[l] / k);
			/* in between (t,t) and (minsup, minsup) */
			setup_experiment(&exps[eix++], (t+minsup)/2, (t+minsup)/2, 0, exponents[i], epsilons[l] / k);
		}

	/* sanity check */
	if (eix != EXTRA * exponents_sz * epsilons_sz)
		die("Error in code, %d %d\n", eix, EXTRA * exponents_sz * epsilons_sz);

	/* confidence-based experiments */
	for (i = 0; i < exponents_sz; i++)
		for (j = 0; j < cs_sz; j++)
			for (l = 0; l < epsilons_sz; l++) {
			/* max */
			setup_experiment(&exps[eix++], t, t * cs[j], cs[j], exponents[i], epsilons[l] / k);
			/* mid */
			setup_experiment(&exps[eix++], t/2, t * cs[j]/2, cs[j], exponents[i], epsilons[l] / k);
			/* in between (t,t) and (minsup, minsup) */
			setup_experiment(&exps[eix++], (t + minsup)/2, (t + minsup) * cs[j]/2, cs[j], exponents[i], epsilons[l] / k);
			/* y=minsup */
			setup_experiment(&exps[eix++], minsup/cs[j], minsup, cs[j], exponents[i], epsilons[l] / k);
			/* mid in valid region */
			setup_experiment(&exps[eix++], (t + minsup/cs[j])/2, (t*cs[j] + minsup)/2, cs[j], exponents[i], epsilons[l] / k);
		}

	/* sanity check */
	if (eix != sz)
		die("Error in code, %d %d\n", eix, sz);

	*exps_sz = eix;
	return exps;
}
#undef EXTRA
#undef PPC

static void run_on_experiments(struct experiment *exps, int exps_sz, void (*expf)(struct experiment *))
{
	int i;

	for (i = 0; i < exps_sz; i++)
		expf(&exps[i]);
}

static void free_experiment(struct experiment *xp)
{
	mpfr_clear(xp->Z);
	mpfr_clear(xp->Z0);
}

static void free_experiments(struct experiment *exps, int exps_sz)
{
	run_on_experiments(exps, exps_sz, free_experiment);
	free(exps);
}

static void print_experiment(struct experiment *xp)
{
	printf("(%d, %d) %4.2f%%, ^%4.2f, | Z0=", xp->x, xp->y, xp->c, xp->exponent);
	mpfr_out_str(stdout, 10, 0, xp->Z0, ROUND_MODE);
	printf(" Z=");
	mpfr_out_str(stdout, 10, 0, xp->Z, ROUND_MODE);
	printf("\n");
}

void dp2d(struct fptree *fp, int minsup, int k)
{
	struct experiment *exps;
	int exps_sz;

	printf("Preparing experiments ... ");
	fflush(stdout);
	exps = prepare_experiments(fp->t, minsup, k, &exps_sz);
	printf("OK (%d experiments)\n", exps_sz);

#if 0
	fpt_node_print(fp->tree);
	fpt_table_print(fp->table, fp->n);
#endif

	//run_on_experiments(exps, exps_sz, print_experiment);

	free_experiments(exps, exps_sz);
}
