#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include <mpfr.h>

#include "dp2d.h"
#include "fp.h"
#include "globals.h"

/** settings for an experiment */
struct experiment {
	/* center of exp-mech invocation */
	int x, y;
	/* min-conf used */
	float c;
	/* exponent */
	float k;
	/* Z values */
	mpfr_t Z0, Z;
};

/**
 * Constant values used to generate experiments and to run them.
 */
#if 0 /* production */
float ks[] = {1, 1.2, 1.4, 1.6, 1.8, 2, 2.0, 2.2, 2.4, 2.6, 2.8, 3};
int ks_sz = sizeof(ks) / sizeof(ks[0]);
float cs[] = {.1, .2, .25, .3, .4, .5, .6, .7, .75, .8, .9};
int cs_sz = sizeof(cs) / sizeof(cs[0]);
#else /* test */
float ks[] = {1, 1.2};
int ks_sz = sizeof(ks) / sizeof(ks[0]);
float cs[] = {.7, .8, .9};
int cs_sz = sizeof(cs) / sizeof(cs[0]);
#endif

#define ROUND_MODE MPFR_RNDN

#define EXTRA 5
#define PPC 5
static struct experiment *prepare_experiments(int t, int minsup, int *exps_sz)
{
	int sz = ks_sz * (cs_sz * PPC + EXTRA), eix = 0, i, k;
	struct experiment *exps = calloc(sz, sizeof(exps[0]));

	/* extra experiments first */
	for (i = 0; i < ks_sz; i++) {
		/* upper corner */
		exps[eix].x = t; exps[eix].y = t;
		exps[eix].c = 0; exps[eix].k = ks[i];
		mpfr_init_set_ui(exps[eix].Z0, 0, ROUND_MODE);
		mpfr_init_set_ui(exps[eix].Z, 0, ROUND_MODE);
		eix++;
		/* lower corner */
		exps[eix].x = t; exps[eix].y = 0;
		exps[eix].c = 0; exps[eix].k = ks[i];
		mpfr_init_set_ui(exps[eix].Z0, 0, ROUND_MODE);
		mpfr_init_set_ui(exps[eix].Z, 0, ROUND_MODE);
		eix++;
		/* line bisector meeting point */
		exps[eix].x = t/2; exps[eix].y = t/2;
		exps[eix].c = 0; exps[eix].k = ks[i];
		mpfr_init_set_ui(exps[eix].Z0, 0, ROUND_MODE);
		mpfr_init_set_ui(exps[eix].Z, 0, ROUND_MODE);
		eix++;
		/* median meeting point */
		exps[eix].x = 2*t/3; exps[eix].y = t/3;
		exps[eix].c = 0; exps[eix].k = ks[i];
		mpfr_init_set_ui(exps[eix].Z0, 0, ROUND_MODE);
		mpfr_init_set_ui(exps[eix].Z, 0, ROUND_MODE);
		eix++;
		/* in between (t,t) and (minsup, minsup) */
		exps[eix].x = (t+minsup)/2; exps[eix].y = (t+minsup)/2;
		exps[eix].c = 0; exps[eix].k = ks[i];
		mpfr_init_set_ui(exps[eix].Z0, 0, ROUND_MODE);
		mpfr_init_set_ui(exps[eix].Z, 0, ROUND_MODE);
		eix++;
	}

	/* sanity check */
	if (eix != EXTRA * ks_sz)
		die("Error in code, %d %d\n", eix, EXTRA * ks_sz);

	/* confidence-based experiments */
	for (i = 0; i < ks_sz; i++)
		for (k = 0; k < cs_sz; k++) {
			/* max */
			exps[eix].x = t; exps[eix].y = t * cs[k];
			exps[eix].c = cs[k]; exps[eix].k = ks[i];
			mpfr_init_set_ui(exps[eix].Z0, 0, ROUND_MODE);
			mpfr_init_set_ui(exps[eix].Z, 0, ROUND_MODE);
			eix++;
			/* mid */
			exps[eix].x = t/2; exps[eix].y = t * cs[k] / 2;
			exps[eix].c = cs[k]; exps[eix].k = ks[i];
			mpfr_init_set_ui(exps[eix].Z0, 0, ROUND_MODE);
			mpfr_init_set_ui(exps[eix].Z, 0, ROUND_MODE);
			eix++;
			/* in between (t,t) and (minsup, minsup) */
			exps[eix].x = (t+minsup)/2; exps[eix].y = (t+minsup) * cs[k]/2;
			exps[eix].c = cs[k]; exps[eix].k = ks[i];
			mpfr_init_set_ui(exps[eix].Z0, 0, ROUND_MODE);
			mpfr_init_set_ui(exps[eix].Z, 0, ROUND_MODE);
			eix++;
			/* y=minsup */
			exps[eix].x = minsup/cs[k]; exps[eix].y = minsup;
			exps[eix].c = cs[k]; exps[eix].k = ks[i];
			mpfr_init_set_ui(exps[eix].Z0, 0, ROUND_MODE);
			mpfr_init_set_ui(exps[eix].Z, 0, ROUND_MODE);
			eix++;
			/* mid in valid region */
			exps[eix].x = (t + minsup / cs[k])/2; exps[eix].y = (t*cs[k] + minsup)/2;
			exps[eix].c = cs[k]; exps[eix].k = ks[i];
			mpfr_init_set_ui(exps[eix].Z0, 0, ROUND_MODE);
			mpfr_init_set_ui(exps[eix].Z, 0, ROUND_MODE);
			eix++;
		}

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
	printf("(%d, %d) %4.2f%%, ^%4.2f, | Z0=", xp->x, xp->y, xp->c, xp->k);
	mpfr_out_str(stdout, 10, 0, xp->Z0, ROUND_MODE);
	printf(" Z=");
	mpfr_out_str(stdout, 10, 0, xp->Z, ROUND_MODE);
	printf("\n");
}

void dp2d(struct fptree *fp, int minsup)
{
	struct experiment *exps;
	int exps_sz;

	printf("Preparing experiments ... ");
	fflush(stdout);
	exps = prepare_experiments(fp->t, minsup, &exps_sz);
	printf("OK (%d experiments)\n", exps_sz);

	run_on_experiments(exps, exps_sz, print_experiment);

	free_experiments(exps, exps_sz);
}
