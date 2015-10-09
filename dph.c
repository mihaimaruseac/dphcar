/**
 * Differentially-private high-confidence association rule extractor.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <gmp.h>
#include <mpfr.h>

#include "dp2d.h"
#include "fp.h"
#include "globals.h"

/* Command line arguments */
static struct {
	/* filename containing the transactions */
	char *tfname;
	/* global value for epsilon */
	double eps;
	/* fraction of epsilon for first step */
	double eps_share;
	/* maximum rule length */
	size_t l_max_r;
	/* min alpha value */
	double minalpha;
	/* max items in generation step */
	size_t mis;
	/* number of trees (0 == no tree) */
	size_t nt;
	/* number of rules to extract */
	size_t k;
	/* random seed */
	long int seed;
} args;

static void usage(const char *prg)
{
	fprintf(stderr, "Usage: %s TFILE EPS EPS_SHARE L_MAX_R MINALPHA MIS NT K [SEED]\n", prg);
	exit(EXIT_FAILURE);
}

static void parse_arguments(int argc, char **argv)
{
	int i;

	printf("Called with: argc=%d\n", argc);
	for (i = 0; i < argc; i++)
		printf("%s ", argv[i]);
	printf("\n");

	if (argc < 9 || argc > 10)
		usage(argv[0]);
	args.tfname = strdup(argv[1]);

	if (sscanf(argv[2], "%lf", &args.eps) != 1 || args.eps < 0)
		usage(argv[0]);
	if (sscanf(argv[3], "%lf", &args.eps_share) != 1 || args.eps_share < 0 || args.eps_share >= 1)
		usage(argv[0]);
	if (sscanf(argv[4], "%lu", &args.l_max_r) != 1 || args.l_max_r < 2 || args.l_max_r > 7)
		usage(argv[0]);
	if (sscanf(argv[5], "%lf", &args.minalpha) != 1)
		usage(argv[0]);
	if (sscanf(argv[6], "%lu", &args.mis) != 1 || args.mis < 2 || args.mis > 7)
		usage(argv[0]);
	if (sscanf(argv[7], "%lu", &args.nt) != 1)
		usage(argv[0]);
	if (sscanf(argv[8], "%lu", &args.k) != 1)
		usage(argv[0]);
	if (argc == 10) {
		if (sscanf(argv[9], "%ld", &args.seed) != 1)
			usage(argv[0]);
	} else
		args.seed = 42;
}

int main(int argc, char **argv)
{
	struct fptree fp;

	parse_arguments(argc, argv);

	fpt_read_from_file(args.tfname, args.l_max_r, &fp);
	printf("data-struct: nodes: %lu, edges: %lu, transactions: %lu, "
			"l_max_t: %lu l_max_r: %lu returns: %d\n",
			fp.n, fp.e, fp.t, fp.l_max_t, fp.l_max_r,
			fp.has_returns);

	dp2d(&fp, args.eps, args.eps_share,
			args.mis, args.nt, args.k, args.minalpha,
			args.seed);

	fpt_cleanup(&fp);
	free(args.tfname);

	return 0;
}
