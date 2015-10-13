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
	/* declared with returns? */
	int has_returns;
	/* global value for epsilon */
	double eps;
	/* fraction of epsilon for first step */
	double eps_share;
	/* maximum rule length */
	size_t l_max_r;
	/* min conf value (c0) */
	double c0;
	/* support threshold */
	double smin;
	/* rule threshold */
	double cmin;
	/* number of rules to extract */
	size_t k;
	/* random seed */
	long int seed;
} args;

static void usage(const char *prg)
{
	fprintf(stderr, "Usage: %s TFILE r/n EPS EPS_SHARE K L_MAX_R C0 SMIN CMIN [SEED]\n", prg);
	exit(EXIT_FAILURE);
}

static void parse_arguments(int argc, char **argv)
{
	int i;

	printf("Called with: argc=%d\n", argc);
	for (i = 0; i < argc; i++)
		printf("%s ", argv[i]);
	printf("\n");

	if (argc < 10 || argc > 11)
		usage(argv[0]);
	args.tfname = strdup(argv[1]);
	if (!strncmp(argv[2], "r", 1)) args.has_returns = 1;
	else if (!strncmp(argv[2], "n", 1)) args.has_returns = 0;
	else usage(argv[0]);
	if (sscanf(argv[3], "%lf", &args.eps) != 1 || args.eps < 0)
		usage(argv[0]);
	if (sscanf(argv[4], "%lf", &args.eps_share) != 1 || args.eps_share < 0 || args.eps_share >= 1)
		usage(argv[0]);
	if (sscanf(argv[5], "%lu", &args.k) != 1)
		usage(argv[0]);
	if (sscanf(argv[6], "%lu", &args.l_max_r) != 1 || args.l_max_r < 2 || args.l_max_r > 7)
		usage(argv[0]);
	if (sscanf(argv[7], "%lf", &args.c0) != 1 || args.c0 < 0 || args.c0 >= 1)
		usage(argv[0]);
	if (sscanf(argv[8], "%lf", &args.smin) != 1 || args.smin < 0 || args.smin >= 1)
		usage(argv[0]);
	if (sscanf(argv[9], "%lf", &args.cmin) != 1 || args.cmin < 0 || args.cmin >= 1)
		usage(argv[0]);
	if (argc == 11) {
		if (sscanf(argv[10], "%ld", &args.seed) != 1)
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
			"l_max_t: %lu l_max_r: %lu returns: %d/%d\n",
			fp.n, fp.e, fp.t, fp.l_max_t, fp.l_max_r,
			fp.has_returns, args.has_returns);

	if (fp.has_returns && !args.has_returns)
		die("File has returns in transaction but declared as with no returns");
	fp.has_returns = args.has_returns;

	dp2d(&fp, args.eps, args.eps_share, args.k, args.c0,
			args.smin * fp.t, args.cmin * fp.t, args.seed);

	fpt_cleanup(&fp);
	free(args.tfname);

	return 0;
}
