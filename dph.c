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
	/* minimum threshold */
	size_t minth;
	/* max spread */
	double mu;
	/* max item spread */
	size_t mis;
	/* number of rules to extract */
	size_t k;
} args;

static void usage(const char *prg)
{
	fprintf(stderr, "Usage: %s TFILE EPS EPS_SHARE MINTH MU MIS K\n", prg);
	exit(EXIT_FAILURE);
}

static void parse_arguments(int argc, char **argv)
{
	int i;

	printf("Called with: argc=%d\n", argc);
	for (i = 0; i < argc; i++)
		printf("%s ", argv[i]);
	printf("\n");

	if (argc != 8)
		usage(argv[0]);
	args.tfname = strdup(argv[1]);
	if (sscanf(argv[2], "%lf", &args.eps) != 1 || args.eps < 0)
		usage(argv[0]);
	if (sscanf(argv[3], "%lf", &args.eps_share) != 1 || args.eps_share < 0 || args.eps_share >= 1)
		usage(argv[0]);
	if (sscanf(argv[4], "%lu", &args.minth) != 1)
		usage(argv[0]);
	if (sscanf(argv[5], "%lf", &args.mu) != 1 || args.mu < 0 || args.mu >= 1)
		usage(argv[0]);
	if (sscanf(argv[6], "%lu", &args.mis) != 1 || args.mis < 2)
		usage(argv[0]);
	if (sscanf(argv[7], "%lu", &args.k) != 1)
		usage(argv[0]);
}

int main(int argc, char **argv)
{
	struct fptree fp;

	parse_arguments(argc, argv);

	fpt_read_from_file(args.tfname, &fp);
	printf("fp-tree: items: %lu, transactions: %lu, nodes: %d, depth: %d\n",
			fp.n, fp.t, fpt_nodes(&fp), fpt_height(&fp));

	dp2d(&fp, args.eps, args.eps_share, args.minth, args.mu, args.mis,
			args.k);

	fpt_cleanup(&fp);
	free(args.tfname);

	return 0;
}
