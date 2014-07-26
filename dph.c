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
	/* minimum confidence level */
	double c;
	/* global value for epsilon */
	double eps;
	/* fraction of epsilon for first step */
	double eps_share;
	/* number of items in the first set */
	int ni;
	/* minimum threshold */
	int minth;
#if 0
	/* k value (top k rules) */
	int k;
#endif
} args;

static void usage(const char *prg)
{
	fprintf(stderr, "Usage: %s TFILE MINCONF EPSILON EPSILON_SHARE NITEMS"
			" MINTHRESHOLD\n", prg);
	exit(EXIT_FAILURE);
}

static void parse_arguments(int argc, char **argv)
{
	/* TODO: use optparse */
	if (argc != 7)
		usage(argv[0]);
	args.tfname = strdup(argv[1]);
	if (sscanf(argv[2], "%lf", &args.c) != 1 || args.c < 0 || args.c >= 1)
		usage(argv[0]);
	if (sscanf(argv[3], "%lf", &args.eps) != 1)
		usage(argv[0]);
	if (sscanf(argv[4], "%lf", &args.eps_share) != 1 || args.eps_share < 0 || args.eps_share >= 1)
		usage(argv[0]);
	if (sscanf(argv[5], "%d", &args.ni) != 1)
		usage(argv[0]);
	if (sscanf(argv[6], "%d", &args.minth) != 1)
		usage(argv[0]);
}

static void final_cleanup(void)
{
	free(args.tfname);
}

int main(int argc, char **argv)
{
	struct fptree fp;

	atexit(final_cleanup);
	parse_arguments(argc, argv);

	fpt_read_from_file(args.tfname, &fp);
	printf("fp-tree: items: %d, transactions: %d, nodes: %d, depth: %d\n",
			fp.n, fp.t, fpt_nodes(&fp), fpt_height(&fp));

	dp2d(&fp, args.c, args.eps, args.eps_share, args.ni, args.minth);

	fpt_cleanup(&fp);

	return 0;
}
