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
struct {
	/* filename containing the transactions */
	char *tfname;
	/* minimum confidence level */
	double c;
	/* global value for epsilon */
	double eps;
	/* fraction of epsilon for first step */
	double eps_share;
#if 0
	/* minimum threshold */
	int minth;
	/* k value (top k rules) */
	int k;
#endif
} args;

static void usage(char *prg)
{
	fprintf(stderr, "Usage: %s TFILE MINCONF EPSILON EPSILON_SHARE\n", prg);
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
}

static void parse_arguments(int argc, char **argv)
{
	/* TODO: use optparse */
	if (argc != 5)
		usage(argv[0]);
	args.tfname = strdup(argv[1]);
	if (sscanf(argv[2], "%lf", &args.c) != 1 || args.c < 0 || args.c >= 1)
		usage(argv[0]);
	if (sscanf(argv[3], "%lf", &args.eps) != 1)
		usage(argv[0]);
	if (sscanf(argv[4], "%lf", &args.eps_share) != 1 || args.eps_share < 0 || args.eps_share >= 1)
		usage(argv[0]);
#if 0
	if (sscanf(argv[2], "%d", &args.minth) != 1)
		usage(argv[0]);
	if (sscanf(argv[3], "%d", &args.k) != 1)
		usage(argv[0]);
#endif
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

	dp2d(&fp, args.c, args.eps, args.eps_share);

	fpt_cleanup(&fp);

	return 0;
}
