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
	/* filename containing the item scores */
	char *ifname;
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
	fprintf(stderr, "Usage: %s TFILE IFILE MINCONF EPSILON EPSILON_SHARE "
			"NITEMS MINTHRESHOLD\n", prg);
	exit(EXIT_FAILURE);
}

static void parse_arguments(int argc, char **argv)
{
    int i;
    printf("Called with: argc=%d\n", argc);
    for (i = 0; i < argc; i++)
        printf("%s ", argv[i]);
    printf("\n");
	/* TODO: use optparse */
	if (argc != 8)
		usage(argv[0]);
	args.tfname = strdup(argv[1]);
	args.ifname = strdup(argv[2]);
	if (sscanf(argv[3], "%lf", &args.c) != 1 || args.c < 0 || args.c >= 1)
		usage(argv[0]);
	if (sscanf(argv[4], "%lf", &args.eps) != 1)
		usage(argv[0]);
	if (sscanf(argv[5], "%lf", &args.eps_share) != 1 || args.eps_share < 0 || args.eps_share >= 1)
		usage(argv[0]);
	if (sscanf(argv[6], "%d", &args.ni) != 1)
		usage(argv[0]);
	if (sscanf(argv[7], "%d", &args.minth) != 1)
		usage(argv[0]);
}

int main(int argc, char **argv)
{
	struct fptree fp;

	parse_arguments(argc, argv);

	fpt_read_from_file(args.tfname, args.ifname, &fp);
	printf("fp-tree: items: %d, transactions: %d, nodes: %d, depth: %d\n",
			fp.n, fp.t, fpt_nodes(&fp), fpt_height(&fp));

	dp2d(&fp, args.c, args.eps, args.eps_share, args.ni, args.minth);

	fpt_cleanup(&fp);
	free(args.tfname);
	free(args.ifname);

	return 0;
}
