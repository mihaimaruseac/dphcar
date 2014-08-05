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
#if 0
	/* number of items in the first set */
	int ni;
	/* minimum threshold */
	int minth;
	/* threshold for S set */
	int thS;
	/* threshold for L set */
	int thL;
	/* weight for M set */
	double wM;
	/* weight for L set */
	double wL;
	/* high item count */
	int hic;
#endif
} args;

static void usage(const char *prg)
{
	fprintf(stderr, "Usage: %s TFILE\n", prg);
	exit(EXIT_FAILURE);
}

static void parse_arguments(int argc, char **argv)
{
	int i;

	printf("Called with: argc=%d\n", argc);
	for (i = 0; i < argc; i++)
		printf("%s ", argv[i]);
	printf("\n");

	if (argc != 4)
		usage(argv[0]);
	args.tfname = strdup(argv[1]);
	if (sscanf(argv[2], "%lf", &args.eps) != 1)
		usage(argv[0]);
	if (sscanf(argv[3], "%lf", &args.eps_share) != 1 || args.eps_share < 0 || args.eps_share >= 1)
		usage(argv[0]);
#if 0
	if (sscanf(argv[6], "%d", &args.ni) != 1)
		usage(argv[0]);
	if (sscanf(argv[7], "%d", &args.minth) != 1)
		usage(argv[0]);
	if (sscanf(argv[8], "%d", &args.thS) != 1)
		usage(argv[0]);
	if (sscanf(argv[9], "%d", &args.thL) != 1)
		usage(argv[0]);
	if (sscanf(argv[10], "%lf", &args.wM) != 1)
		usage(argv[0]);
	if (sscanf(argv[11], "%lf", &args.wL) != 1)
		usage(argv[0]);
	if (sscanf(argv[12], "%d", &args.hic) != 1)
		usage(argv[0]);
#endif
}

int main(int argc, char **argv)
{
	struct fptree fp;

	parse_arguments(argc, argv);

	fpt_read_from_file(args.tfname, &fp);
	printf("fp-tree: items: %lu, transactions: %lu, nodes: %d, depth: %d\n",
			fp.n, fp.t, fpt_nodes(&fp), fpt_height(&fp));

	dp2d(&fp, /*args.c,*/ args.eps, args.eps_share/*, args.ni, args.minth,*/
			/*args.ifname, args.hic*/);

	fpt_cleanup(&fp);
	free(args.tfname);

	return 0;
}
