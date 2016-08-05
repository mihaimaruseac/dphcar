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
	/* mode: private or nonprivate */
	enum {NONPRIVATE=0, PRIVATE} priv_mode;
	/* filename containing the transactions */
	char *tfname;
	/* global value for epsilon */
	double eps;
	/* fraction of epsilon for first step */
	double eps_ratio1;
	/* confidence threshold */
	double c0;
	/* max number of items in rule */
	size_t lmax;
	/* number of rules to extract */
	size_t k;
	/* random seed */
	long int seed;
} args;

static void usage(const char *prg)
{
	fprintf(stderr, "Usage: %s PRIV_MODE(n|p) TFILE EPS EPS_RATIO_1 C0 RLEN K [SEED]\n", prg);
	fprintf(stderr, "\tWhen PRIV_MODE is n (non-private) EPS*, K and SEED are ignored\n");
	exit(EXIT_FAILURE);
}

static void parse_arguments(int argc, char **argv)
{
	int i;

	printf("Called with: argc=%d\n", argc);
	for (i = 0; i < argc; i++)
		printf("%s ", argv[i]);
	printf("\n");

	if (argc < 8 || argc > 9)
		usage(argv[0]);
	if (!strncmp(argv[1], "n", 1)) args.priv_mode = NONPRIVATE;
	else if (!strncmp(argv[1], "p", 1)) args.priv_mode = PRIVATE;
	else usage(argv[0]);
	args.tfname = strdup(argv[2]);
	if (sscanf(argv[3], "%lf", &args.eps) != 1 || args.eps < 0)
		usage(argv[0]);
	if (sscanf(argv[4], "%lf", &args.eps_ratio1) != 1 || args.eps_ratio1 < 0 || args.eps_ratio1 >= 1)
		usage(argv[0]);
	if (sscanf(argv[5], "%lf", &args.c0) != 1 || args.c0 < 0 || args.c0 >= 1)
		usage(argv[0]);
	if (sscanf(argv[6], "%lu", &args.lmax) != 1 || args.lmax < 2 || args.lmax > 7)
		usage(argv[0]);
	if (sscanf(argv[7], "%lu", &args.k) != 1)
		usage(argv[0]);
	if (argc == 9) {
		if (sscanf(argv[8], "%ld", &args.seed) != 1)
			usage(argv[0]);
	} else
		args.seed = 42;
}

int main(int argc, char **argv)
{
	struct fptree fp;

	parse_arguments(argc, argv);

	fpt_read_from_file(args.tfname, &fp);
	printf("fp-tree: items: %lu, transactions: %lu, nodes: %d, depth: %d\n",
			fp.n, fp.t, fpt_nodes(&fp), fpt_height(&fp));

	dp2d(&fp, args.eps, args.eps_ratio1, args.c0, args.lmax, args.k,
			args.seed, args.priv_mode);

	fpt_cleanup(&fp);
	free(args.tfname);

	return 0;
}
