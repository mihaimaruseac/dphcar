/**
 * Differentially-private high-confidence association rule extractor.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dp2d.h"
#include "fp.h"
#include "itstree.h"

/* Command line arguments */
static struct {
	/* filename containing the transactions */
	char *tfname;
	/* filename containing the recall */
	char *rfname;
	/* global value for epsilon */
	double eps;
	/* fraction of epsilon for first step */
	double er1;
	/* confidence threshold */
	double c0;
	/* max number of items in rule */
	size_t lmax;
	/* num items (to be removed later) */
	size_t ni;
	/* branching factor */
	size_t cspl;
	/* random seed */
	long int seed;
} args;

static void usage(const char *prg)
{
	fprintf(stderr, "Usage: %s TFILE IFILE EPS EPS_RATIO_1 C0 RLEN NI BF [SEED]\n", prg);
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
	args.rfname = strdup(argv[2]);
	if (sscanf(argv[3], "%lf", &args.eps) != 1 || args.eps < 0)
		usage(argv[0]);
	if (sscanf(argv[4], "%lf", &args.er1) != 1 || args.er1 < 0 || args.er1 >= 1)
		usage(argv[0]);
	if (sscanf(argv[5], "%lf", &args.c0) != 1 || args.c0 < 0 || args.c0 >= 1)
		usage(argv[0]);
	if (sscanf(argv[6], "%lu", &args.lmax) != 1 || args.lmax < 2 || args.lmax > 7)
		usage(argv[0]);
	if (sscanf(argv[7], "%lu", &args.ni) != 1)
		usage(argv[0]);
	if (sscanf(argv[8], "%lu", &args.cspl) != 1)
		usage(argv[0]);
	if (argc == 10) {
		if (sscanf(argv[9], "%ld", &args.seed) != 1)
			usage(argv[0]);
	} else
		args.seed = 42;
}

int main(int argc, char **argv)
{
	struct itstree_node *itst;
	struct fptree fp;

	parse_arguments(argc, argv);

	fpt_read_from_file(args.tfname, &fp);
	printf("fp-tree: items: %lu, transactions: %lu, nodes: %d, depth: %d\n",
			fp.n, fp.t, fpt_nodes(&fp), fpt_height(&fp));

	if (!strncmp(args.rfname, "-", 1))
		itst = init_empty_itstree();
	else
		itst = load_its(args.rfname, args.lmax, args.ni);
	dp2d(&fp, itst, args.eps, args.er1, args.c0, args.lmax,
			args.ni, args.cspl, args.seed);

	free_itstree(itst);
	fpt_cleanup(&fp);
	free(args.tfname);
	free(args.rfname);

	return 0;
}
