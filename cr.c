/**
 * Differentially-private high-confidence association rule extractor.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dp2d.h"
#include "fp.h"
#include "itstree.h"
#include "recall.h"

/* Command line arguments */
static struct {
	/* filename containing the transactions */
	char *tfname;
	/* max number of items in rule */
	size_t lmax;
	/* num items (to be removed later) */
	size_t ni;
} args;

static void usage(const char *prg)
{
	fprintf(stderr, "Usage: %s TFILE RMAX NI\n", prg);
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
	if (sscanf(argv[2], "%lu", &args.lmax) != 1 || args.lmax < 2 || args.lmax > 7)
		usage(argv[0]);
	if (sscanf(argv[3], "%lu", &args.ni) != 1)
		usage(argv[0]);
}

int main(int argc, char **argv)
{
	struct itstree_node *itst;
	struct fptree fp;

	parse_arguments(argc, argv);

	fpt_read_from_file(args.tfname, &fp);
	printf("fp-tree: items: %lu, transactions: %lu, nodes: %d, depth: %d\n",
			fp.n, fp.t, fpt_nodes(&fp), fpt_height(&fp));

	itst = build_recall_tree(&fp, args.lmax, args.ni);
	save_its(itst, args.tfname, args.lmax, args.ni);

	free_itstree(itst);
	fpt_cleanup(&fp);
	free(args.tfname);

	return 0;
}
