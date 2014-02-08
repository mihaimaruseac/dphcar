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
	/* minimum threshold */
	int minth;
} args;

static void usage(char *prg)
{
	fprintf(stderr, "Usage: %s TFILE MINSUP\n", prg);
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
}

static void parse_arguments(int argc, char **argv)
{
	/* TODO: use optparse */
	if (argc != 3)
		usage(argv[0]);
	args.tfname = strdup(argv[1]);
	if (sscanf(argv[2], "%d", &args.minth) != 1)
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

	dp2d(&fp, args.minth);

	fpt_cleanup(&fp);

	return 0;
}
