/**
 * Non-private.
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <gmp.h>
#include <mpfr.h>

#include "globals.h"
#include "histogram.h"

#define LINELEN 1024

struct experiment {
	const size_t num_items;
	const int items[3];
	struct histogram *h;
	int valid;
};

#define RETAIL 1
#define KOSARAK 2
#define MUSHROOM 3
#define PUMSB 4
#define DATASET_OPT PUMSB
#if DATASET_OPT == RETAIL
#define DATASET "retail"
static struct experiment exps [] = {
	{.num_items = 1, .items = {39}},
	{.num_items = 1, .items = {48}},
	{.num_items = 1, .items = {38}},
	{.num_items = 1, .items = {32}},
	{.num_items = 1, .items = {41}},
	{.num_items = 1, .items = {65}},
	{.num_items = 1, .items = {89}},
	{.num_items = 1, .items = {225}},
	{.num_items = 1, .items = {237}},
	{.num_items = 1, .items = {170}},
	{.num_items = 2, .items = {48, 170}},
	{.num_items = 2, .items = {48, 65}},
	{.num_items = 2, .items = {32, 89}},
	{.num_items = 2, .items = {48, 32}},
	{.num_items = 2, .items = {38, 65}},
	{.num_items = 2, .items = {39, 38}},
	{.num_items = 2, .items = {39, 48}},
	{.num_items = 2, .items = {225, 170}},
	{.num_items = 2, .items = {41, 225}},
	{.num_items = 2, .items = {38, 170}},
	{.num_items = 3, .items = {32, 65, 225}},
	{.num_items = 3, .items = {38, 32, 170}},
	{.num_items = 3, .items = {39, 38, 65}},
	{.num_items = 3, .items = {32, 237, 170}},
	{.num_items = 3, .items = {32, 41, 65}},
	{.num_items = 3, .items = {39, 48, 38}},
	{.num_items = 3, .items = {38, 32, 65}},
	{.num_items = 3, .items = {39, 38, 170}},
	{.num_items = 3, .items = {48, 89, 237}},
	{.num_items = 3, .items = {48, 38, 237}},
};
#elif DATASET_OPT == KOSARAK
#define DATASET "kosarak"
static struct experiment exps [] = {
	{.num_items = 1, .items = {6}},
	{.num_items = 1, .items = {3}},
	{.num_items = 1, .items = {11}},
	{.num_items = 1, .items = {1}},
	{.num_items = 1, .items = {218}},
	{.num_items = 1, .items = {7}},
	{.num_items = 1, .items = {4}},
	{.num_items = 1, .items = {27}},
	{.num_items = 1, .items = {148}},
	{.num_items = 1, .items = {55}},
	{.num_items = 2, .items = {3, 55}},
	{.num_items = 2, .items = {3, 7}},
	{.num_items = 2, .items = {1, 4}},
	{.num_items = 2, .items = {3, 1}},
	{.num_items = 2, .items = {11, 7}},
	{.num_items = 2, .items = {6, 11}},
	{.num_items = 2, .items = {6, 3}},
	{.num_items = 2, .items = {27, 55}},
	{.num_items = 2, .items = {218, 27}},
	{.num_items = 2, .items = {11, 55}},
	{.num_items = 3, .items = {1, 7, 27}},
	{.num_items = 3, .items = {11, 1, 55}},
	{.num_items = 3, .items = {6, 11, 7}},
	{.num_items = 3, .items = {1, 148, 55}},
	{.num_items = 3, .items = {1, 218, 7}},
	{.num_items = 3, .items = {6, 3, 11}},
	{.num_items = 3, .items = {11, 1, 7}},
	{.num_items = 3, .items = {6, 11, 55}},
	{.num_items = 3, .items = {3, 4, 148}},
	{.num_items = 3, .items = {3, 11, 148}},
};
#elif DATASET_OPT == MUSHROOM
#define DATASET "mushroom"
static struct experiment exps [] = {
	{.num_items = 1, .items = {85}},
	{.num_items = 1, .items = {86}},
	{.num_items = 1, .items = {34}},
	{.num_items = 1, .items = {90}},
	{.num_items = 1, .items = {36}},
	{.num_items = 1, .items = {39}},
	{.num_items = 1, .items = {59}},
	{.num_items = 1, .items = {63}},
	{.num_items = 1, .items = {24}},
	{.num_items = 1, .items = {53}},
	{.num_items = 2, .items = {86, 53}},
	{.num_items = 2, .items = {86, 39}},
	{.num_items = 2, .items = {90, 59}},
	{.num_items = 2, .items = {86, 90}},
	{.num_items = 2, .items = {34, 39}},
	{.num_items = 2, .items = {85, 34}},
	{.num_items = 2, .items = {85, 86}},
	{.num_items = 2, .items = {63, 53}},
	{.num_items = 2, .items = {36, 63}},
	{.num_items = 2, .items = {34, 53}},
	{.num_items = 3, .items = {90, 39, 63}},
	{.num_items = 3, .items = {34, 90, 53}},
	{.num_items = 3, .items = {85, 34, 39}},
	{.num_items = 3, .items = {90, 24, 53}},
	{.num_items = 3, .items = {90, 36, 39}},
	{.num_items = 3, .items = {85, 86, 34}},
	{.num_items = 3, .items = {34, 90, 39}},
	{.num_items = 3, .items = {85, 34, 53}},
	{.num_items = 3, .items = {86, 59, 24}},
	{.num_items = 3, .items = {86, 34, 24}},
};
#elif DATASET_OPT == PUMSB
#define DATASET "pumsb"
static struct experiment exps [] = {
	{.num_items = 1, .items = {7072}},
	{.num_items = 1, .items = {161}},
	{.num_items = 1, .items = {197}},
	{.num_items = 1, .items = {84}},
	{.num_items = 1, .items = {4499}},
	{.num_items = 1, .items = {4502}},
	{.num_items = 1, .items = {168}},
	{.num_items = 1, .items = {4933}},
	{.num_items = 1, .items = {4937}},
	{.num_items = 1, .items = {4496}},
	{.num_items = 2, .items = {161, 4496}},
	{.num_items = 2, .items = {161, 4502}},
	{.num_items = 2, .items = {84, 168}},
	{.num_items = 2, .items = {161, 84}},
	{.num_items = 2, .items = {197, 4502}},
	{.num_items = 2, .items = {7072, 197}},
	{.num_items = 2, .items = {7072, 161}},
	{.num_items = 2, .items = {4933, 4496}},
	{.num_items = 2, .items = {4499, 4933}},
	{.num_items = 2, .items = {197, 4496}},
	{.num_items = 3, .items = {84, 4502, 4933}},
	{.num_items = 3, .items = {197, 84, 4496}},
	{.num_items = 3, .items = {7072, 197, 4502}},
	{.num_items = 3, .items = {84, 4937, 4496}},
	{.num_items = 3, .items = {84, 4499, 4502}},
	{.num_items = 3, .items = {7072, 161, 197}},
	{.num_items = 3, .items = {197, 84, 4502}},
	{.num_items = 3, .items = {7072, 197, 4496}},
	{.num_items = 3, .items = {161, 168, 4937}},
	{.num_items = 3, .items = {161, 197, 4937}},
};
#endif


/* Command line arguments */
static struct {
	/* filename containing the transactions */
	char *tfname;
	/* dirname for the outputs */
	char *dirname;
} args;

static void usage(const char *prg)
{
	fprintf(stderr, "Usage: %s iFILE OUTDIRNAME\n", prg);
	exit(EXIT_FAILURE);
}

static void parse_arguments(int argc, char **argv)
{
	int i;

	printf("Called with: argc=%d\n", argc);
	for (i = 0; i < argc; i++)
		printf("%s ", argv[i]);
	printf("\n");

	if (argc != 3)
		usage(argv[0]);
	args.tfname = strdup(argv[1]);

	if (!strstr(args.tfname, DATASET))
		die("Code was not compiled for this dataset!");

	args.dirname = strdup(argv[2]);
}

struct itemset {
	unsigned int size;
	int *elems;
	int count;
};

struct file_data {
	int transaction_count;
	struct itemset *itemsets;
	unsigned int item_count;
};

static int itemset_cmp(const void *isa, const void *isb)
{
	const struct itemset *a = isa;
	const struct itemset *b = isb;
	unsigned int i;

	if (a->size != b->size)
		return a->size - b->size;

	for (i = 0; i < a->size; i++)
		if (a->elems[i] != b->elems[i])
			return a->elems[i] - b->elems[i];

	return 0;
}

static void free_itemsets(struct file_data *data)
{
	if (data->itemsets) {
		unsigned int i;

		for (i = 0; i < data->item_count; i++)
			free(data->itemsets[i].elems);
		free(data->itemsets);
	}

	free(data);
}

static struct file_data *read_input(char *fname)
{
	struct file_data *data = calloc(1, sizeof(*data));
	FILE *f = fopen(fname, "r");
	char buffer[LINELEN], *p;
	unsigned int i, j;

	printf("Reading the input...\n");

	if (!fgets(buffer, LINELEN, f) ||
			sscanf(buffer, "(%d)", &(data->transaction_count)) != 1)
		goto err;

	data->item_count = data->transaction_count;
	data->itemsets = calloc(data->item_count, sizeof(*data->itemsets));

	for (i = 0;; i++) {
		if (i == data->item_count) {
			data->item_count *= 2;
			data->itemsets = realloc(data->itemsets, data->item_count * sizeof(*data->itemsets));
		}

		if (!fgets(buffer, LINELEN, f))
			break;

		p = strchr(buffer, ' ');
		data->itemsets[i].size = 0;
		while (p) {
			p++;
			data->itemsets[i].size++;
			p = strchr(p, ' ');
		}

		data->itemsets[i].elems = calloc(data->itemsets[i].size, sizeof(data->itemsets[i].elems[0]));

		p = buffer;
		for (j = 0; j < data->itemsets[i].size; j++) {
			if (sscanf(p, "%d", &data->itemsets[i].elems[j]) != 1)
				goto err;
			p = strchr(p, ' ') + 1;
		}

		qsort(data->itemsets[i].elems, data->itemsets[i].size, sizeof(data->itemsets[i].elems[0]), int_cmp);

		if (sscanf(p, "(%d)", &data->itemsets[i].count) != 1)
			goto err;
	}

	data->item_count = i;
	qsort(data->itemsets, data->item_count, sizeof(data->itemsets[0]), itemset_cmp);
	printf("Reading done: %d itemsets in %d transactions.\n", data->item_count, data->transaction_count);
out:
	fclose(f);
	return data;
err:
	printf("Invalid data in input file |%s|\n", buffer);
	data->item_count = i;
	free_itemsets(data);
	data = NULL;
	goto out;
}

int main(int argc, char **argv)
{
	size_t num_exps = sizeof(exps) / sizeof(exps[0]);
	size_t i, j, k, l, ve, i1, i2;
	struct file_data *data = NULL;
	int status = EXIT_FAILURE;
	struct itemset tmp, *res;
	int *tmp_items = NULL;
	int *masks = NULL;
	double c;
	FILE *f;

	parse_arguments(argc, argv);

	data = read_input(argv[1]);
	if (!data)
		goto out;

	for (i = 0; i < num_exps; i++) {
		exps[i].h = init_histogram();
		exps[i].valid = 0;
	}

	tmp.size = data->itemsets[data->item_count - 1].size;

	if (tmp.size > 30) {
		printf("Too long itemset, will not continue\n");
		goto out;
	}

	tmp_items = calloc(tmp.size, sizeof(tmp_items[0]));
	tmp.elems = tmp_items;

	masks = calloc(tmp.size, sizeof(masks[0]));
	for (i = 0; i < tmp.size; i++)
		masks[i] = 1 << i;
	
	/* read dataset and then save to histogram */
	for (i = 0; i < data->item_count; i++) {
		if (data->itemsets[i].size < 2)
			continue;

		/* filter only valid exps */
		ve = 0;
		for (j = 0; j < num_exps; j++) {
			l = 0;
			for (i1 = 0; i1 < exps[j].num_items && !l; i1++)
				for (i2 = 0; i2 < data->itemsets[i].size && !l; i2++)
					if (exps[j].items[i1] == data->itemsets[i].elems[i2])
						l++;

			if (l) {
				ve++;
				exps[j].valid = 1;
			}
		}

		if (!ve)
			continue;

		k = 1 << data->itemsets[i].size;
		k--;

		while (--k) {
			l = 0;
			for (j = 0; j < data->itemsets[i].size; j++)
				if (k & masks[j])
					tmp_items[l++] = data->itemsets[i].elems[j];
			tmp.size = l;

			res = bsearch(&tmp, data->itemsets, i + 1, sizeof(data->itemsets[0]), itemset_cmp);
			if (!res) {
				printf("BADD!\n");
				continue;
			}

			c = (data->itemsets[i].count + 0.0) / (res->count + 0.0);

			for (j = 0; j < num_exps; j++) {
				if (!exps[j].valid)
					continue;

				histogram_register(exps[j].h, c);
			}
		}

		for (j = 0; j < num_exps; j++)
			exps[j].valid = 0;
	}

	/* TODO: save histograms to files */
	for (j = 0; j < num_exps; j++) {
		char *fname = NULL, *tmp = NULL;

		asprintf(&fname, "%s/%s_%lu", args.dirname, DATASET, exps[j].num_items);
		for (i = 0; i < exps[j].num_items; i++) {
			asprintf(&tmp, "%s_%d", fname, exps[j].items[i]);
			free(fname);
			fname = strdup(tmp);
			free(tmp);
		}

		printf("Saving exp %lu to %s\n", j, fname);
		f = fopen(fname, "w");
		if (!f) {
			perror("Unable to save output");
			die("Unable to save output");
		}

		fprintf(f, "%lu\n%d", exps[j].num_items, exps[j].items[0]);
		for (i = 1; i < exps[j].num_items; i++)
			fprintf(f, " %d", exps[j].items[i]);
		fprintf(f, "\n\n");
		histogram_dump(f, exps[j].h, 1, "\t");

		fclose(f);
		free(fname);
	}

	status = EXIT_SUCCESS;
out:
	free(args.tfname);
	free(args.dirname);
	if (data) free_itemsets(data);
	if (tmp_items) free(tmp_items);
	if (masks) free(masks);
	for (i = 0; i < num_exps; i++)
		if (exps[i].h)
			free_histogram(exps[i].h);

	exit(status);
}
