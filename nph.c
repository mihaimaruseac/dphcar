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

#define NUM 50
#define MAX_SIZE 5
#define TOP 100

struct experiment {
	size_t num_items;
	int items[MAX_SIZE];
	struct histogram *h;
	int valid;
};

/* Command line arguments */
static struct {
	/* filename containing the transactions */
	char *tfname;
	/* dirname for the outputs */
	char *dirname;
	/* dataset */
	char *dataset;
} args;

static void usage(const char *prg)
{
	fprintf(stderr, "Usage: %s iFILE DATASET OUTDIRNAME\n", prg);
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
	args.dataset = strdup(argv[2]);
	args.dirname = strdup(argv[3]);
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

static int itemset_cnt_cmp(const void *isa, const void *isb)
{
	const struct itemset *a = isa;
	const struct itemset *b = isb;

	return b->count - a->count;
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
	struct experiment *exps = calloc(NUM * MAX_SIZE, sizeof(exps[0]));
	size_t i, j, k, l, ve, i1, i2, num_exps = 0;
	struct drand48_data randbuffer;
	struct file_data *data = NULL;
	int status = EXIT_FAILURE;
	struct itemset tmp, *res;
	int *tmp_items = NULL;
	int *masks = NULL;
	double c;
	FILE *f;

	parse_arguments(argc, argv);
	srand48_r(42, &randbuffer);

	data = read_input(argv[1]);
	if (!data)
		goto out;

	qsort(data->itemsets, data->item_count, sizeof(data->itemsets[0]), itemset_cnt_cmp);
	for (i = 2; i <= MAX_SIZE; i++) {
		k = 0;
		for (j = 0; j < NUM; j++) {
			exps[num_exps].num_items = i;
			for (; k < data->item_count; k++) {
				if (data->itemsets[k].size != i)
					continue;
				for (l = 0; l < i; l++)
					exps[num_exps].items[l] = data->itemsets[k].elems[l];
				num_exps++;
				break;
			}

			if (k++ == data->item_count)
				break;
		}
	}
	qsort(data->itemsets, data->item_count, sizeof(data->itemsets[0]), itemset_cmp);

	printf("Exps: %lu\n", num_exps);

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
			for (i1 = 0; i1 < exps[j].num_items; i1++)
				for (i2 = 0; i2 < data->itemsets[i].size; i2++)
					if (exps[j].items[i1] == data->itemsets[i].elems[i2]) {
						l++;
						break;
					}

			if (l == exps[j].num_items) {
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

	/* save histograms to files */
	for (j = 0; j < num_exps; j++) {
		char *fname = NULL, *tmp = NULL;

		asprintf(&fname, "%s/%s_%lu", args.dirname, args.dataset, exps[j].num_items);
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
