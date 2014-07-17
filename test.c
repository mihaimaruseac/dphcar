#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LINELEN 1024

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

static int int_cmp(const void *a, const void *b)
{
	const int *ia = a, *ib = b;
	return *ia - *ib;
}


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
	double ts[] = {0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0};
	long long *cs = calloc(sizeof(ts)/sizeof(ts[0]), sizeof(*cs));
	struct file_data *data = NULL;
	int status = EXIT_FAILURE;
	struct itemset tmp, *res;
	unsigned int i, j, k, l;
	int *tmp_items = NULL;
	int *masks = NULL;
	double c;

	if (argc != 2) {
		printf("Need filename\n");
		goto out;
	}

	data = read_input(argv[1]);
	if (!data)
		goto out;

	tmp.size = data->itemsets[data->item_count - 1].size;
	tmp_items = calloc(tmp.size, sizeof(tmp_items[0]));
	tmp.elems = tmp_items;
	
	if (tmp.size > 30) {
		printf("Too long itemset, will not continue\n");
		goto out;
	}

	masks = calloc(tmp.size, sizeof(masks[0]));
	for (i = 0; i < tmp.size; i++)
		masks[i] = 1 << i;

	for (i = 0; i < data->item_count; i++) {
		if (data->itemsets[i].size < 2)
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
			for (j = 0; j < sizeof(ts) / sizeof(ts[0]); j++)
				if (c > ts[j]) {
					cs[j]++;
					break;
				}
		}
	}

	for (i = 0; i < sizeof(ts) / sizeof(ts[0]); i++)
		printf(">%lf: %lld\n", ts[i], cs[i]);

	status = EXIT_SUCCESS;
out:
	if (data) free_itemsets(data);
	if (tmp_items) free(tmp_items);
	if (masks) free(masks);
	exit(status);
}
