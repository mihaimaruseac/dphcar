#include <stdio.h>
#include <stdlib.h>

#include "rule_data.h"

struct seen_data {
	int item;
	struct seen *iptr;
};

struct seen {
	struct seen_data *data;
	size_t sp;
	size_t sz;
};

static int seen_data_cmp(const void *a, const void *b)
{
	const struct seen_data *pa = a, *pb = b;
	return pa->item - pb->item;
}

#define INITIALSZ 10
struct seen *init_seen_node()
{
	struct seen *ret = calloc(1, sizeof(*ret));
	ret->sp = INITIALSZ;
	ret->data = calloc(ret->sp, sizeof(ret->data[0]));
	return ret;
}
#undef INITIALSZ

#define FILLFACTOR 2
void record_new_seen(struct seen *seen, const int *cf, size_t sz)
{
	struct seen_data k, *p;

	if (!sz)
		return;

	k.item = cf[0];
	p = bsearch(&k, seen->data, seen->sz, sizeof(k), seen_data_cmp);

	if (!p) {
		if (seen->sz == seen->sp) {
			seen->sp *= FILLFACTOR;
			seen->data = realloc(seen->data,
					seen->sp * sizeof(seen->data[0]));
		}
		seen->data[seen->sz].item = cf[0];
		seen->data[seen->sz].iptr = init_seen_node();
		seen->sz++;
		qsort(seen->data, seen->sz, sizeof(k), seen_data_cmp);
		p = bsearch(&k,seen->data,seen->sz, sizeof(k), seen_data_cmp);
	}

	record_new_seen(p->iptr, cf+1, sz-1);
}
#undef FILLFACTOR

int search_seen(const struct seen *seen, const int *cf, size_t sz)
{
	struct seen_data k, *p;

	if (!sz)
		return 1;

	k.item = cf[0];
	p = bsearch(&k, seen->data, seen->sz, sizeof(k), seen_data_cmp);

	if (!p)
		return 0;

	return search_seen(p->iptr, cf+1, sz-1);
}

void free_seen_node(struct seen *seen)
{
	size_t i;

	for (i = 0; i < seen->sz; i++)
		free_seen_node(seen->data[i].iptr);
	free(seen->data);
	free(seen);
}

