#include <stdio.h>
#include <stdlib.h>

#include "rule_data.h"

struct rdnode {
	int item;
	struct rule_data *iptr;
};

struct rule_data {
	struct rdnode *data;
	size_t sp;
	size_t sz;
};

static int rdcmp(const void *a, const void *b)
{
	const struct rdnode *pa = a, *pb = b;
	return pa->item - pb->item;
}

#define INITIALSZ 10
struct rule_data *init_rule_data()
{
	struct rule_data *ret = calloc(1, sizeof(*ret));
	ret->sp = INITIALSZ;
	ret->data = calloc(ret->sp, sizeof(ret->data[0]));
	return ret;
}
#undef INITIALSZ

#define FILLFACTOR 2
void record_new_rule(struct rule_data *rd, const int *cf, size_t sz)
{
	struct rdnode k, *p;

	if (!sz)
		return;

	k.item = cf[0];
	p = bsearch(&k, rd->data, rd->sz, sizeof(k), rdcmp);

	if (!p) {
		if (rd->sz == rd->sp) {
			rd->sp *= FILLFACTOR;
			rd->data = realloc(rd->data,
					rd->sp * sizeof(rd->data[0]));
		}
		rd->data[rd->sz].item = cf[0];
		rd->data[rd->sz].iptr = init_rule_data();
		rd->sz++;
		qsort(rd->data, rd->sz, sizeof(k), rdcmp);
		p = bsearch(&k,rd->data,rd->sz, sizeof(k), rdcmp);
	}

	record_new_rule(p->iptr, cf+1, sz-1);
}
#undef FILLFACTOR

int search_rule_data(const struct rule_data *rd, const int *cf, size_t sz)
{
	struct rdnode k, *p;

	if (!sz)
		return 1;

	k.item = cf[0];
	p = bsearch(&k, rd->data, rd->sz, sizeof(k), rdcmp);

	if (!p)
		return 0;

	return search_rule_data(p->iptr, cf+1, sz-1);
}

void free_rule_data(struct rule_data *rd)
{
	size_t i;

	for (i = 0; i < rd->sz; i++)
		free_rule_data(rd->data[i].iptr);
	free(rd->data);
	free(rd);
}
