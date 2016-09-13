#include <stdio.h>
#include <stdlib.h>

#include "itstree.h"

struct node_info {
	int item;
	struct itstree_node *iptr;
};

struct itstree_node {
	struct node_info *data;
	size_t sp;
	size_t sz;
};

static int cmp(const void *a, const void *b)
{
	const struct node_info *pa = a, *pb = b;
	return pa->item - pb->item;
}

#define INITIALSZ 10
struct itstree_node *init_empty_itstree()
{
	struct itstree_node *ret = calloc(1, sizeof(*ret));
	ret->sp = INITIALSZ;
	ret->data = calloc(ret->sp, sizeof(ret->data[0]));
	return ret;
}
#undef INITIALSZ

#define FILLFACTOR 2
void record_new_rule(struct itstree_node *itst, const int *cf, size_t sz)
{
	struct node_info k, *p;

	if (!sz)
		return;

	k.item = cf[0];
	p = bsearch(&k, itst->data, itst->sz, sizeof(k), cmp);

	if (!p) {
		if (itst->sz == itst->sp) {
			itst->sp *= FILLFACTOR;
			itst->data = realloc(itst->data,
					itst->sp * sizeof(itst->data[0]));
		}
		itst->data[itst->sz].item = cf[0];
		itst->data[itst->sz].iptr = init_empty_itstree();
		itst->sz++;
		qsort(itst->data, itst->sz, sizeof(k), cmp);
		p = bsearch(&k, itst->data, itst->sz, sizeof(k), cmp);
	}

	record_new_rule(p->iptr, cf+1, sz-1);
}
#undef FILLFACTOR

int search_rule(const struct itstree_node *itst, const int *cf, size_t sz)
{
	struct node_info k, *p;

	if (!sz)
		return 1;

	k.item = cf[0];
	p = bsearch(&k, itst->data, itst->sz, sizeof(k), cmp);

	if (!p)
		return 0;

	return search_rule(p->iptr, cf+1, sz-1);
}

void free_itstree(struct itstree_node *itst)
{
	size_t i;

	for (i = 0; i < itst->sz; i++)
		free_itstree(itst->data[i].iptr);
	free(itst->data);
	free(itst);
}
