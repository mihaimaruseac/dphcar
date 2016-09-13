#include <stdio.h>
#include <stdlib.h>

#include "itstree.h"

struct children_info {
	int item;
	struct itstree_node *iptr;
};

struct itstree_node {
	/* children data */
	struct children_info *children;
	size_t sp;
	size_t sz;
	/* record if rule has been seen in difpriv mining */
	int dpseen;
	/* rule counters (for recall) */
	size_t rc25, rc50;
};

static int cmp(const void *a, const void *b)
{
	const struct children_info *pa = a, *pb = b;
	return pa->item - pb->item;
}

#define INITIALSZ 10
struct itstree_node *init_empty_itstree()
{
	struct itstree_node *ret = calloc(1, sizeof(*ret));
	ret->sp = INITIALSZ;
	ret->children = calloc(ret->sp, sizeof(ret->children[0]));
	return ret;
}
#undef INITIALSZ

#define FILLFACTOR 2
static void do_record_new_rule(struct itstree_node *itst, const int *its,
		size_t sz, int private, size_t rc25, size_t rc50)
{
	struct children_info k, *p;

	if (!sz) {
		if (private)
			itst->dpseen = 1;
		return;
	}

	k.item = its[0];
	p = bsearch(&k, itst->children, itst->sz, sizeof(k), cmp);

	if (!p) {
		if (itst->sz == itst->sp) {
			itst->sp *= FILLFACTOR;
			itst->children = realloc(itst->children,
					itst->sp * sizeof(itst->children[0]));
		}
		itst->children[itst->sz].item = its[0];
		itst->children[itst->sz].iptr = init_empty_itstree();
		itst->sz++;
		qsort(itst->children, itst->sz, sizeof(k), cmp);
		p = bsearch(&k, itst->children, itst->sz, sizeof(k), cmp);
	}

	do_record_new_rule(p->iptr, its+1, sz-1, private, rc25, rc50);
}
#undef FILLFACTOR

void record_its_private(struct itstree_node *itst, const int *its, size_t sz)
{
	do_record_new_rule(itst, its, sz, 1, 0, 0);
}

void record_its(struct itstree_node *itst, const int *its, size_t sz,
		size_t rc25, size_t rc50)
{
	do_record_new_rule(itst, its, sz, 0, rc25, rc50);
}

int search_its_private(const struct itstree_node *itst, const int *its,
		size_t sz)
{
	struct children_info k, *p;

	if (!sz)
		return itst->dpseen;

	k.item = its[0];
	p = bsearch(&k, itst->children, itst->sz, sizeof(k), cmp);

	if (!p)
		return 0;

	return search_its_private(p->iptr, its+1, sz-1);
}

void free_itstree(struct itstree_node *itst)
{
	size_t i;

	for (i = 0; i < itst->sz; i++)
		free_itstree(itst->children[i].iptr);
	free(itst->children);
	free(itst);
}
