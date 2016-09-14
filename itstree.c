#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

#include "globals.h"
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
	size_t rc30, rc50, rc70;
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
		size_t sz, int private, size_t rc30, size_t rc50, size_t rc70)
{
	struct children_info k, *p;

	if (!sz) {
		if (private)
			itst->dpseen = 1;
		else {
			itst->rc30 = rc30;
			itst->rc50 = rc50;
			itst->rc70 = rc70;
		}
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

	do_record_new_rule(p->iptr, its+1, sz-1, private, rc30, rc50, rc70);
}
#undef FILLFACTOR

void record_its_private(struct itstree_node *itst, const int *its, size_t sz)
{
	do_record_new_rule(itst, its, sz, 1, 0, 0, 0);
}

void record_its(struct itstree_node *itst, const int *its, size_t sz,
		size_t rc30, size_t rc50, size_t rc70)
{
	do_record_new_rule(itst, its, sz, 0, rc30, rc50, rc70);
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

static void save_its_node(FILE *f, const struct itstree_node *n)
{
	size_t i;
	int item;

	fwrite(&n->sp,     sizeof(n->sp),     1, f);
	fwrite(&n->sz,     sizeof(n->sz),     1, f);
	fwrite(&n->dpseen, sizeof(n->dpseen), 1, f);
	fwrite(&n->rc30,   sizeof(n->rc30),   1, f);
	fwrite(&n->rc50,   sizeof(n->rc50),   1, f);
	fwrite(&n->rc70,   sizeof(n->rc70),   1, f);

	for (i = 0; i < n->sz; i++) {
		item = n->children[i].item;
		fwrite(&item, sizeof(item), 1, f);
		save_its_node(f, n->children[i].iptr);
	}
}

static struct itstree_node *read_its_node(FILE *f)
{
	struct itstree_node *ret = init_empty_itstree();
	size_t i;
	int item;

	fread(&ret->sp,     sizeof(ret->sp),     1, f);
	fread(&ret->sz,     sizeof(ret->sz),     1, f);
	fread(&ret->dpseen, sizeof(ret->dpseen), 1, f);
	fread(&ret->rc30,   sizeof(ret->rc30),   1, f);
	fread(&ret->rc50,   sizeof(ret->rc50),   1, f);
	fread(&ret->rc70,   sizeof(ret->rc70),   1, f);
	ret->children = realloc(ret->children, ret->sp*sizeof(ret->children[0]));

	for (i = 0; i < ret->sz; i++) {
		fread(&item, sizeof(item), 1, f);
		ret->children[i].item = item;
		ret->children[i].iptr = read_its_node(f);
	}

	return ret;
}

void save_its(const struct itstree_node *itst, const char *fname,
		size_t lmax, size_t ni)
{
	char *filename = NULL;
	FILE *f;

	asprintf(&filename, "%s_%lu_%lu", fname, lmax, ni);
	f = fopen(filename, "w");

	if (!f)
		die("Unable to save file %s", filename);

	printf("Saving its to %s ... ", filename);
	fwrite(&lmax, sizeof(lmax), 1, f);
	fwrite(&ni,   sizeof(ni),   1, f);
	save_its_node(f, itst);
	printf("OK\n");

	fclose(f);
	free(filename);
}

struct itstree_node *load_its(const char *fname, size_t lmax, size_t ni)
{
	FILE *f = fopen(fname, "r");
	struct itstree_node *ret;
	size_t lmaxc, nic;

	if (!f)
		die("Unable to read itemset tree from %s", fname);

	printf("Loading its ... ");
	fread(&lmaxc, sizeof(lmaxc), 1, f);
	fread(&nic,   sizeof(nic),   1, f);
	if (lmaxc != lmax || nic != ni)
		die("Itemset tree input filename %s for wrong settings", fname);

	ret = read_its_node(f);
	printf("OK\n");

	fclose(f);
	return ret;
}

static void do_count(const struct itstree_node *itst,
		size_t *n30, size_t *n50, size_t *n70, int private)
{
	size_t i;

	if (!private || itst->dpseen) {
		*n30 += itst->rc30;
		*n50 += itst->rc50;
		*n70 += itst->rc70;
	}

	for (i = 0; i < itst->sz; i++)
		do_count(itst->children[i].iptr, n30, n50, n70, private);
}

void itstree_count_real(const struct itstree_node *itst,
		size_t *n30, size_t *n50, size_t *n70)
{
	do_count(itst, n30, n50, n70, 0);
}

void itstree_count_priv(const struct itstree_node *itst,
		size_t *p30, size_t *p50, size_t *p70)
{
	do_count(itst, p30, p50, p70, 1);
}
