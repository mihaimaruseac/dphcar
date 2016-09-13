#include <stdio.h>
#include <stdlib.h>

#include "fp.h"
#include "globals.h"
#include "itstree.h"
#include "recall.h"

struct item_count {
	int value;
	int real_count;
};

static int cmp(const void *a, const void *b)
{
	const struct item_count *ia = a, *ib = b;
	return int_cmp_r(&ia->real_count, &ib->real_count);
}

static void build_items_table(const struct fptree *fp, struct item_count *ic)
{
	size_t i;

	for (i = 0; i < fp->n; i++) {
		ic[i].value = i + 1;
		ic[i].real_count = fpt_item_count(fp, i);
	}

	qsort(ic, fp->n, sizeof(ic[0]), cmp);
}

static size_t ic_search(const struct item_count *ic, size_t n, int x)
{
	size_t i;

	for (i = 0; i < n; i++)
		if (ic[i].value == x)
			return i;

	die("Invalid value in ic_search");
}

static void generate_rules_from_itemset(const int *AB, size_t ab_length,
		const struct fptree *fp, struct itstree_node *itst)
{
	size_t i, j, max, a_length, rc30, rc50, rc70;
	int *cf = calloc(ab_length, sizeof(cf[0]));
	int *A = calloc(ab_length, sizeof(A[0]));
	int sup_ab, sup_a;
	double c;

	max = (1 << ab_length) - 1;
	rc30 = rc50 = rc70 = 0;
	sup_ab = fpt_itemset_count(fp, AB, ab_length);
	for (i = 1; i < max; i++) {
		a_length = 0;
		for (j = 0; j < ab_length; j++)
			if (i & (1 << j))
				A[a_length++] = AB[j];

		sup_a = fpt_itemset_count(fp, A, a_length);
		c = div_or_zero(sup_ab, sup_a);
		if (c >= .3) rc30++;
		if (c >= .5) rc50++;
		if (c >= .7) rc70++;
	}

	for (i = 0; i < ab_length; i++)
		cf[i] = AB[i];
	qsort(cf, ab_length, sizeof(cf[0]), int_cmp);
	record_its(itst, cf, ab_length, rc30, rc50, rc70);

	free(cf);
	free(A);
}

static void generate(const struct fptree *fp, const struct item_count *ic,
		struct itstree_node *itst, size_t ni, size_t lmax,
		int *AB, size_t ab_length)
{
	size_t i, j, st, found, ix = ab_length - 1;

	st = ix > 0 ? ic_search(ic, ni, AB[ix - 1]) : 0;
	for (i = st; i < ni; i++) {
		AB[ix] = ic[i].value;

		found = 0;
		for (j = 0; j < ix && !found; j++)
			if (AB[j] == AB[ix])
				found = 1;
		if (found)
			continue;

		if (ab_length > 1)
			generate_rules_from_itemset(AB, ab_length, fp, itst);
		if (ab_length < lmax)
			generate(fp, ic, itst, ni, lmax, AB, ab_length + 1);
	}
}

static void generate_itemsets(const struct fptree *fp,
		const struct item_count *ic, struct itstree_node *itst,
		size_t ni, size_t lmax)
{
	int *AB = calloc(lmax, sizeof(AB[0]));
	generate(fp, ic, itst, ni, lmax, AB, 1);
	free(AB);
}

struct itstree_node * build_recall_tree(const struct fptree *fp,
		size_t lmax, size_t ni)
{
	struct item_count *ic = calloc(fp->n, sizeof(ic[0]));
	struct itstree_node *itst = init_empty_itstree();

	printf("Building the recall tree ... ");
	build_items_table(fp, ic);
	generate_itemsets(fp, ic, itst, ni, lmax);
	printf("OK\n");

	free(ic);
	return itst;
}
