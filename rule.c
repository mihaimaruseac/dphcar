#include <search.h>
#include <stdlib.h>
#include <stdio.h>

#include "globals.h"
#include "rule.h"

void print_rule(const struct rule *r)
{
	size_t i;

	for (i = 0; i < r->A->length; i++)
		printf("%d ", r->A->items[i]);
	printf("-> ");
	for (i = 0; i < r->B->length; i++)
		if (r->B->items[i])
			printf("%d ", r->B->items[i]);
}

struct rule *build_rule_A_B(const struct itemset *A, const struct itemset *B)
{
	struct rule *r = calloc(1, sizeof(*r));
	size_t i;

	r->A = calloc(1, sizeof(*(r->A)));
	r->B = calloc(1, sizeof(*(r->B)));

	r->A->length = A->length;
	r->B->length = B->length;

	r->A->items = calloc(A->length, sizeof(r->A->items[0]));
	r->B->items = calloc(B->length, sizeof(r->B->items[0]));

	for (i = 0; i < A->length; i++)
		r->A->items[i] = A->items[i];
	for (i = 0; i < B->length; i++)
		r->B->items[i] = B->items[i];

	return r;
}

struct rule *build_rule_A_AB(const struct itemset *A, const struct itemset *AB)
{
	struct rule *r = calloc(1, sizeof(*r));
	size_t l = AB->length, i, j = 0;
	int *p;

	r->A = calloc(1, sizeof(*(r->A)));
	r->B = calloc(1, sizeof(*(r->B)));

	r->A->length = A->length;
	r->B->length = AB->length - A->length;

	r->A->items = calloc(r->A->length, sizeof(r->A->items[0]));
	r->B->items = calloc(r->B->length, sizeof(r->B->items[0]));

	for (i = 0; i < r->A->length; i++)
		r->A->items[i] = A->items[i];

	for (i = 0; i < AB->length; i++) {
		p = lfind(&AB->items[i], r->A->items, &l,
				sizeof(r->A->items[0]), int_cmp);
		if (!p)
			r->B->items[j++] = AB->items[i];
	}

	if (j != r->B->length)
		die("Length of B: expected %lu got %lu", r->B->length, j);

	if (l != AB->length)
		die("Length changed from %lu to %lu", AB->length, l);

	return r;
}

struct itemset *build_itemset(const int *items, int length)
{
	struct itemset *its = calloc(1, sizeof(*its));
	size_t i;

	its->length = length;
	its->items = calloc(its->length, sizeof(its->items[0]));

	for (i = 0; i < its->length; i++)
		its->items[i] = items[i];

	return its;
}

void free_rule(struct rule *r)
{
	if (!r)
		return;

	free_itemset(r->A);
	free_itemset(r->B);
	free(r);
}

void free_itemset(struct itemset *its)
{
	if (!its)
		return;

	free(its->items);
	free(its);
}
