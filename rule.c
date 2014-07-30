#include <search.h>
#include <stdlib.h>
#include <stdio.h>

#include "globals.h"
#include "rule.h"

static int itemset_cmp(const void *a, const void *b)
{
	const struct itemset *ia = a, *ib = b;
	size_t i;

	if (ia->length < ib->length)
		return -1;
	if (ia->length > ib->length)
		return 1;

	for (i = 0; i < ia->length; i++) {
		if (ia->items[i] < ib->items[i])
			return -1;
		if (ia->items[i] > ib->items[i])
			return 1;
	}

	return 0;
}

static int rule_cmp(const void *a, const void *b)
{
	struct rule * const *ra = a, * const *rb = b;
	struct rule *rra = *ra, *rrb = *rb;
	int ret;

	ret = itemset_cmp(rra->A, rrb->A);
	if (ret) return ret;
	return itemset_cmp(rra->B, rrb->B);
}

struct rule_table *init_rule_table()
{
	struct rule_table *rt = calloc(1, sizeof(*rt));

	rt->sz = 0;
#define INITIAL 10
	rt->av = INITIAL;
#undef INITIAL
	rt->rules = calloc(rt->av, sizeof(rt->rules[0]));
	rt->supA = calloc(rt->av, sizeof(rt->supA[0]));
	rt->supAB = calloc(rt->av, sizeof(rt->supAB[0]));
	rt->c = calloc(rt->av, sizeof(rt->c[0]));

	return rt;
}

void save_rule(struct rule_table *rt, const struct rule *r,
		int supA, int supAB)
{
	struct rule **p;
	size_t ix;

	if (rt->sz == rt->av) {
#define INCREASE_RATE 2
		rt->av *= INCREASE_RATE;
#undef INCREASE_RATE
		rt->rules = realloc(rt->rules, rt->av * sizeof(rt->rules[0]));
		rt->supA = realloc(rt->supA, rt->av * sizeof(rt->supA[0]));
		rt->supAB = realloc(rt->supAB, rt->av * sizeof(rt->supAB[0]));
		rt->c = realloc(rt->c, rt->av * sizeof(rt->c[0]));
	}

	p = lsearch(&r, rt->rules, &(rt->sz), sizeof(r), rule_cmp);
	ix = p - &rt->rules[0];

	if (ix != rt->sz - 1) {
		if (rt->supA[ix] != supA || rt->supAB[ix] != supAB)
			die("Same rule with different sup %d %d / %d %d",
				rt->supA[ix], rt->supAB[ix], supA, supAB);
	} else {
		rt->supA[ix] = supA;
		rt->supAB[ix] = supAB;
		rt->c[ix] = supAB / (supA + 0.0);
	}
}

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
	size_t l = A->length, i, j = 0;
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

	if (l != A->length)
		die("Length changed from %lu to %lu", AB->length, l);

	return r;
}

struct itemset *build_itemset(const int *items, size_t length)
{
	struct itemset *its = calloc(1, sizeof(*its));
	size_t i;

	its->length = length;
	its->items = calloc(its->length, sizeof(its->items[0]));
	its->length = 0;

	for (i = 0; i < length; i++)
		if (items[i])
			its->items[its->length++] = items[i];

	its->items = realloc(its->items, its->length * sizeof(its->items[0]));

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

void free_rule_table(struct rule_table *rt)
{
	size_t i;

	if (!rt)
		return;

	for (i = 0; i < rt->sz; i++) {
		// TODO: clean
		print_rule(rt->rules[i]);
		printf(" (%d, %d) %lf\n", rt->supA[i], rt->supAB[i], rt->c[i]);
		free_rule(rt->rules[i]);
	}

	free(rt->supA);
	free(rt->supAB);
	free(rt->c);
	free(rt->rules);
	free(rt);
}
