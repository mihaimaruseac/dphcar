#include <stdio.h>
#include <stdlib.h>

#include "rs.h"

struct reservoir_item {
	void *item_ptr;
	double w;
	double u;
};

struct reservoir {
	struct reservoir_item *its;
	size_t actual;
	size_t sz;
};

struct reservoir *init_reservoir(size_t sz)
{
	struct reservoir *ret = calloc(1, sizeof(*ret));
	ret->sz = sz;
	ret->actual = 0;
	ret->its = calloc(ret->sz, sizeof(ret->its[0]));
	return ret;
}

void free_reservoir(struct reservoir *r)
{
	size_t i;

	for (i = 0;  i < r->actual; i++)
		free(r->its[i].item_ptr);
	free(r->its);
	free(r);
}
