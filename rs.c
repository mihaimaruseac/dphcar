#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "rs.h"

struct reservoir_item {
	void *item_ptr;
	double w;
	double v;
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

static inline double generate_random_uniform(struct drand48_data *randbuffer)
{
	double u;
	drand48_r(randbuffer, &u);
	return u;
}

static void store_item(struct reservoir *r, void *it, double w, double v)
{
	/* TODO */
}

void add_to_reservoir(struct reservoir *r, void *it, double w,
		struct drand48_data *randbuffer)
{
	double u = generate_random_uniform(randbuffer);
	double v = -log(u)/w;
	store_item(r, it, w, v);
}

void add_to_reservoir_log(struct reservoir *r, void *it, double logw,
		struct drand48_data *randbuffer)
{
	double u = generate_random_uniform(randbuffer);
	double v = log(log(1/u)) - logw;
	store_item(r, it, logw, v);
}
