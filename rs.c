#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "globals.h"
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
	/* utility functions */
#if PRINT_RS_TRACE
	void (*print_fun)(void *it);
#endif
	void *(*clone_fun)(const void *it);
	void (*free_fun)(void *it);
};

struct reservoir *init_reservoir(size_t sz,
#if PRINT_RS_TRACE
		void (*print_fun)(void *it),
#endif
		void *(*clone_fun)(const void *it),
		void (*free_fun)(void *it))
{
	struct reservoir *ret = calloc(1, sizeof(*ret));
	ret->its = calloc(sz, sizeof(ret->its[0]));
	ret->actual = 0;
	ret->sz = sz;
#if PRINT_RS_TRACE
	ret->print_fun = print_fun;
#endif
	ret->clone_fun = clone_fun;
	ret->free_fun = free_fun;
	return ret;
}

void free_reservoir(struct reservoir *r)
{
	size_t i;

	for (i = 0;  i < r->actual; i++)
		r->free_fun(r->its[i].item_ptr);
	free(r->its);
	free(r);
}

static inline double generate_random_uniform(struct drand48_data *randbuffer)
{
	double u;
	drand48_r(randbuffer, &u);
	return u;
}

static void store_item_at(struct reservoir *r, size_t ix,
		const void *it, double w, double v)
{
	r->its[ix].item_ptr = r->clone_fun(it);
	r->its[ix].w = w;
	r->its[ix].v = v;
}

static int reservoir_cmp(const void *a, const void *b)
{
	const struct reservoir_item *ra = a, *rb = b;
	return double_cmp(&ra->v, &rb->v);
}

#if PRINT_RS_TRACE
static void print_reservoir(struct reservoir *r)
{
	size_t i;

	printf("Reservoir now:\n");
	for (i = 0; i < r->actual; i++) {
		printf("\t|");
		r->print_fun(&r->its[i].item_ptr);
		printf("| w=%lf v=%lf\n", r->its[i].w, r->its[i].v);
	}
}
#endif

static void store_item(struct reservoir *r, const void *it,
		double w, double v)
{
	/* not a full reservoir yet */
	if (r->actual < r->sz) {
		store_item_at(r, r->actual, it, w, v);
		r->actual++;
		goto end;
	}

	/* no changes to the reservoir */
	if (v >= r->its[r->sz - 1].v)
		return;

	r->free_fun(r->its[r->sz-1].item_ptr);
	store_item_at(r, r->sz - 1, it, w, v);

end:
	if (r->actual == r->sz) {
		qsort(r->its, r->sz, sizeof(r->its[0]), reservoir_cmp);
#if PRINT_RS_TRACE
		print_reservoir(r);
#endif
	}
}

void add_to_reservoir(struct reservoir *r, const void *it, double w,
		struct drand48_data *randbuffer)
{
	double u = generate_random_uniform(randbuffer);
	double v = -log(u)/w;
	store_item(r, it, w, v);
}

void add_to_reservoir_log(struct reservoir *r, const void *it, double logw,
		struct drand48_data *randbuffer)
{
	double u = generate_random_uniform(randbuffer);
	double v = log(log(1/u)) - logw;
	store_item(r, it, logw, v);
}
