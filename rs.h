/**
 * Generic reservoir sampling.
 */

#ifndef _RS_H
#define _RS_H

/* print actions to the reservoir */
#ifndef PRINT_RS_TRACE
#define PRINT_RS_TRACE 1
#endif

struct reservoir;
struct drand48_data;

/* Notice one extra parameter when tracing the reservoir */
struct reservoir *init_reservoir(size_t sz,
#if PRINT_RS_TRACE
		void (*print_fun)(void *it),
#endif
		void *(*clone_fun)(const void *it, size_t sz),
		void (*free_fun)(void *it));
void free_reservoir(struct reservoir *r);

void add_to_reservoir(struct reservoir *r, const void *it, size_t sz,
		double w, struct drand48_data *randbuffer);
void add_to_reservoir_log(struct reservoir *r, const void *it, size_t sz,
		double logw, struct drand48_data *randbuffer);

/* utility functions */
void *shallow_clone(const void *it, size_t sz);
void no_print(void *it);

#endif
