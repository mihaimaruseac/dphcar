/**
 * Generic reservoir sampling.
 */

#ifndef _RS_H
#define _RS_H

struct reservoir;
struct drand48_data;

struct reservoir *init_reservoir(size_t sz);
void free_reservoir(struct reservoir *r);

/* Notice one extra parameter when tracing the reservoir */
void add_to_reservoir(struct reservoir *r, void *it, double w,
		struct drand48_data *randbuffer);
void add_to_reservoir_log(struct reservoir *r, void *it, double logw,
		struct drand48_data *randbuffer);

#endif
