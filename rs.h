/**
 * Generic reservoir sampling.
 */

#ifndef _RS_H
#define _RS_H

/* print actions to the reservoir */
#ifndef PRINT_RS_TRACE
#define PRINT_RS_TRACE 0
#endif
/* detailed trace: print all items, including those ignored, very verbose!! */
#ifndef DETAILED_RS_TRACE
#define DETAILED_RS_TRACE 1
#endif

struct reservoir;
struct reservoir_iterator;
struct drand48_data;

/* Notice one extra parameter when tracing the reservoir */
struct reservoir *init_reservoir(size_t sz,
		void (*print_fun)(const void *it, size_t nmemb, const void *data),
		void *(*clone_fun)(const void *it, size_t nmemb, size_t sz),
		void (*free_fun)(void *it));
void free_reservoir(struct reservoir *r);

/**
 * Add item to reservoir using weight (log weight).
 * To add a single element (struct) set nmemb to 1, sz to sizeof struct.
 * To add an array of elements set nmemb to size of array, sz to sizeof struct.
 */
void add_to_reservoir(struct reservoir *r, const void *it, const void *data,
		size_t nmemb, size_t sz,
		double w, struct drand48_data *randbuffer);
void add_to_reservoir_log(struct reservoir *r, const void *it, const void *data,
		size_t nmemb, size_t sz,
		double logw, struct drand48_data *randbuffer);

struct reservoir_iterator *init_reservoir_iterator(struct reservoir *r);
void free_reservoir_iterator(struct reservoir_iterator *ri);

/**
 * Returns next item in reservoir or NULL if no more items can be found.
 * Updates nmemb and sz if not NULL.
 * Do not free the returned pointer as it is still held on by the reservoir.
 */
const void *next_item(struct reservoir_iterator *ri, size_t *nmemb, size_t *sz);

/* utility functions */
void *shallow_clone(const void *it, size_t nmemb, size_t sz);
void no_print(const void *it, size_t nmemb, const void *data);
void print_int_array(const void *it, size_t nmemb, const void *data);
void print_int(const void *it, size_t unused, const void *data);
void print_size_t_array(const void *it, size_t nmemb, const void *data);
void print_size_t(const void *it, size_t unused, const void *data);

#endif
