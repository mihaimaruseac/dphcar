/**
 * Global functions and utilities.
 */

#ifndef _GLOBALS_H
#define _GLOBALS_H

#define die(s, ...) \
	do {\
		fprintf(stderr, "[%s: %s %d] "s"\n", __FILE__, \
				__func__, __LINE__, ##__VA_ARGS__); \
		exit(EXIT_FAILURE); \
	} while (0)

struct drand48_data;

/* qsort functions for integer comparisons */
int int_cmp(const void *a, const void *b);
int int_cmp_r(const void *a, const void *b);
int double_cmp(const void *a, const void *b);
int double_cmp_r(const void *a, const void *b);

void init_rng(struct drand48_data *buffer);

/* Laplace mechanism */
double laplace_mechanism(double x, double eps, double sens,
		struct drand48_data *buffer);

/* version of bsearch which returns the rightmost insertion index
 * (the first index for which the element is at least equal to the key)
 */
int bsearch_i(const void *key, const void *base, size_t nmemb, size_t size,
		int (*compar)(const void *, const void *));

#endif
