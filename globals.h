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

/* qsort function for integer comparisons */
int int_cmp(const void *a, const void *b);

#endif
