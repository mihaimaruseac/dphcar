#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

#include "globals.h"

void init_rng(struct drand48_data *buffer)
{
	struct timeval now;

	gettimeofday(&now, NULL);
	srand48_r(now.tv_usec + now.tv_sec, buffer);
}

int int_cmp(const void *a, const void *b)
{
	const int *ia = a, *ib = b;
	return *ia - *ib;
}

static double laplace(double lambda, struct drand48_data *buffer)
{
	double rnd;

	drand48_r(buffer, &rnd); /* rnd \in [0, 1)      */
	rnd = 0.5 - rnd;         /* rnd \in (-0.5, 0.5] */

	if (signbit(rnd)) /* rnd < 0 */
		return lambda * log(1 + 2 * rnd);
	return -lambda * log(1 - 2 * rnd);
}

double laplace_mechanism(double x, double eps, double sens,
		struct drand48_data *buffer)
{
	return x + laplace(sens/eps, buffer);
}

int bsearch_i(const void *key, const void *base, size_t nmemb, size_t size,
		int (*compar)(const void *, const void *))
{
	int low = 0, high = nmemb - 1, mid, test;

	while (low <= high) {
		mid = low + ((high - low) >> 1);
		test = compar(base + mid * size, key);
		if (test > 0)
			high = mid - 1;
		else
			low = mid + 1;
	}

	return low;
}
