/**
 * Differential Privacy using 2D grid.
 */

#ifndef _DP2D_H
#define _DP2D_H

struct fptree;

void dp2d(const struct fptree *fp, const char *npfile,
		double eps, double eps_share, int minth, size_t mis, size_t k,
		double minalpha, long int seed);

#endif
