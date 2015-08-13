/**
 * Differential Privacy using 2D grid.
 */

#ifndef _DP2D_H
#define _DP2D_H

struct fptree;

enum bin_mode {
	NONE,
	RANDOM,
	EQUIWIDTH,
	EQUIDENSITY,
};

void dp2d(const struct fptree *fp,
		size_t shelves, size_t bins, enum bin_mode bin_mode,
		double eps, double eps_share, int minth, size_t mis, size_t k,
		double minalpha, long int seed);

#endif
