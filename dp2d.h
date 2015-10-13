/**
 * Differential Privacy using 2D grid.
 */

#ifndef _DP2D_H
#define _DP2D_H

struct fptree;

void dp2d(const struct fptree *fp, double eps, double eps_share,
		size_t k, double c0, double sigma_max, double cmin,
		long int seed);

#endif
