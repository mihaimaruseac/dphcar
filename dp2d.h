/**
 * Differential Privacy using 2D grid.
 */

#ifndef _DP2D_H
#define _DP2D_H

struct fptree;

void dp2d(const struct fptree *fp, double eps, double eps_ratio1,
		double c0, size_t lmax, size_t k, long int seed);

/* non private version */
void dp2d_np(const struct fptree *fp, double c0, size_t lmax, long int seed);

#endif
