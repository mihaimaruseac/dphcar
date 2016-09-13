/**
 * Differential Privacy using 2D grid.
 */
#ifndef _DP2D_H
#define _DP2D_H

struct fptree;
struct itstree_node;

void dp2d(const struct fptree *fp, struct itstree_node *itst,
		double eps, double eps_ratio1, double c0, size_t lmax,
		size_t ni, size_t cspl, long int seed);

#endif
