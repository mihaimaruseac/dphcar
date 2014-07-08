/**
 * Differential Privacy using 2D grid.
 */

#ifndef _DP2D_H
#define _DP2D_H

struct fptree;

void dp2d(struct fptree *fp, double c, double eps, double eps_share, int ni);

#endif
