/**
 * Build structure to be used for recall
 */
#ifndef _RECALL_H
#define _RECALL_H

struct fptree;
struct itstree_node;

struct itstree_node * build_recall_tree(const struct fptree *fp,
		size_t lmax, size_t ni);

#endif
