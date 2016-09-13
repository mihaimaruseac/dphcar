/**
 * Tree of itemsets (for recall and duplicate removal).
 */
#ifndef _ITSTREE_H
#define _ITSTREE_H

struct itstree_node;

struct itstree_node *init_empty_itstree();
void record_new_rule(struct itstree_node *itst, const int *cf, size_t sz);
int search_rule(const struct itstree_node *itst, const int *cf, size_t sz);
void free_itstree(struct itstree_node *itst);

#endif

