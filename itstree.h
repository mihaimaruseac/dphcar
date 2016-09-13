/**
 * Tree of itemsets (for recall and duplicate removal).
 */
#ifndef _ITSTREE_H
#define _ITSTREE_H

struct itstree_node;

struct itstree_node *init_empty_itstree();
void free_itstree(struct itstree_node *itst);

void record_its_private(struct itstree_node *itst, const int *its, size_t sz);
void record_its(struct itstree_node *itst, const int *its, size_t sz,
		size_t rc30, size_t rc50, size_t rc70);

int search_its_private(const struct itstree_node *itst, const int *its,
		size_t sz);

void save_its(const struct itstree_node *itst, const char *fname,
		size_t lmax, size_t ni);

#endif

