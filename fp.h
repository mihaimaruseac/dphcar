/**
 * Fp-tree creation and manipulation.
 */
#ifndef _FP_H
#define _FP_H

struct table;
struct fptree_node;

/**
 * A fp-tree structure.
 *
 * Contains all the information needed to reconstruct a transaction file.
 */
struct fptree {
	/* number of items */
	size_t n;
	/* number of transactions */
	size_t t;
	/* header table for the tree, opaque */
	struct table *table;
	/* root of the tree, opaque */
	struct fptree_node *tree;
};

/**
 * Read a transaction file and construct a fp-tree from it.
 */
void fpt_read_from_file(const char *fname, struct fptree *fp);

/**
 * Cleanup the data structures used in a fp-tree.
 */
void fpt_cleanup(const struct fptree *fp);

int fpt_height(const struct fptree *fp);
int fpt_nodes(const struct fptree *fp);

int fpt_item_count(const struct fptree *fp, int it);
int fpt_itemset_count(const struct fptree *fp, const int *its, int itslen);

/** Debug printing. */
void fpt_tree_print(const struct fptree *fp);
void fpt_table_print(const struct fptree *fp);

#endif
