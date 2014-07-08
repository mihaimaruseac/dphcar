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
	int n;
	/* number of transactions */
	int t;
	/* header table for the tree, opaque */
	struct table *table;
	/* root of the tree, opaque */
	struct fptree_node *tree;
};

/**
 * Read a transaction file and construct a fp-tree from it.
 */
void fpt_read_from_file(char *fname, struct fptree *fp);

/**
 * Cleanup the data structures used in a fp-tree.
 */
void fpt_cleanup(struct fptree *fp);

int fpt_height(struct fptree *fp);
int fpt_nodes(struct fptree *fp);

/** Debug printing. */
void fpt_tree_print(struct fptree *fp);
void fpt_table_print(struct fptree *fp);

#endif
