/**
 * Fp-tree creation and manipulation.
 */
#ifndef _FP_H
#define _FP_H

#if 0 /* moving to graphs */
struct table;
struct fptree_node;
#endif
struct graph;

/**
 * A fp-tree structure.
 *
 * Contains all the information needed to reconstruct a transaction file.
 */
struct fptree {
	/* number of nodes */
	size_t n;
	/* number of edges */
	size_t e;
	/* number of transactions */
	size_t t;
#if 0 /* moving to graphs */
	/* header table for the tree, opaque */
	struct table *table;
	/* root of the tree, opaque */
	struct fptree_node *tree;
#else
	/* number of items in adjacency list */
	size_t alc;
	/* adjacency list for the graph, opaque */
	struct graph *graph;
	/* documents, opaque */
	struct docs *docs;
#endif
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

size_t fpt_item_count(const struct fptree *fp, size_t it);
size_t fpt_itemset_count(const struct fptree *fp,
		const size_t *its, size_t itslen);

/** Debug printing. */
void fpt_tree_print(const struct fptree *fp);
void fpt_table_print(const struct fptree *fp);

size_t *fp_grph_children(const struct fptree *fp, size_t node, size_t *sz);

#endif
