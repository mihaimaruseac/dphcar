#include <ctype.h>
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fp.h"
#include "globals.h"
#include "histogram.h"

static const size_t end_of_transaction = 0;

struct graph {
	/* node */
	size_t n;
	/* count of children */
	size_t cnt;
	/* list of children */
	size_t *next;
};

struct node_data {
	/* label of edge */
	size_t label;
	/* count of transactions passing through edge */
	size_t count;
	/* children */
	struct node *child;
};

struct node {
	/* number of children/labels/counts */
	size_t sz;
	/* allocated space for children/labels/counts */
	size_t asz;
	/* node data */
	struct node_data *data;
};

struct fptree_private {
	/* number of items in adjacency list */
	size_t alc;
	/* adjacency list for the graph */
	struct graph *graph;
	/* documents as trie */
	struct node *root;
	/* unique document occurences udcnt[i] = #_s(i) */
	size_t *udcnt;
};

/* sort graph adjacency list in ascending order of nodes */
static int fpgraph_cmp(const void *a, const void *b)
{
	const struct graph *ga=a, *gb=b;
	return ga->n - gb->n;
}

static int node_data_cmp(const void *a, const void *b)
{
	const struct node_data *nda=a, *ndb=b;
	return nda->label - ndb->label;
}

#define LINELENGTH 4096
#define INITIAL_SIZE 100
static size_t* read_line(FILE *f, size_t *pisz)
{
	size_t i, l, save=0, isz=0, isp=INITIAL_SIZE, *items;
	char line[LINELENGTH], *p, *q;

	items = calloc(isp, sizeof(items[0]));

	while (fgets(line, LINELENGTH, f)) {
		l = strlen(line) - 1;
		p = line;

		/* reform number split between reads */
		if (save) {
			while (isdigit(*p)) {
				save = save * 10 + (*p - '0');
				p++;
			}
			if (isz >= isp) {
				isp *= 2;
				items = realloc(items, isp * sizeof(items[0]));
			}
			items[isz++] = save;
			save = 0;
		}

		/* see if another number is split */
		if (isdigit(line[l])) {
			q = strrchr(line, ' ');
			if (q) {
				save = strtol(q, NULL, 0);
				*q = 0; /* remove last number */
			}
		}

		while ((i = strtol(p, &q, 10))) {
			if (isz >= isp) {
				isp *= 2;
				items = realloc(items, isp * sizeof(items[0]));
			}
			items[isz++] = i;
			p = q;
		}

		if (line[l] == '\n')
			break; /* line completely read */
	}

	*pisz = isz;
	return items;
}
#undef LINELENGTH
#undef INITIAL_SIZE

#define INITIAL_SIZE 10
static struct node *init_node()
{
	struct node *n = calloc(1, sizeof(*n));

	n->asz = INITIAL_SIZE;
	n->data = calloc(n->asz, sizeof(n->data[0]));

	return n;
}
#undef INITIAL_SIZE

static struct node *record_transaction_element(size_t elm,
		struct node *parent)
{
	struct node_data elmd = {elm, 0, 0};
	struct node_data *ndp = bsearch(&elmd, parent->data, parent->sz,
			sizeof(parent->data[0]), node_data_cmp);
	struct node *n;

	if (!ndp) {
		if (parent->sz == parent->asz) {
			parent->asz *= 2;
			parent->data = realloc(parent->data,
					parent->asz * sizeof(parent->data[0]));
		}

		n = init_node();
		parent->data[parent->sz].label = elm;
		parent->data[parent->sz].count = 1;
		parent->data[parent->sz].child = n;
		parent->sz++;

		qsort(parent->data, parent->sz, sizeof(parent->data[0]), node_data_cmp);
		return n;
	}

	ndp->count++;
	return ndp->child;
}

static size_t fpt_search_path(const struct node *n,
		const size_t *pth, size_t pthlen, size_t p)
{
	struct node_data elmd = {pth[p], 0, 0};
	struct node_data *ndp;

	ndp = bsearch(&elmd, n->data, n->sz, sizeof(n->data[0]), node_data_cmp);
	if (!ndp)
		return 0;

	p++;
	if (p == pthlen)
		return ndp->count;

	return fpt_search_path(ndp->child, pth, pthlen, p);
}

static void fpt_add_transaction(size_t *items, size_t st, size_t c, size_t sz,
		const struct fptree *fp, struct node *fpn)
{
	size_t elm = (st + c >= sz) ? end_of_transaction : items[st + c];
	struct node *nn;

	if (c >= fp->l_max_r)
		return;

	nn = record_transaction_element(elm, fpn);
	fpt_add_transaction(items, st, c + 1, sz, fp, nn);
}

static void fpt_mine_path_with(struct histogram *h,
		const struct node *n, size_t x)
{
	size_t i;
	double v;

	for (i = 0; i < n->sz; i++) {
		if (n->data[i].label == end_of_transaction)
			continue;
		v = (n->data[i].count + 0.0) / (x + 0.0);
		histogram_register(h, v);
		fpt_mine_path_with(h, n->data[i].child, x);
	}
}

static void fpt_mine_path(struct histogram *h,
		const struct node *n, size_t x)
{
	size_t i;

	fpt_mine_path_with(h, n, x);
	for (i = 0; i < n->sz; i++)
		if (n->data[i].label != end_of_transaction)
			fpt_mine_path(h, n->data[i].child,
					n->data[i].count);
}

static void free_doc_tree(struct node *n)
{
	size_t i;

	for (i = 0; i < n->sz; i++)
		free_doc_tree(n->data[i].child);

	free(n->data);
	free(n);
}

static void read_edges(FILE *f, struct fptree *fp)
{
	size_t j, *items, sz = 0;

	items = read_line(f, &sz); /* ignore remainder of first line */
	if (sz)
		die("First line in transaction file has bogus data");
	free(items);

	fp->fpt->graph = calloc(fp->n, sizeof(fp->fpt->graph[0]));

	for (fp->fpt->alc = 0; fp->fpt->alc < fp->n; fp->fpt->alc++) {
		items = read_line(f, &sz);
		if (!sz)
			die("Invalid graph line in transaction file");

		fp->fpt->graph[fp->fpt->alc].n = items[0];
		fp->fpt->graph[fp->fpt->alc].cnt = sz - 1;
		fp->fpt->graph[fp->fpt->alc].next = calloc(sz, sizeof(fp->fpt->graph[fp->fpt->alc].next[0]));

		for (j = 1; j < sz; j++)
			fp->fpt->graph[fp->fpt->alc].next[j-1] = items[j];

		free(items);
	}

	qsort(fp->fpt->graph, fp->fpt->alc, sizeof(fp->fpt->graph[0]), fpgraph_cmp);
}

static void read_docs(FILE *f, struct fptree *fp)
{
	size_t i, j, k, sz = 0, *items;

	items = read_line(f, &sz); /* read the -- separator */
	if (sz)
		die("Missing separator between graph and transactions (or too few transactions)");
	free(items);

	fp->l_max_t = 0;
	fp->fpt->root = init_node();

	for (i = 0; i < fp->t; i++) {
		items = read_line(f, &sz);
		if (!sz)
			die("Invalid graph line in transaction file");

		if (sz > fp->l_max_t)
			fp->l_max_t = sz;

		for (j = 0; j < sz; j++)
			fpt_add_transaction(items, j, 0, sz, fp, fp->fpt->root);

		for (j = 0; j < sz; j++) {
			for (k = 0; k < j; k++)
				if (items[j] == items[k])
					break;

			if (k == j)
				fp->fpt->udcnt[items[j]]++;
		}

		free(items);
	}

	/* check for returns */
	fp->has_returns = 1;
	for (i = 1; i <= fp->n && fp->has_returns; i++)
		if (fp->fpt->udcnt[i] != fp->fpt->root->data[i-1].count)
			fp->has_returns = 0;
	fp->has_returns = !fp->has_returns;
}

void fpt_read_from_file(const char *fname, size_t lmax, struct fptree *fp)
{
	FILE *f = fopen(fname, "r");

	if (!f)
		die("Invalid transaction filename %s", fname);

	if (fscanf(f, "%lu%lu%lu", &fp->n, &fp->e, &fp->t) != 3)
		die("Invalid header in transaction file %s", fname);

	fp->fpt = calloc(1, sizeof(*fp->fpt));

	printf("Reading graph from file ... ");
	fflush(stdout);
	read_edges(f, fp);
	printf("OK\n");

	fp->l_max_r = lmax;
	fp->fpt->udcnt = calloc(fp->n + 1, sizeof(fp->fpt->udcnt[0]));

	printf("Reading docs from file ... ");
	fflush(stdout);
	read_docs(f, fp);
	printf("OK\n");

	fclose(f);
}

void fpt_cleanup(const struct fptree *fp)
{
	size_t i;

	for (i = 0; i < fp->n; i++)
		free(fp->fpt->graph[i].next);

	free_doc_tree(fp->fpt->root);
	free(fp->fpt->graph);
	free(fp->fpt->udcnt);
	free(fp->fpt);
}

size_t fpt_item_count(const struct fptree *fp, size_t it)
{
	return fp->fpt->udcnt[it];
}

size_t fpt_itemset_count(const struct fptree *fp,
		const size_t *its, size_t itslen)
{
	return fpt_search_path(fp->fpt->root, its, itslen, 0);
}

size_t *fp_grph_children(const struct fptree *fp, size_t node, size_t *sz)
{
	struct graph *n;
	size_t *ret, i;

	n = bsearch(&node, fp->fpt->graph, fp->fpt->alc, sizeof(fp->fpt->graph[0]), fpgraph_cmp);
	ret = calloc(n->cnt, sizeof(ret[0]));
	*sz = n->cnt;

	for (i = 0; i < n->cnt; i++)
		ret[i] = n->next[i];

	return ret;
}

void fpt_mine(const struct fptree *fp, struct histogram *h)
{
	size_t i;

	for (i = 0; i < fp->fpt->root->sz; i++)
		fpt_mine_path(h, fp->fpt->root->data[i].child,
				fp->fpt->root->data[i].count);
}
