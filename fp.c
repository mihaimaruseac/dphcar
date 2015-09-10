#include <ctype.h>
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fp.h"
#include "globals.h"

#if 0 /* moving to graphs */
struct fptree_node {
	/* item value */
	int val;
	/* count of item on this path */
	int cnt;
	/* next node in item-chain in tree */
	struct fptree_node *next;
	/* number of children nodes */
	int num_children;
	/* vector of children nodes */
	struct fptree_node **children;
	/* size of said vector */
	int sz_children;
	/* parent in tree */
	struct fptree_node *parent;
};

struct table {
	/* item value */
	size_t val;
	/* count of item */
	size_t cnt;
	/* pointer to first node in item-chain in tree */
	struct fptree_node *fst;
	/* pointer to last node in item-chain in tree */
	struct fptree_node *lst;
	/* reverse permutation index */
	size_t rpi;
};
#endif

struct graph {
	/* node */
	size_t n;
	/* count of children */
	size_t cnt;
	/* list of children */
	size_t *next;
};

struct docs {
	/* length */
	size_t sz;
	/* values */
	size_t *vals;
};

#if 0 /* moving to graphs */
/* sort table entries in descending order */
static int fptable_cmp(const void *a, const void *b)
{
	const struct table *ta = a, *tb = b;
	return tb->cnt - ta->cnt;
}
#endif

#define LINELENGTH 4096
#define INITIAL_SIZE 100

static int* read_line(FILE *f, int *pisz)
{
	int i, l, save=0, isz=0, isp=INITIAL_SIZE, *items;
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

#if 0 /* moving to graphs */
static void single_item_stat(FILE *f, struct fptree *fp)
{
	size_t x, sa = INITIAL_SIZE, *xs, i;
	char line[LINELENGTH];

	fp->t = fp->n = 0;

	while (fgets(line, LINELENGTH, f))
		if (line[strlen(line) - 1] == '\n')
			fp->t++;

	fseek(f, 0, SEEK_SET);
	xs = calloc(sa, sizeof(xs[0]));

	while (fscanf(f, "%lu", &x) == 1) {
		if ((size_t)x > fp->n)
			fp->n = x;
		if (x >= sa) {
			i = sa;
			while (x >= sa)
				sa *= 2;
			xs = realloc(xs, sa * sizeof(xs[0]));
			for (; i < sa; i++)
				xs[i] = 0;
		}
		xs[x]++;
	}

	fp->table = calloc(fp->n, sizeof(fp->table[0]));
	for (i = 0; i < fp->n; i++) {
		fp->table[i].val = i + 1;
		fp->table[i].cnt = xs[i + 1];
		fp->table[i].fst = NULL;
		fp->table[i].lst = NULL;
		fp->table[i].rpi = i;
	}
	qsort(fp->table, fp->n, sizeof(fp->table[0]), fptable_cmp);

	for (i = 0; i < fp->n; i++)
		for (x = 0; x < fp->n; x++)
			if (fp->table[x].val == i + 1)
				fp->table[i].rpi = x;

	free(xs);
	fseek(f, 0, SEEK_SET);
}
#endif

#if 0 /* moving to graphs */
static void fpt_add_transaction(const int *t, int c, int sz,
		struct fptree_node *fpn, struct table *tb);
static void read_transactions(FILE *f, const struct fptree *fp)
{
	int l, i, save = 0, newt, *items, isz = 0, isp = INITIAL_SIZE;
	char line[LINELENGTH], *p, *q;

	items = calloc(isp, sizeof(items[0]));

	while (fgets(line, LINELENGTH, f)) {
		l = strlen(line) - 1;
		newt = line[l] == '\n';
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

		if (newt) {
			for (i = 0; i < isz; i++)
				items[i] = fp->table[items[i]-1].rpi;
			qsort(items, isz, sizeof(items[0]), int_cmp);
			for (i = 0; i < isz; i++)
				items[i] = fp->table[items[i]].val;
			fpt_add_transaction(items, 0, isz, fp->tree, fp->table);
			isz = 0;
		}
	}

	free(items);
}
#endif

#undef LINELENGTH
#undef INITIAL_SIZE
#define INITIAL_SIZE 10

#if 0 /* moving to graphs */
static struct fptree_node *fpt_node_new()
{
	struct fptree_node *ret = calloc(1, sizeof(ret[0]));
	ret->sz_children = INITIAL_SIZE;
	ret->children = calloc(ret->sz_children, sizeof(ret->children[0]));
	return ret;
}
#endif

#undef INITIAL_SIZE

static void read_edges(FILE *f, struct fptree *fp)
{
	int j, *items, sz = 0;
	size_t i;

	items = read_line(f, &sz); /* ignore remainder of first line */
	if (sz)
		die("First line in transaction file has bogus data");
	free(items);

	fp->graph = calloc(fp->n, sizeof(fp->graph[0]));

	for (i = 0; i < fp->n; i++) {
		items = read_line(f, &sz);
		if (!sz)
			die("Invalid graph line in transaction file");

		fp->graph[i].n = items[0];
		fp->graph[i].cnt = sz - 1;
		fp->graph[i].next = calloc(sz, sizeof(fp->graph[i].next[0]));

		for (j = 1; j < sz; j++)
			fp->graph[i].next[j-1] = items[j];

		free(items);
	}
}

static void read_docs(FILE *f, struct fptree *fp)
{
	int j, *items, sz = 0;
	size_t i;

	items = read_line(f, &sz); /* read the -- separator */
	if (sz)
		die("Missing separator between graph and transactions (or too few transactions)");
	free(items);

	fp->docs = calloc(fp->t, sizeof(fp->docs[0]));

	for (i = 0; i < fp->n; i++) {
		items = read_line(f, &sz);
		if (!sz)
			die("Invalid graph line in transaction file");

		/* TODO: check validity of doc path? */
		fp->docs[i].sz = sz;
		fp->docs[i].vals = calloc(sz, sizeof(fp->docs[i].vals[0]));
		for (j = 0; j < sz; j++)
			fp->docs[i].vals[j] = items[j];

		free(items);
	}
}

#if 0 /* moving to graphs */
static void fpt_add_transaction(const int *t, int c, int sz,
		struct fptree_node *fpn, struct table *tb)
{
	struct fptree_node *n;
	int i, elem = t[c];

	if (c >= sz)
		return;

	for (i = 0; i < fpn->num_children; i++)
		if (fpn->children[i]->val == elem) {
			fpn->children[i]->cnt++;
			fpt_add_transaction(t, c + 1, sz, fpn->children[i], tb);
			return;
		}

	if (fpn->num_children == fpn->sz_children) {
		fpn->sz_children *= 2;
		fpn->children = realloc(fpn->children, fpn->sz_children * sizeof(fpn->children));
	}
	n = fpt_node_new();
	n->val = elem;
	n->cnt = 1;
	i = tb[elem-1].rpi;
	if (tb[i].fst == NULL)
		tb[i].fst = tb[i].lst = n;
	else {
		tb[i].lst->next = n;
		tb[i].lst = n;
	}
	fpn->children[fpn->num_children++] = n;
	n->parent = fpn;
	fpt_add_transaction(t, c + 1, sz, n, tb);
}
#endif

#if 0 /* moving to graphs */
static void fpt_node_free(const struct fptree_node *r)
{
	int i;

	for (i = 0; i < r->num_children; i++)
		fpt_node_free(r->children[i]);
	free(r->children);
	free((void*)r);
}
#endif

#if 0 /* moving to graphs */
static int fpt_get_height(const struct fptree_node *r)
{
	int ret = 0, i, x;

	for (i = 0; i < r->num_children; i++) {
		x = fpt_get_height(r->children[i]);
		if (ret < x)
			ret = x;
	}

	return ret + 1;
}
#endif

#if 0 /* moving to graphs */
static int fpt_get_nodes(const struct fptree_node *r)
{
	int ret = 1, i;

	for (i = 0; i < r->num_children; i++)
		ret += fpt_get_nodes(r->children[i]);

	return ret;
}
#endif

#if 0 /* moving to graphs */
static void fpt_node_print(const struct fptree_node *r, int gap)
{
	int i, j;

	if (gap)
		printf("%*c", 2 * gap, ' ');
	printf("%p %d %d %p %d <", r, r->val, r->cnt, r->next, r->num_children);
	for (j = 0; j < r->num_children; j++)
		printf("%p ", r->children[j]);
	printf("> %d\n", r->sz_children);
	for (i = 0; i < r->num_children; i++) {
		fpt_node_print(r->children[i], gap + 1);
	}
}
#endif

void fpt_tree_print(const struct fptree *fp)
{
#if 0 /* moving to graphs */
	fpt_node_print(fp->tree, 0);
#else
	(void)fp;
#endif
}

void fpt_table_print(const struct fptree *fp)
{
#if 0 /* moving to graphs */
	struct table *table = fp->table;
	struct fptree_node *p;
	int n = fp->n;
	int i;

	for (i = 0; i < n; i++) {
		printf("%d] %lu %lu %lu | %p -> %p |", i, table[i].val, table[i].cnt, table[i].rpi, table[i].fst, table[i].lst);
		p = table[i].fst;
		while (p != table[i].lst) {
			printf(" %p", p);
			p = p->next;
		}
		printf(" %p\n", p);
	}
#else
	(void)fp;
#endif
}

void fpt_read_from_file(const char *fname, struct fptree *fp)
{
	FILE *f = fopen(fname, "r");

	if (!f)
		die("Invalid transaction filename %s", fname);

#if 0 /* moving to graphs */
	printf("Reading file to determine counts of items ... ");
	fflush(stdout);
	single_item_stat(f, fp);
	printf("OK\n");

	fp->tree = fpt_node_new();
	printf("Reading file to build fp-tree ... ");
	fflush(stdout);
	read_transactions(f, fp);
#else
	if (fscanf(f, "%lu%lu%lu", &fp->n, &fp->e, &fp->t) != 3)
		die("Invalid header in transaction file %s", fname);

	printf("Reading graph from file ... ");
	fflush(stdout);
	read_edges(f, fp);
	printf("OK\n");

	printf("Reading docs from file ... ");
	fflush(stdout);
	read_docs(f, fp);
#endif
	printf("OK\n");

	fclose(f);
}

void fpt_cleanup(const struct fptree *fp)
{
#if 0 /* moving to graphs */
	free(fp->table);
	fpt_node_free(fp->tree);
#else
	size_t i;

	for (i = 0; i < fp->n; i++)
		free(fp->graph[i].next);
	for (i = 0; i < fp->n; i++)
		free(fp->docs[i].vals);

	free(fp->graph);
	free(fp->docs);
#endif
}

int fpt_height(const struct fptree *fp)
{
#if 0 /* moving to graphs */
	return fpt_get_height(fp->tree);
#else
	(void)fp;
	return 0;
#endif
}

int fpt_nodes(const struct fptree *fp)
{
#if 0 /* moving to graphs */
	return fpt_get_nodes(fp->tree);
#else
	(void)fp;
	return 0;
#endif
}

size_t fpt_item_count(const struct fptree *fp, size_t it)
{
	if (it >= fp->n)
		return 0;

#if 0 /* moving to graphs */
	return fp->table[fp->table[it].rpi].cnt;
#else
	size_t i, j, cnt = 0;

	it++;
	for (i = 0; i < fp->t; i++)
		for (j = 0; j < fp->docs[i].sz; j++)
			if (fp->docs[i].vals[j] == it)
				cnt++;
	return cnt;
#endif
}

#if 0 /* moving to graphs */
static int search_on_path(const struct fptree_node *n,
		const int *key, int keylen)
{
	int i = keylen - 2;
	struct fptree_node *p = n->parent;

	while (p && i >= 0) {
		/* cut */
		if (i > 0 && p->val == key[i-1])
			return 0;
		/* found */
		if (p->val == key[i])
			i--;
		p = p->parent;
	}

	if (!p)
		return 0;

	return n->cnt;
}
#endif

size_t fpt_itemset_count(const struct fptree *fp,
		const size_t *its, size_t itslen)
{
#if 0 /* moving to graphs */
	int *search_key = calloc(itslen, sizeof(search_key[0]));
	int i, count = 0, key_len = 0;
	struct fptree_node *p, *l;

	for (i = 0; i < itslen; i++)
		if (its[i] > 0)
			search_key[key_len++] = fp->table[its[i] - 1].rpi;
	qsort(search_key, key_len, sizeof(search_key[0]), int_cmp);
	for (i = 0; i < key_len; i++)
		search_key[i] = fp->table[search_key[i]].val;

	i = fp->table[search_key[key_len - 1] - 1].rpi;
	p = fp->table[i].fst;
	l = fp->table[i].lst;

	while (p && p != l) {
		count += search_on_path(p, search_key, key_len);
		p = p->next;
	}
	if (p)
		count += search_on_path(p, search_key, key_len);

	free(search_key);
	return count;
#else
	(void)fp; (void)its; (void)itslen;
	return 0;
#endif
}
