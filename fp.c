#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fp.h"
#include "globals.h"

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
};

struct table {
	/* item value */
	int val;
	/* count of item */
	int cnt;
	/* pointer to first node in item-chain in tree */
	struct fptree_node *fst;
	/* pointer to last node in item-chain in tree */
	struct fptree_node *lst;
	/* reverse permutation index */
	int rpi;
};

/* sort table entries in descending order */
static int fptable_cmp(const void *a, const void *b)
{
	const struct table *ta = a, *tb = b;
	return tb->cnt - ta->cnt;
}

#define LINELENGTH 4 //096
#define INITIAL_SIZE 100

static void single_item_stat(FILE *f, struct fptree *fp)
{
	int x, sa = INITIAL_SIZE, *xs, i;
	char line[LINELENGTH];

	fp->t = fp->n = 0;

	while (fgets(line, LINELENGTH, f))
		if (line[strlen(line) - 1] == '\n')
			fp->t++;

	fseek(f, 0, SEEK_SET);
	xs = calloc(sa, sizeof(xs[0]));

	while (fscanf(f, "%d", &x) == 1) {
		if (x > fp->n)
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

static void fpt_add_transaction(int *t, int c, int sz, struct fptree *fp);
static void read_transactions(FILE *f, struct fptree *fp)
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
			save = strtol(q, NULL, 0);
			*q = 0; /* remove last number */
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
			fpt_add_transaction(items, 0, isz, fp);
			isz = 0;
		}
	}

	free(items);
}

#undef LINELENGTH
#undef INITIAL_SIZE
#define INITIAL_SIZE 10

static struct fptree_node *fpt_node_new()
{
	struct fptree_node *ret = calloc(1, sizeof(ret[0]));
	ret->sz_children = INITIAL_SIZE;
	ret->children = calloc(ret->sz_children, sizeof(ret->children[0]));
	return ret;
}

#undef INITIAL_SIZE

static void fpt_add_transaction(int *t, int c, int sz, struct fptree *fp)
{
	int i, elem = t[c];

#if 0
	for (i = 0; i < fp->num_children; i++)
		if (fp->children[i].val == t)
#else
	for (i = 0; i < sz; i++)
#endif
		printf("%d ", t[i]);
	printf("\n");
}

static void fpt_node_free(struct fptree_node *r)
{
	int i;

	for (i = 0; i < r->num_children; i++)
		fpt_node_free(r->children[i]);
	free(r->children);
	free(r);
}

void fpt_read_from_file(char *fname, struct fptree *fp)
{
	FILE *f = fopen(fname, "r");

	if (!f)
		die("Invalid transaction filename %s", fname);

	printf("Reading file to determine counts of items ... ");
	fflush(stdout);
	single_item_stat(f, fp);
	printf("OK\n");

	fp->tree = fpt_node_new();
	printf("Reading file to build fp-tree ... ");
	fflush(stdout);
	read_transactions(f, fp);
	printf("OK\n");

	fclose(f);
}

void fpt_cleanup(struct fptree *fp)
{
	free(fp->table);
	fpt_node_free(fp->tree);
}
