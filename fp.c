#include <ctype.h>
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fp.h"
#include "globals.h"

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

/* sort graph adjacency list in ascending order of nodes */
static int fpgraph_cmp(const void *a, const void *b)
{
	const struct graph *ga=a, *gb=b;
	return ga->n - gb->n;
}

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

#undef LINELENGTH
#undef INITIAL_SIZE

static void read_edges(FILE *f, struct fptree *fp)
{
	int j, *items, sz = 0;

	items = read_line(f, &sz); /* ignore remainder of first line */
	if (sz)
		die("First line in transaction file has bogus data");
	free(items);

	fp->graph = calloc(fp->n, sizeof(fp->graph[0]));

	for (fp->alc = 0; fp->alc < fp->n; fp->alc++) {
		items = read_line(f, &sz);
		if (!sz)
			die("Invalid graph line in transaction file");

		fp->graph[fp->alc].n = items[0];
		fp->graph[fp->alc].cnt = sz - 1;
		fp->graph[fp->alc].next = calloc(sz, sizeof(fp->graph[fp->alc].next[0]));

		for (j = 1; j < sz; j++)
			fp->graph[fp->alc].next[j-1] = items[j];

		free(items);
	}

	qsort(fp->graph, fp->alc, sizeof(fp->graph[0]), fpgraph_cmp);
}

static void read_docs(FILE *f, struct fptree *fp)
{
	int j, *items, sz = 0;
	size_t i;

	items = read_line(f, &sz); /* read the -- separator */
	if (sz)
		die("Missing separator between graph and transactions (or too few transactions)");
	free(items);

	fp->l_max_t = 0;
	fp->docs = calloc(fp->t, sizeof(fp->docs[0]));

	for (i = 0; i < fp->t; i++) {
		items = read_line(f, &sz);
		if (!sz)
			die("Invalid graph line in transaction file");

		if (sz > fp->l_max_t)
			fp->l_max_t = sz;

		/* TODO: check validity of doc path? */
		fp->docs[i].sz = sz;
		fp->docs[i].vals = calloc(sz, sizeof(fp->docs[i].vals[0]));
		for (j = 0; j < sz; j++)
			fp->docs[i].vals[j] = items[j];

		free(items);
	}
}

void fpt_read_from_file(const char *fname, struct fptree *fp)
{
	FILE *f = fopen(fname, "r");

	if (!f)
		die("Invalid transaction filename %s", fname);

	if (fscanf(f, "%lu%lu%lu", &fp->n, &fp->e, &fp->t) != 3)
		die("Invalid header in transaction file %s", fname);

	printf("Reading graph from file ... ");
	fflush(stdout);
	read_edges(f, fp);
	printf("OK\n");

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
		free(fp->graph[i].next);
	for (i = 0; i < fp->n; i++)
		free(fp->docs[i].vals);

	free(fp->graph);
	free(fp->docs);
}

size_t fpt_item_count(const struct fptree *fp, size_t it)
{
	size_t i, j, cnt = 0;

	if (it >= fp->n)
		return 0;

	it++;
	for (i = 0; i < fp->t; i++)
		for (j = 0; j < fp->docs[i].sz; j++)
			if (fp->docs[i].vals[j] == it)
				cnt++;
	return cnt;
}

size_t fpt_itemset_count(const struct fptree *fp,
		const size_t *its, size_t itslen)
{
	size_t i, j, ret=0, found;
	struct docs const *doc;

	for (i = 0; i < fp->t; i++) {
		doc = &fp->docs[i];
		if (doc->sz < itslen)
			continue;

		for (j = 0, found = 0; j < doc->sz; j++) {
			if (doc->vals[j] == its[found]) {
				found++;
				if (found == itslen) {
					found = 0;
					ret++;
				}
			} else
				found = 0;
		}
	}

	return ret;
}

size_t *fp_grph_children(const struct fptree *fp, size_t node, size_t *sz)
{
	struct graph *n;
	size_t *ret, i;

	n = bsearch(&node, fp->graph, fp->alc, sizeof(fp->graph[0]), fpgraph_cmp);
	ret = calloc(n->cnt, sizeof(ret[0]));
	*sz = n->cnt;

	for (i = 0; i < n->cnt; i++)
		ret[i] = n->next[i];

	return ret;
}
