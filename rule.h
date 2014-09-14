#ifndef _RULE_H
#define _RULE_H

struct itemset {
	size_t length;
	int *items;
};

struct rule {
	struct itemset *A, *B;
};

struct rule_table {
	size_t sz;
	size_t av;
	struct rule **rules;
	int *supA;
	int *supAB;
	double *c;
};

struct rule_table *init_rule_table();

void save_rule(struct rule_table *rt, const struct rule *r,
		int supA, int supAB);
void save_rule2(struct rule_table *rt, const struct rule *r, double c);

void print_rule(const struct rule *r);

/* build rule from A and B */
struct rule *build_rule_A_B(const struct itemset *A, const struct itemset *B);
/* build rule from A and AB */
struct rule *build_rule_A_AB(const struct itemset *A, const struct itemset *AB);

struct itemset *build_itemset(const int *items, size_t length);
struct itemset *build_itemset_add_items(const struct itemset *base, int *item, size_t length);
struct itemset *build_itemset_del_items(const struct itemset *base, int *item, size_t length);

void free_rule(struct rule *r);
void free_itemset(struct itemset *its);
void free_rule_table(struct rule_table *rt);
void free_rule_table2(struct rule_table *rt);

#endif

