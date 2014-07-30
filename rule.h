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
	int sz;
	int av;
	struct rule *rules;
};

void print_rule(const struct rule *r);

/* build rule from A and B */
struct rule *build_rule_A_B(const struct itemset *A, const struct itemset *B);
/* build rule from A and AB */
struct rule *build_rule_A_AB(const struct itemset *A, const struct itemset *AB);

struct itemset *build_itemset(const int *items, int length);

void free_rule(struct rule *r);
void free_itemset(struct itemset *its);

#endif

