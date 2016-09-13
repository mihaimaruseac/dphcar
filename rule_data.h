/**
 * Data structure to keep data about rules (for recall and to eliminate
 * duplicates).
 */
#ifndef _RULE_DATA_H
#define _RULE_DATA_H

struct rule_data;

struct rule_data *init_rule_data();
void record_new_rule(struct rule_data *rd, const int *cf, size_t sz);
int search_rule_data(const struct rule_data *rd, const int *cf, size_t sz);
void free_rule_data(struct rule_data *rd);

#endif

