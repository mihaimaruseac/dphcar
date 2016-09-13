/**
 * Data structure to keep data about rules (for recall and to eliminate
 * duplicates).
 */
#ifndef _RULE_DATA_H
#define _RULE_DATA_H

struct seen_data;

struct seen *init_seen_node();
void record_new_seen(struct seen *seen, const int *cf, size_t sz);
int search_seen(const struct seen *seen, const int *cf, size_t sz);
void free_seen_node(struct seen *seen);

#endif

