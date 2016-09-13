/**
 * Histogram data structure.
 */
#ifndef _HISTOGRAM_H
#define _HISTOGRAM_H

struct histogram;

struct histogram *init_histogram();

void histogram_register(struct histogram *h, double val);

size_t histogram_get_bin(const struct histogram *h, int bin);
double histogram_bin_bound(const struct histogram *h, int bin);
size_t histogram_get_all(const struct histogram *h);
int histogram_get_count_bins(const struct histogram *h);

void histogram_dump(FILE *f, const struct histogram *h, int cumulative,
		const char *header);
void histogram_load(FILE *f, struct histogram *h, int cumulative,
		const char *header);

void free_histogram(struct histogram *h);

#endif


