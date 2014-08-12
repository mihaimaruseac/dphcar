#ifndef _HISTOGRAM_H
#define _HISTOGRAM_H

struct histogram;

struct histogram *init_histogram();

void histogram_register(struct histogram *h, double val);

size_t histogram_get_bin(struct histogram *h, int bin);
double histogram_bin_bound(struct histogram *h, int bin);
size_t histogram_get_all(struct histogram *h);
int histogram_get_count_bins(struct histogram *h);

void histogram_dump(FILE *f, struct histogram *h, int cumulative,
		const char *header);
void histogram_load(FILE *f, struct histogram *h, int cumulative,
		const char *header);

void free_histogram(struct histogram *h);

#endif


