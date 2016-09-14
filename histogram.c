#include <stdio.h>
#include <stdlib.h>

#include "globals.h"
#include "histogram.h"

static const double c_values[] =
	{0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0};
static const int c_num_values = sizeof(c_values) / sizeof(c_values[0]);

struct histogram {
	size_t total;
	size_t *values;
};

struct histogram *init_histogram()
{
	struct histogram *h = calloc(1, sizeof(*h));
	h->values = calloc(c_num_values, sizeof(h->values[0]));
	return h;
}

void histogram_register(struct histogram *h, double val)
{
	int i;

	/* get rid of nan's  and perfect 0's */
	if (val != val || val == 0)
		return;

	h->total++;
	for (i = 0; i < c_num_values; i++)
		if (val > c_values[i]) {
			h->values[i]++;
			return;
		}
}

size_t histogram_get_bin_c(const struct histogram *h, int bin)
{
	size_t ret = 0;

	if (bin < 0 || bin >= c_num_values)
		die("Invalid bin requested");

	for (; bin >= 0; bin--)
		ret += h->values[bin];

	return ret;
}

size_t histogram_get_bin(const struct histogram *h, int bin)
{
	if (bin < 0 || bin >= c_num_values)
		die("Invalid bin requested");
	return h->values[bin];
}

double histogram_bin_bound(const struct histogram *h, int bin)
{
	if (bin < 0 || bin >= c_num_values)
		die("Invalid bin requested");
	h = h;
	return c_values[bin];
}

size_t histogram_get_all(const struct histogram *h)
{
	return h->total;
}

int histogram_get_count_bins(const struct histogram *h)
{
	h = h;
	return c_num_values;
}

void histogram_dump(FILE *f, const struct histogram *h, int cumulative,
		const char *header)
{
	size_t s = 0;
	int i;

	for (i = 0; i < c_num_values; i++) {
		fprintf(f, "%s%3.2f\t%10lu\t%3.2f", header,
				c_values[i], h->values[i],
				(0.0 + h->values[i]) / h->total);
		s += h->values[i];
		if (cumulative)
			fprintf(f, "\t%12lu\t%3.2f", s, (0.0 + s) / h->total);
		fprintf(f, "\n");
	}
}

void histogram_load(FILE *f, struct histogram *h, int cumulative,
		const char *header)
{
	float c, v;
	size_t hv;
	int i;

	for (i = 0; i < c_num_values; i++) {
		if (fscanf(f, header) != 0)
			die("Unable to read header line %d!", i);
		if (fscanf(f, "%f%lu%f", &c, &hv, &v) != 3)
			die("Unable to read sample data line %d!", i);
		if ((int)(100 *c) != (int)(100 * c_values[i]))
			die("Unmatched bin value %d %3.2lf %3.2lf!", i, c, c_values[i]);
		h->values[i] = hv;
		h->total += hv;
		if (cumulative && fscanf(f, "%lu%f", &hv, &c) != 2)
			die("Unable to read cumulative data line %d!", i);
		if (fscanf(f, "\n") != 0)
			die("Extra output on line!");
	}
}

void free_histogram(struct histogram *h)
{
	free(h->values);
	free(h);
}
