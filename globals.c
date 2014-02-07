#include "globals.h"

int int_cmp(const void *a, const void *b)
{
	const int *ia = a, *ib = b;
	return *ia - *ib;
}
