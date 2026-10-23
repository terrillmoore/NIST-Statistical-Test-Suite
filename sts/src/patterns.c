/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
              W I N D O W  C O U N T S  O F  T H E  S T R E A M
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <string.h>
#include "../include/externs.h"
#include "../include/utilities.h"

/* Adds the m-bit windows starting at start .. start+count-1 to counts (1 <= m <= 32) */
void
windowCounts(unsigned int *counts, long start, long count, int m)
{
	int			j, chunk = 65 - m;
	uint64_t	x;

	while ( count > 0 ) {
		if ( chunk > count )
			chunk = (int)count;
		x = bitWindow(epsilonWords, start);
		for ( j=0; j<chunk; j++ ) {
			counts[x >> (64 - m)]++;
			x <<= 1;
		}
		start += chunk;
		count -= chunk;
	}
}

/* counts[v] = number of i in 0..n-1 whose m bits i, i+1, ... (taken mod n) read v, first bit
   on top (1 <= m <= 32) */
void
cyclicWindowCounts(unsigned int *counts, int m, int n)
{
	uint64_t		tail;
	long			i, pos;
	unsigned int	v;
	int				j;

	memset(counts, 0, ((size_t)1 << m)*sizeof(unsigned int));
	if ( n < m ) {
		for ( i=0; i<n; i++ ) {
			v = 0;
			for ( j=0; j<m; j++ ) {
				pos = (i + j) % n;
				v = (v << 1) | (int)((epsilonWords[pos >> 6] >> (63 - (pos & 63))) & 1);
			}
			counts[v]++;
		}
		return;
	}
	windowCounts(counts, 0, (long)n - m + 1, m);
	if ( m > 1 ) {
		/* the last m-1 windows wrap: bits n-m+1 .. n-1, then 0 .. m-2 */
		tail = bitWindow(epsilonWords, (long)n - m + 1) >> (65 - m) << (65 - m) | epsilonWords[0] >> (m - 1);
		for ( j=0; j<m-1; j++ )
			counts[(tail << j) >> (64 - m)]++;
	}
}
