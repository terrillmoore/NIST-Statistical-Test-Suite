#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                    B L O C K  F R E Q U E N C Y  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* Ones before bit p of the stream, for a p that never decreases between calls */
typedef struct {
	long	words;		/* whole words counted so far */
	long	ones;		/* ones in them */
} RunningCount;

static inline long
onesBefore(RunningCount *c, long p)
{
	while ( c->words < (p >> 6) )
		c->ones += popcount64(epsilonWords[c->words++]);
	if ( p & 63 )
		return c->ones + popcount64(epsilonWords[c->words] >> (64 - (p & 63)));
	return c->ones;
}

void
BlockFrequency(int M, int n)
{
	int				i, N, blockSum;
	long			before = 0, after;
	double			p_value, sum, pi, v, chi_squared, *square = NULL;
	RunningCount	count = { 0, 0 };
	
	N = n/M; 		/* # OF SUBSTRING BLOCKS      */
	sum = 0.0;
	
	/* (blockSum/M - 0.5)^2 depends on blockSum only; tabulate it when there are many blocks */
	if ( (N > 0) && (M <= N) )
		square = (double *)malloc((M+1)*sizeof(double));
	if ( square != NULL ) {
		for ( i=0; i<=M; i++ ) {
			pi = (double)i/(double)M;
			v = pi - 0.5;
			square[i] = v*v;
		}
		for ( i=0; i<N; i++ ) {
			after = onesBefore(&count, (long)(i+1)*M);
			sum += square[after - before];
			before = after;
		}
		free(square);
	}
	else {
		for ( i=0; i<N; i++ ) {
			after = onesBefore(&count, (long)(i+1)*M);
			blockSum = (int)(after - before);
			before = after;
			pi = (double)blockSum/(double)M;
			v = pi - 0.5;
			sum += v*v;
		}
	}
	chi_squared = 4.0 * M * sum;
	p_value = cephes_igamc(N/2.0, chi_squared/2.0);

	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t\tBLOCK FREQUENCY TEST\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(a) Chi^2           = %f\n", chi_squared);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(b) # of substrings = %d\n", N);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(c) block length    = %d\n", M);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(d) Note: %d bits were discarded.\n", n % M);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_BLOCK_FREQUENCY], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
	fprintf(results[TEST_BLOCK_FREQUENCY], "%f\n", p_value);
}
