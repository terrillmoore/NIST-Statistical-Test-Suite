#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"  

double psi2(int m, int n);

/* Counts of the cyclic windows of the current stream, for window length countsLength */
static unsigned int	*counts = NULL;
static long			countsSize = 0;
static int			countsLength = -1;

void
Serial(int m, int n)
{
	double	p_value1, p_value2, psim0, psim1, psim2, del1, del2;
	
	countsLength = -1;
	psim0 = psi2(m, n);
	psim1 = psi2(m-1, n);
	psim2 = psi2(m-2, n);
	del1 = psim0 - psim1;
	del2 = psim0 - 2.0*psim1 + psim2;
	p_value1 = cephes_igamc(pow(2, m-1)/2, del1/2.0);
	p_value2 = cephes_igamc(pow(2, m-2)/2, del2/2.0);
	
	fprintf(stats[TEST_SERIAL], "\t\t\t       SERIAL TEST\n");
	fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_SERIAL], "\t\t COMPUTATIONAL INFORMATION:		  \n");
	fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_SERIAL], "\t\t(a) Block length    (m) = %d\n", m);
	fprintf(stats[TEST_SERIAL], "\t\t(b) Sequence length (n) = %d\n", n);
	fprintf(stats[TEST_SERIAL], "\t\t(c) Psi_m               = %f\n", psim0);
	fprintf(stats[TEST_SERIAL], "\t\t(d) Psi_m-1             = %f\n", psim1);
	fprintf(stats[TEST_SERIAL], "\t\t(e) Psi_m-2             = %f\n", psim2);
	fprintf(stats[TEST_SERIAL], "\t\t(f) Del_1               = %f\n", del1);
	fprintf(stats[TEST_SERIAL], "\t\t(g) Del_2               = %f\n", del2);
	fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_SERIAL], "%s\t\tp_value1 = %f\n", p_value1 < ALPHA ? "FAILURE" : "SUCCESS", p_value1);
	fprintf(results[TEST_SERIAL], "%f\n", p_value1);

	fprintf(stats[TEST_SERIAL], "%s\t\tp_value2 = %f\n\n", p_value2 < ALPHA ? "FAILURE" : "SUCCESS", p_value2);
	fprintf(results[TEST_SERIAL], "%f\n", p_value2);
}

/*
 * The counts of the m-bit patterns are those of the (m+1)-bit patterns summed over the last
 * bit, so after the first call the next shorter length needs no pass over the stream.
 * The squares are exact integers: while their total stays below 2^53 every partial sum of
 * the original double loop is exact too, so the integer total converted once is that sum.
 */
double
psi2(int m, int n)
{
	long		i;
	uint64_t	total;
	double		sum;
	
	if ( (m == 0) || (m == -1) )
		return 0.0;
	if ( (m >= 1) && (countsLength == m + 1) ) {
		for ( i=0; i<(1L << m); i++ )
			counts[i] = counts[2*i] + counts[2*i+1];
	}
	else {
		if ( (m >= 1) && (m <= 30) && ((1L << m) > countsSize) ) {
			free(counts);
			countsSize = 1L << m;
			if ( (counts = (unsigned int *) malloc(countsSize*sizeof(unsigned int))) == NULL )
				countsSize = 0;
		}
		if ( (m < 1) || (m > 30) || (counts == NULL) ) {
			fprintf(stats[TEST_SERIAL], "Serial Test:  Insufficient memory available.\n");
			countsLength = -1;
			return 0.0;
		}
		cyclicWindowCounts(counts, m, n);
	}
	countsLength = m;
	total = 0;
	for ( i=0; i<(1L << m); i++ )
		total += (uint64_t)counts[i] * counts[i];
	if ( total < (1ULL << 53) )
		sum = (double)total;
	else {
		sum = 0.0;
		for ( i=0; i<(1L << m); i++ )
			sum += pow(counts[i], 2);
	}
	sum = (sum * pow(2, m)/(double)n) - (double)n;
	
	return sum;
}
