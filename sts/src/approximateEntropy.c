#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"  

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                A P P R O X I M A T E  E N T R O P Y   T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define TERMS	(1 << 14)

/*
 * One pass counts the cyclic (m+1)-bit patterns; the m-bit counts are their sums over the
 * last bit. The terms P*log(P/n) of the sums are memoized per count P for the current n,
 * and added in the original order.
 */
void
ApproximateEntropy(int m, int n)
{
	static unsigned int	*counts = NULL;
	static long			countsSize = 0;
	static double		*term = NULL;
	static int			termN = -1;
	int				i, r, blockSize, seqLength;
	double			sum, numOfBlocks, ApEn[2], apen, chi_squared, p_value;
	unsigned int	*P, c;
	long			size, v;
	
	fprintf(stats[TEST_APEN], "\t\t\tAPPROXIMATE ENTROPY TEST\n");
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_APEN], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_APEN], "\t\t(a) m (block length)    = %d\n", m);

	seqLength = n;
	r = 0;
	
	if ( (m < -1) || (m > 29) ) {
		fprintf(stats[TEST_APEN], "ApEn:  Insufficient memory available.\n");
		return;
	}
	/* counts of the (m+1)-bit patterns at 2^m.., then those of the m-bit patterns at 0.. */
	size = m >= 0 ? 3L << m : 0;
	if ( size > countsSize ) {
		free(counts);
		countsSize = size;
		if ( (counts = (unsigned int *) malloc(size*sizeof(unsigned int))) == NULL ) {
			countsSize = 0;
			fprintf(stats[TEST_APEN], "ApEn:  Insufficient memory available.\n");
			return;
		}
	}
	if ( term == NULL ) {
		if ( (term = (double *) malloc(TERMS*sizeof(double))) == NULL ) {
			fprintf(stats[TEST_APEN], "ApEn:  Insufficient memory available.\n");
			return;
		}
	}
	if ( termN != n ) {
		/* P*log(P/n) <= 0 for P <= n, so 1 marks a term not computed yet */
		for ( i=0; i<TERMS; i++ )
			term[i] = 1;
		termN = n;
	}
	if ( m >= 0 ) {
		cyclicWindowCounts(counts + (1L << m), m+1, n);
		for ( v=0; v<(1L << m); v++ )
			counts[v] = counts[(1L << m) + 2*v] + counts[(1L << m) + 2*v + 1];
	}
	for ( blockSize=m; blockSize<=m+1; blockSize++ ) {
		if ( blockSize == 0 ) {
			ApEn[0] = 0.00;
			r++;
		}
		else {
			numOfBlocks = (double)seqLength;
			P = blockSize == m ? counts : counts + (1L << m);
			sum = 0.0;
			for ( v=0; v<(blockSize > 0 ? 1L << blockSize : 0); v++ ) {
				c = P[v];
				if ( c > 0 ) {
					if ( c >= TERMS )
						sum += c*log(c/numOfBlocks);
					else {
						if ( term[c] > 0 )
							term[c] = c*log(c/numOfBlocks);
						sum += term[c];
					}
				}
			}
			sum /= numOfBlocks;
			ApEn[r] = sum;
			r++;
		}
	}
	apen = ApEn[0] - ApEn[1];
	
	chi_squared = 2.0*seqLength*(log(2) - apen);
	p_value = cephes_igamc(pow(2, m-1), chi_squared/2.0);
	
	fprintf(stats[TEST_APEN], "\t\t(b) n (sequence length) = %d\n", seqLength);
	fprintf(stats[TEST_APEN], "\t\t(c) Chi^2               = %f\n", chi_squared);
	fprintf(stats[TEST_APEN], "\t\t(d) Phi(m)	       = %f\n", ApEn[0]);
	fprintf(stats[TEST_APEN], "\t\t(e) Phi(m+1)	       = %f\n", ApEn[1]);
	fprintf(stats[TEST_APEN], "\t\t(f) ApEn                = %f\n", apen);
	fprintf(stats[TEST_APEN], "\t\t(g) Log(2)              = %f\n", log(2.0));
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");

	if ( m > (int)(log(seqLength)/log(2)-5) ) {
		fprintf(stats[TEST_APEN], "\t\tNote: The blockSize = %d exceeds recommended value of %d\n", m,
			MAX(1, (int)(log(seqLength)/log(2)-5)));
		fprintf(stats[TEST_APEN], "\t\tResults are inaccurate!\n");
		fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	}
	
	fprintf(stats[TEST_APEN], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
	fprintf(results[TEST_APEN], "%f\n", p_value);
}
