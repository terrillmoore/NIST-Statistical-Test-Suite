#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                              R U N S  T E S T 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* Number of k in 1..n-1 with epsilon[k] != epsilon[k-1] */
static long
transitions(int n)
{
	const uint64_t	*w = epsilonWords;
	long			k, words = ((long)n + 63) / 64, count = 0;

	/* bit i of x ^ (x shifted by one) compares bits i and i+1; past bit n-1 the words are zero */
	for ( k=0; k<words; k++ )
		count += popcount64(w[k] ^ ((w[k] << 1) | (w[k+1] >> 63)));
	return count - (long)((w[(n-1) >> 6] >> (63 - (n-1) % 64)) & 1);
}

void
Runs(int n)
{
	int		S;
	double	pi, V, erfc_arg, p_value;

	S = sequenceOnes();
	pi = (double)S / (double)n;

	if ( fabs(pi - 0.5) > (2.0 / sqrt(n)) ) {
		fprintf(stats[TEST_RUNS], "\t\t\t\tRUNS TEST\n");
		fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
		fprintf(stats[TEST_RUNS], "\t\tPI ESTIMATOR CRITERIA NOT MET! PI = %f\n", pi);
		p_value = 0.0;
	}
	else {

		V = (double)(1 + transitions(n));
	
		erfc_arg = fabs(V - 2.0 * n * pi * (1-pi)) / (2.0 * pi * (1-pi) * sqrt(2*n));
		p_value = erfc(erfc_arg);
		
		fprintf(stats[TEST_RUNS], "\t\t\t\tRUNS TEST\n");
		fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
		fprintf(stats[TEST_RUNS], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
		fprintf(stats[TEST_RUNS], "\t\t(a) Pi                        = %f\n", pi);
		fprintf(stats[TEST_RUNS], "\t\t(b) V_n_obs (Total # of runs) = %d\n", (int)V);
		fprintf(stats[TEST_RUNS], "\t\t(c) V_n_obs - 2 n pi (1-pi)\n");
		fprintf(stats[TEST_RUNS], "\t\t    -----------------------   = %f\n", erfc_arg);
		fprintf(stats[TEST_RUNS], "\t\t      2 sqrt(2n) pi (1-pi)\n");
		fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
		if ( isNegative(p_value) || isGreaterThanOne(p_value) )
			fprintf(stats[TEST_RUNS], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

		fprintf(stats[TEST_RUNS], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
	}

	fprintf(results[TEST_RUNS], "%f\n", p_value);
}
