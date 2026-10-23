#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
         D I S C R E T E  F O U R I E R  T R A N S F O R M  T E S T 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void	__ogg_fdrffti(int n, double *wsave, int *ifac);
double	*__ogg_fdrfftfbits(int n, const uint64_t *bits, double *r, double *wsave, int *ifac);

/*
 * The work arrays and the factors and twiddles of the transform depend only on n, so they
 * are set up once. sqrt is correctly rounded and never decreasing, so sqrt(s) < bound holds
 * exactly when s is below the least s0 with sqrt(s0) >= bound; the count compares the
 * squared magnitudes with s0.
 */
void
DiscreteFourierTransform(int n)
{
	static double	*X = NULL, *wsave = NULL, limit;
	static int		ifac[15], prepared = 0;
	double	p_value, upperBound, percentile, N_l, N_o, d, s, *Y;
	int		i, count;

	upperBound = sqrt(2.995732274*n);
	if ( n != prepared ) {
		free(X);
		free(wsave);
		X = (double*) malloc(n*sizeof(double));
		wsave = (double *)calloc(2*(size_t)n, sizeof(double));
		if ( (X == NULL) || (wsave == NULL) ) {
			free(X);
			free(wsave);
			X = wsave = NULL;
			prepared = 0;
			fprintf(stats[7],"\t\tUnable to allocate working arrays for the DFT.\n");
			return;
		}
		__ogg_fdrffti(n, wsave, ifac);		/* INITIALIZE WORK ARRAYS */
		s = upperBound*upperBound;
		while ( sqrt(s) < upperBound )
			s = nextafter(s, HUGE_VAL);
		while ( sqrt(nextafter(s, 0)) >= upperBound )
			s = nextafter(s, 0);
		limit = s;
		prepared = n;
	}
	
	Y = __ogg_fdrfftfbits(n, epsilonWords, X, wsave, ifac);	/* APPLY FORWARD FFT */
	
	/* m[0] = |Y[0]|, m[i] = sqrt(Y[2i-1]^2 + Y[2i]^2); count the m[i] < upperBound, i < n/2 */
	count = 0;				       /* CONFIDENCE INTERVAL */
	if ( (n/2 > 0) && (Y[0]*Y[0] < limit) )
		count++;
	for ( i=1; i<n/2; i++ )
		if ( Y[2*i-1]*Y[2*i-1] + Y[2*i]*Y[2*i] < limit )
			count++;
	percentile = (double)count/(n/2)*100;
	N_l = (double) count;       /* number of peaks less than h = sqrt(3*n) */
	N_o = (double) 0.95*n/2.0;
	d = (N_l - N_o)/sqrt(n/4.0*0.95*0.05);
	p_value = erfc(fabs(d)/sqrt(2.0));

	fprintf(stats[TEST_FFT], "\t\t\t\tFFT TEST\n");
	fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_FFT], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_FFT], "\t\t(a) Percentile = %f\n", percentile);
	fprintf(stats[TEST_FFT], "\t\t(b) N_l        = %f\n", N_l);
	fprintf(stats[TEST_FFT], "\t\t(c) N_o        = %f\n", N_o);
	fprintf(stats[TEST_FFT], "\t\t(d) d          = %f\n", d);
	fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");

	fprintf(stats[TEST_FFT], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
	fprintf(results[TEST_FFT], "%f\n", p_value);
}
