#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"  

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
               O V E R L A P P I N G  T E M P L A T E  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
double	Pr(int u, double eta);

static uint64_t	*runBuffer = NULL, *hitBuffer = NULL;
static long		bufferWords = 0;

/* dst[t] = a[t] & the 64 bits of b that start at bit 64t+s, for t < words */
static void
andShifted(uint64_t *dst, const uint64_t *a, const uint64_t *b, long words, long s)
{
	const uint64_t	*p = b + (s >> 6);
	int				r = (int)(s & 63);
	long			t;

	if ( r == 0 )
		for ( t=0; t<words; t++ )
			dst[t] = a[t] & p[t];
	else
		for ( t=0; t<words; t++ )
			dst[t] = a[t] & ((p[t] << r) | (p[t+1] >> (64 - r)));
}

/* Packed array with bit p set when bits p .. p+m-1 of the stream are all ones (m >= 1) */
static const uint64_t *
onesWindows(int n, int m)
{
	long			words = ((long)n + 63) / 64, need = words + m / 64 + 3, step, covered = 0;
	const uint64_t	*run = epsilonWords, *hit = NULL;

	if ( need > bufferWords ) {
		free(runBuffer);
		free(hitBuffer);
		runBuffer = (uint64_t *)malloc(need * sizeof(uint64_t));
		hitBuffer = (uint64_t *)malloc(need * sizeof(uint64_t));
		if ( (runBuffer == NULL) || (hitBuffer == NULL) ) {
			bufferWords = 0;
			return NULL;
		}
		bufferWords = need;
	}
	memset(runBuffer + words, 0, (need - words) * sizeof(uint64_t));
	/* run holds starts of runs of length step, hit of length covered */
	for ( step=1; step<=m; step*=2 ) {
		if ( m & step ) {
			if ( hit == NULL ) {
				if ( run == epsilonWords )
					hit = run;
				else {
					memcpy(hitBuffer, run, words * sizeof(uint64_t));
					hit = hitBuffer;
				}
			}
			else {
				andShifted(hitBuffer, hit, run, words, covered);
				hit = hitBuffer;
			}
			covered += step;
		}
		if ( 2*step <= m ) {
			andShifted(runBuffer, run, run, words, step);
			run = runBuffer;
		}
	}
	return hit;
}

void
OverlappingTemplateMatchings(int m, int n)
{
	int				i;
	double			W_obs, eta, sum, chi2, p_value, lambda;
	int				M, N, K = 5;
	unsigned int	nu[6] = { 0, 0, 0, 0, 0, 0 };
	//double			pi[6] = { 0.143783, 0.139430, 0.137319, 0.124314, 0.106209, 0.348945 };
	double			pi[6] = { 0.364091, 0.185659, 0.139381, 0.100571, 0.0704323, 0.139865 };
	const uint64_t	*hit = NULL;

	M = 1032;
	N = n/M;
	
	if ( (m < 0) || ((m >= 1) && (m <= M) && (N > 0) && ((hit = onesWindows(n, m)) == NULL)) ) {
		fprintf(stats[TEST_OVERLAPPING], "\t\t    OVERLAPPING TEMPLATE OF ALL ONES TEST\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\tTEMPLATE DEFINITION:  Insufficient memory, Overlapping Template Matchings test aborted!\n");
	}
	
	lambda = (double)(M-m+1)/pow(2,m);
	eta = lambda/2.0;
	sum = 0.0;
	for ( i=0; i<K; i++ ) {			/* Compute Probabilities */
		pi[i] = Pr(i, eta);
		sum += pi[i];
	}
	pi[K] = 1 - sum;

	for ( i=0; i<N; i++ ) {
		if ( m < 1 )
			W_obs = M-m+1;
		else if ( (m > M) || (hit == NULL) )
			W_obs = 0;
		else
			W_obs = rangeOnes(hit, (long)i*M, M-m+1);
		if ( W_obs <= 4 )
			nu[(int)W_obs]++;
		else
			nu[K]++;
	}
	sum = 0;
	chi2 = 0.0;                                   /* Compute Chi Square */
	for ( i=0; i<K+1; i++ ) {
		chi2 += pow((double)nu[i] - (double)N*pi[i], 2)/((double)N*pi[i]);
		sum += nu[i];
	}
	p_value = cephes_igamc(K/2.0, chi2/2.0);

	fprintf(stats[TEST_OVERLAPPING], "\t\t    OVERLAPPING TEMPLATE OF ALL ONES TEST\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t(a) n (sequence_length)      = %d\n", n);
	fprintf(stats[TEST_OVERLAPPING], "\t\t(b) m (block length of 1s)   = %d\n", m);
	fprintf(stats[TEST_OVERLAPPING], "\t\t(c) M (length of substring)  = %d\n", M);
	fprintf(stats[TEST_OVERLAPPING], "\t\t(d) N (number of substrings) = %d\n", N);
	fprintf(stats[TEST_OVERLAPPING], "\t\t(e) lambda [(M-m+1)/2^m]     = %f\n", lambda);
	fprintf(stats[TEST_OVERLAPPING], "\t\t(f) eta                      = %f\n", eta);
	fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t   F R E Q U E N C Y\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t  0   1   2   3   4 >=5   Chi^2   P-value  Assignment\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
	fprintf(stats[TEST_OVERLAPPING], "\t\t%3d %3d %3d %3d %3d %3d  %f ",
		nu[0], nu[1], nu[2], nu[3], nu[4], nu[5], chi2);

	if ( isNegative(p_value) || isGreaterThanOne(p_value) )
		fprintf(stats[TEST_OVERLAPPING], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

	fprintf(stats[TEST_OVERLAPPING], "%f %s\n\n", p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS");
	fprintf(results[TEST_OVERLAPPING], "%f\n", p_value);
}

double
Pr(int u, double eta)
{
	int		l;
	double	sum, p;
	
	if ( u == 0 )
		p = exp(-eta);
	else {
		sum = 0.0;
		for ( l=1; l<=u; l++ )
			sum += exp(-eta-u*log(2)+l*log(eta)-cephes_lgam(l+1)+cephes_lgam(u)-cephes_lgam(l)-cephes_lgam(u-l+1));
		p = sum;
	}
	return p;
}
