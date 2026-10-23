#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"  

/*
 * Berlekamp-Massey runs on 64 blocks at once. Word i of an array holds coefficient (or bit) i
 * of all 64 blocks, block g+r in bit 63-r. Each step updates every block with the same word
 * operations; per block it is the original algorithm:
 *   d = s[N] + sum_{i=1..L} C[i] s[N-i];  if d: C += x^(N-m) B, and if L <= N/2 then
 *   L = N+1-L, m = N, B = old C.
 * Two facts make that possible. deg C <= L at every step, so summing over all i gives the
 * same d. And with D = x^(N-m) B kept instead of B, every block shifts D by one each step:
 * D' = x (switch ? old C : D). The arrays never go past M+1 coefficients, so the original
 * code's arrays of M entries lose nothing that a later step reads.
 */

/* Element (r, c) of a 64 x 64 bit matrix is bit 63-c of a[r]; it moves to bit 63-r of a[c] */
static void
transpose64(uint64_t *a)
{
	int			j, k;
	uint64_t	m = 0x00000000FFFFFFFFULL, t;

	for ( j=32; j!=0; j>>=1, m^=m<<j ) {
		for ( k=0; k<64; k=(k+j+1)&~j ) {
			t = (a[k] ^ (a[k+j] >> j)) & m;
			a[k] ^= t;
			a[k+j] ^= t << j;
		}
	}
}

/* Linear complexity of blocks first .. first+count-1 (count <= 64) of M bits each */
static void
groupComplexity(int M, long first, int count, int *L, uint64_t *R, uint64_t *restrict C, uint64_t *B)
{
	uint64_t			a[64], plane[34], d, dn, sw, c, b, cn, carry, old;
	uint64_t *restrict	D;
	const uint64_t		*Rn;
	long				t0, t;
	long long			delta;
	int					r, i, k, N_, hi, np;

	np = 2;
	while ( ((1LL << (np-1)) - 1) < M )
		np++;
	/* R[M-1-t] = bit t of every block, so the window s[N-i], i = 0, 1, ... runs forward.
	   Unused lanes repeat the first block. */
	for ( t0=0; t0<M; t0+=64 ) {
		for ( r=0; r<64; r++ )
			a[r] = bitWindow(epsilonWords, (first + (r < count ? r : 0))*M + t0);
		transpose64(a);
		for ( t=t0; (t<t0+64) && (t<M); t++ )
			R[M-1-t] = a[t-t0];
	}
	memset(C, 0, (M+2)*sizeof(uint64_t));
	memset(B, 0, (M+4)*sizeof(uint64_t));
	memset(plane, 0, sizeof(plane));
	C[0] = ~0ULL;
	D = B + M + 2;
	D[1] = ~0ULL;
	hi = 1;
	d = R[M-1];
	for ( N_=0; N_<M; N_++ ) {
		/* plane[] holds N - 2L per block in two's complement; the top plane is its sign */
		sw = d & ~plane[np-1];
		Rn = R + (M - 2 - N_);
		dn = 0;
		for ( i=0; i<=hi; i++ ) {
			c = C[i];
			b = D[i];
			cn = c ^ (d & b);
			C[i] = cn;
			D[i] = b ^ ((c ^ b) & sw);
			dn ^= cn & Rn[i];
		}
		D--;
		carry = ~0ULL;
		for ( k=0; k<np; k++ ) {
			old = plane[k];
			plane[k] = ((old ^ carry) & ~sw) | (~old & sw);
			carry &= old;
		}
		hi++;
		while ( (hi > 0) && (C[hi] == 0) && (D[hi] == 0) )
			hi--;
		d = dn;
	}
	for ( r=0; r<count; r++ ) {
		delta = 0;
		for ( k=0; k<np; k++ )
			delta |= (long long)((plane[k] >> (63 - r)) & 1) << k;
		if ( delta & (1LL << (np-1)) )
			delta -= 1LL << np;
		L[r] = (int)((M - delta) / 2);
	}
}

void
LinearComplexity(int M, int n)
{
	static uint64_t	*work = NULL;
	static int		*complexity = NULL;
	static long		workSize = 0, complexitySize = 0;
	int       i, ii, N, L, parity, sign, K = 6;
	double    p_value, T_, mean, nu[7], chi2;
	double    pi[7] = { 0.01047, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833 };
	long      need;
	
	N = (int)floor(n/M);
	if ( M < 1 ) {
		printf("Insufficient Memory for Work Space:: Linear Complexity Test\n");
		return;
	}
	need = 3*(long)M + 7;
	if ( need > workSize ) {
		free(work);
		if ( (work = (uint64_t *) malloc(need*sizeof(uint64_t))) == NULL ) {
			workSize = 0;
			printf("Insufficient Memory for Work Space:: Linear Complexity Test\n");
			return;
		}
		workSize = need;
	}
	if ( N > complexitySize ) {
		free(complexity);
		if ( (complexity = (int *) malloc(N*sizeof(int))) == NULL ) {
			complexitySize = 0;
			printf("Insufficient Memory for Work Space:: Linear Complexity Test\n");
			return;
		}
		complexitySize = N;
	}

	fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "\tL I N E A R  C O M P L E X I T Y\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "\tM (substring length)     = %d\n", M);
	fprintf(stats[TEST_LINEARCOMPLEXITY], "\tN (number of substrings) = %d\n", N);
	fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "        F R E Q U E N C Y                            \n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "  C0   C1   C2   C3   C4   C5   C6    CHI2    P-value\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
	fprintf(stats[TEST_LINEARCOMPLEXITY], "\tNote: %d bits were discarded!\n", n%M);

	/* R: one zero word, then M words of block bits; C: M+2 words; B: M+4 words */
	work[0] = 0;
	for ( ii=0; ii<N; ii+=64 )
		groupComplexity(M, ii, N - ii < 64 ? N - ii : 64, complexity + ii,
						work + 1, work + M + 1, work + 2*(long)M + 3);

	if ( (parity = (M+1)%2) == 0 )
		sign = -1;
	else
		sign = 1;
	mean = M/2.0 + (9.0+sign)/36.0 - 1.0/pow(2, M) * (M/3.0 + 2.0/9.0);
	if ( (parity = M%2) == 0 )
		sign = 1;
	else
		sign = -1;
	for ( i=0; i<K+1; i++ )
		nu[i] = 0.00;
	for ( ii=0; ii<N; ii++ ) {
		L = complexity[ii];
		T_ = sign * (L - mean) + 2.0/9.0;
		
		if ( T_ <= -2.5 )
			nu[0]++;
		else if ( T_ > -2.5 && T_ <= -1.5 )
			nu[1]++;
		else if ( T_ > -1.5 && T_ <= -0.5 )
			nu[2]++;
		else if ( T_ > -0.5 && T_ <= 0.5 )
			nu[3]++;
		else if ( T_ > 0.5 && T_ <= 1.5 )
			nu[4]++;
		else if ( T_ > 1.5 && T_ <= 2.5 )
			nu[5]++;
		else
			nu[6]++;
	}
	chi2 = 0.00;
	for ( i=0; i<K+1; i++ ) 
		fprintf(stats[TEST_LINEARCOMPLEXITY], "%4d ", (int)nu[i]);
	for ( i=0; i<K+1; i++ )
		chi2 += pow(nu[i]-N*pi[i], 2) / (N*pi[i]);
	p_value = cephes_igamc(K/2.0, chi2/2.0);

	fprintf(stats[TEST_LINEARCOMPLEXITY], "%9.6f%9.6f\n", chi2, p_value);
	fprintf(results[TEST_LINEARCOMPLEXITY], "%f\n", p_value);
}
