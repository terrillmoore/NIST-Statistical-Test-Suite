/* Notes from RFB: 

Looks like the user-level routines are:

Real FFT

void __ogg_fdrffti(int n, double *wsave, int *ifac)
double *__ogg_fdrfftfbits(int n, const uint64_t *bits, double *r, double *wsave, int *ifac)

__ogg_fdrffti == initialization
__ogg_fdrfftfbits == forward transform of the +1/-1 values of packed bits

Parameters are
n == length of sequence
bits == sequence to be transformed (input)
r == work array of length n (allocated by caller)
wsave == work array of length 2n (allocated by caller)
returned == r or wsave, whichever holds the transformed sequence (output)
ifac == work array of length 15 (allocated by caller)
*/

/********************************************************************
*                                                                  *
* THIS FILE IS PART OF THE OggSQUISH SOFTWARE CODEC SOURCE CODE.   *
*                                                                  *
********************************************************************

file: fft.c
function: Fast discrete Fourier and cosine transforms and inverses
author: Monty <xiphmont@mit.edu>
modifications by: Monty
last modification date: Jul 1 1996

********************************************************************/

/* These Fourier routines were originally based on the Fourier
routines of the same names from the NETLIB bihar and fftpack
fortran libraries developed by Paul N. Swarztrauber at the National
Center for Atmospheric Research in Boulder, CO USA.  They have been
reimplemented in C and optimized in a few ways for OggSquish. */

/* As the original fortran libraries are public domain, the C Fourier
routines in this file are hereby released to the public domain as
well.  The C routines here produce output exactly equivalent to the
original fortran routines.  Of particular interest are the facts
that (like the original fortran), these routines can work on
arbitrary length vectors that need not be powers of two in
length. */

#include <math.h>
#include <stdint.h>
#define STIN static

static void drfti1(int n, double *wa, int *ifac)
{
	static int		ntryh[4] = { 4,2,3,5 };
	static double	tpi = 6.28318530717958647692528676655900577;
	double			arg, argh, argld, fi;
	int				ntry=0,i,j=-1;
	int				k1, l1, l2, ib;
	int				ld, ii, ip, is, nq, nr;
	int				ido, ipm, nfm1;
	int				nl=n;
	int				nf=0;
	int				e, e0, ido0;
	long			f0, i0;
	
L101:
	j++;
	if ( j < 4 )
		ntry = ntryh[j];
	else
		ntry += 2;
	
L104:
	nq = nl/ntry;
	nr = nl-ntry*nq;
	if ( nr != 0 )
		goto L101;
	
	nf++;
	ifac[nf+1] = ntry;
	nl = nq;
	if ( ntry != 2 )
		goto L107;
	if ( nf == 1 )
		goto L107;
	
	for ( i=1; i<nf; i++ ) {
		ib = nf-i+1;
		ifac[ib+1] = ifac[ib];
	}
	ifac[2] = 2;
	
L107:
	if( nl != 1 )
		goto L104;
	ifac[0] = n;
	ifac[1] = nf;
	argh = tpi/n;
	is = 0;
	nfm1 = nf-1;
	l1 = 1;
	
	if ( nfm1 == 0 )
		return;
	
	/*
	 * With ld = 2^e o, o odd, (double)ld*argh is exactly 2^e times (double)o*argh, so
	 * fi*argld rounds the same product as f0*argld0 for ld0 = 2^e0 o and f0 = fi 2^(e-e0).
	 * When ld0 is a row of the first factor (l1 = 1, ld0 < ifac[2]) and f0 is in its range,
	 * the entry is copied from that row once the row is filled.
	 */
	ido0 = n/ifac[2];
	for ( k1=0; k1<nfm1; k1++ ) {
		ip = ifac[k1+2];
		ld = 0;
		l2 = l1*ip;
		ido = n/l2;
		ipm = ip-1;
		
		for ( j=0; j<ipm; j++ ) {
			ld += l1;
			i = is;
			argld = (double)ld*argh;
			for ( e=0; ((ld >> e) & 1) == 0; e++ )
				;
			for ( e0=e; e0>=0; e0-- )
				if ( ((ld >> e << e0) < ifac[2]) && ((k1 > 0) || ((ld >> e << e0) < ld)) )
					break;
			fi = 0.0;
			for ( ii=2; ii<ido; ii+=2 ) {
				fi += 1.0;
				f0 = (long)(ii/2) << (e - e0);
				if ( (e0 >= 0) && (f0 <= (ido0-1)/2) ) {
					i0 = ((long)(ld >> e << e0) - 1)*ido0 + 2*(f0-1);
					wa[i++] = wa[i0];
					wa[i++] = wa[i0+1];
				}
				else {
					arg = fi*argld;
					wa[i++] = cos(arg);
					wa[i++] = sin(arg);
				}
			}
			is += ido;
		}
		l1 = l2;
	}
}

void __ogg_fdrffti(int n, double *wsave, int *ifac)
{	
	if ( n == 1 )
		return;
	drfti1(n, wsave+n, ifac);
}


STIN void dradf2(int ido, int l1, double *cc, double *ch, double *wa1)
{
	int		i , k;
	double	ti2 , tr2;
	int		t0 , t1, t2, t3, t4, t5, t6;
	
	t1 = 0;
	t0 = (t2 = l1*ido);
	t3 = ido<<1;
	for ( k=0; k<l1; k++ ) {
		ch[t1<<1] = cc[t1]+cc[t2];
		ch[(t1<<1)+t3-1] = cc[t1]-cc[t2];
		t1 += ido;
		t2 += ido;
	}
	
	if ( ido < 2 )
		return;
	if ( ido == 2 )
		goto L105;
	
	t1 = 0;
	t2 = t0;
	for ( k=0; k<l1; k++ ) {
		t3 = t2;
		t4 = (t1<<1)+(ido<<1);
		t5 = t1;
		t6 = t1+t1;
		for ( i=2; i<ido; i+=2 ) {
			t3 += 2;
			t4 -= 2;
			t5 += 2;
			t6 += 2;
			tr2 = wa1[i-2]*cc[t3-1]+wa1[i-1]*cc[t3];
			ti2 = wa1[i-2]*cc[t3]-wa1[i-1]*cc[t3-1];
			ch[t6] = cc[t5]+ti2;
			ch[t4] = ti2-cc[t5];
			ch[t6-1] = cc[t5-1]+tr2;
			ch[t4-1] = cc[t5-1]-tr2;
		}
		t1 += ido;
		t2 += ido;
	}
	
	if ( ido%2 == 1 )
		return;
	
L105:
	t3 = (t2 = (t1 = ido)-1);
	t2 += t0;
	for ( k=0; k<l1; k++ ) {
		ch[t1] = -cc[t2];
		ch[t1-1] = cc[t3];
		t1 += ido<<1;
		t2 += ido;
		t3 += ido;
	}
}

STIN void dradf4(int ido, int l1, double *cc, double *ch, double *wa1,
				 double *wa2, double *wa3)
{
	static double	hsqt2 = .70710678118654752440084436210485;
	int				i, k, t0, t1, t2, t3, t4, t5, t6;
	double			ci2, ci3, ci4, cr2, cr3, cr4;
	double			ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
	
	t0 = l1*ido;
	t1 = t0;
	t4 = t1<<1;
	t2 = t1+(t1<<1);
	t3 = 0;
	
	for ( k=0; k<l1; k++ ) {
		tr1 = cc[t1]+cc[t2];
		tr2 = cc[t3]+cc[t4];
		ch[t5 = t3<<2] = tr1+tr2;
		ch[(ido<<2)+t5-1] = tr2-tr1;
		ch[(t5 += (ido<<1))-1] = cc[t3]-cc[t4];
		ch[t5] = cc[t2]-cc[t1];
		
		t1 += ido;
		t2 += ido;
		t3 += ido;
		t4 += ido;
	}
	
	if ( ido < 2 )
		return;
	if ( ido == 2 )
		goto L105;
	
	t1 = 0;
	for ( k=0; k<l1; k++ ) {
		t2 = t1;
		t4 = t1<<2;
		t5 = (t6 = ido<<1)+t4;
		for ( i=2; i<ido; i+=2 ) {
			t3 = (t2 += 2);
			t4 += 2;
			t5 -= 2;
			
			t3 += t0;
			cr2 = wa1[i-2]*cc[t3-1]+wa1[i-1]*cc[t3];
			ci2 = wa1[i-2]*cc[t3]-wa1[i-1]*cc[t3-1];
			t3 += t0;
			cr3 = wa2[i-2]*cc[t3-1]+wa2[i-1]*cc[t3];
			ci3 = wa2[i-2]*cc[t3]-wa2[i-1]*cc[t3-1];
			t3 += t0;
			cr4 = wa3[i-2]*cc[t3-1]+wa3[i-1]*cc[t3];
			ci4 = wa3[i-2]*cc[t3]-wa3[i-1]*cc[t3-1];
			
			tr1 = cr2+cr4;
			tr4 = cr4-cr2;
			ti1 = ci2+ci4;
			ti4 = ci2-ci4;
			ti2 = cc[t2]+ci3;
			ti3 = cc[t2]-ci3;
			tr2 = cc[t2-1]+cr3;
			tr3 = cc[t2-1]-cr3;
			
			
			ch[t4-1] = tr1+tr2;
			ch[t4] = ti1+ti2;
			
			ch[t5-1] = tr3-ti4;
			ch[t5] = tr4-ti3;
			
			ch[t4+t6-1] = ti4+tr3;
			ch[t4+t6] = tr4+ti3;
			
			ch[t5+t6-1] = tr2-tr1;
			ch[t5+t6] = ti1-ti2;
		}
		t1 += ido;
	}
	if ( ido%2 == 1 )
		return;
	
L105:
		
	t2 = (t1 = t0+ido-1)+(t0<<1);
	t3 = ido<<2;
	t4 = ido;
	t5 = ido<<1;
	t6 = ido;
	
	for ( k=0; k<l1; k++ ) {
		ti1 = -hsqt2*(cc[t1]+cc[t2]);
		tr1 = hsqt2*(cc[t1]-cc[t2]);
		ch[t4-1] = tr1+cc[t6-1];
		ch[t4+t5-1] = cc[t6-1]-tr1;
		ch[t4] = ti1-cc[t1+t0];
		ch[t4+t5] = ti1+cc[t1+t0];
		t1 += ido;
		t2 += ido;
		t4 += t3;
		t6 += ido;
	}
}

STIN void dradfg(int ido, int ip, int l1, int idl1, double *cc, double *c1,
				 double *c2, double *ch, double *ch2, double *wa)
{
	static double	tpi = 6.28318530717958647692528676655900577;
	int				idij, ipph, i, j, k, l, ic, ik, is;
	int				t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
	double			dc2, ai1, ai2, ar1, ar2, ds2;
	int				nbd;
	double			dcp, arg, dsp, ar1h, ar2h;
	int				idp2, ipp2;
	
	arg = tpi/(double)ip;
	dcp = cos(arg);
	dsp = sin(arg);
	ipph = (ip+1)>>1;
	ipp2 = ip;
	idp2 = ido;
	nbd = (ido-1)>>1;
	t0 = l1*ido;
	t10 = ip*ido;
	
	if ( ido == 1 )
		goto L119;
	for ( ik=0; ik<idl1; ik++ )
		ch2[ik] = c2[ik];
	
	t1 = 0;
	for ( j=1; j<ip; j++ ) {
		t1 += t0;
		t2 = t1;
		for ( k=0; k<l1; k++ ) {
			ch[t2] = c1[t2];
			t2 += ido;
		}
	}
	
	is = -ido;
	t1 = 0;
	if ( nbd > l1 ) {
		for ( j=1; j<ip; j++ ) {
			t1 += t0;
			is += ido;
			t2 = -ido+t1;
			for ( k=0; k<l1; k++ ) {
				idij = is-1;
				t2 += ido;
				t3 = t2;
				for ( i=2; i<ido; i+=2 ) {
					idij += 2;
					t3 += 2;
					ch[t3-1] = wa[idij-1]*c1[t3-1]+wa[idij]*c1[t3];
					ch[t3] = wa[idij-1]*c1[t3]-wa[idij]*c1[t3-1];
				}
			}
		}
	}
	else {
		for ( j=1; j<ip; j++ ) {
			is += ido;
			idij = is-1;
			t1 += t0;
			t2 = t1;
			for ( i=2; i<ido; i+=2 ) {
				idij += 2;
				t2 += 2;
				t3 = t2;
				for ( k=0; k<l1; k++ ) {
					ch[t3-1] = wa[idij-1]*c1[t3-1]+wa[idij]*c1[t3];
					ch[t3] = wa[idij-1]*c1[t3]-wa[idij]*c1[t3-1];
					t3 += ido;
				}
			}
		}
	}
	
	t1 = 0;
	t2 = ipp2*t0;
	if( nbd < l1 ) {
		for ( j=1; j<ipph; j++ ) {
			t1 += t0;
			t2 -= t0;
			t3 = t1;
			t4 = t2;
			for ( i=2; i<ido; i+=2 ) {
				t3 += 2;
				t4 += 2;
				t5 = t3-ido;
				t6 = t4-ido;
				for ( k=0; k<l1; k++ ) {
					t5 += ido;
					t6 += ido;
					c1[t5-1] = ch[t5-1]+ch[t6-1];
					c1[t6-1] = ch[t5]-ch[t6];
					c1[t5] = ch[t5]+ch[t6];
					c1[t6] = ch[t6-1]-ch[t5-1];
				}
			}
		}
	}
	else {
		for ( j=1; j<ipph; j++ ) {
			t1 += t0;
			t2 -= t0;
			t3 = t1;
			t4 = t2;
			for ( k=0; k<l1; k++ ) {
				t5 = t3;
				t6 = t4;
				for ( i=2; i<ido; i+=2 ) {
					t5 += 2;
					t6 += 2;
					c1[t5-1] = ch[t5-1]+ch[t6-1];
					c1[t6-1] = ch[t5]-ch[t6];
					c1[t5] = ch[t5]+ch[t6];
					c1[t6] = ch[t6-1]-ch[t5-1];
				}
				t3 += ido;
				t4 += ido;
			}
		}
	}
	
L119:
	for ( ik=0; ik<idl1; ik++ )
		c2[ik] = ch2[ik];
	
	t1 = 0;
	t2 = ipp2*idl1;
	for ( j=1; j<ipph; j++ ) {
		t1 += t0;
		t2 -= t0;
		t3 = t1-ido;
		t4 = t2-ido;
		for ( k=0; k<l1; k++ ) {
			t3 += ido;
			t4 += ido;
			c1[t3] = ch[t3]+ch[t4];
			c1[t4] = ch[t4]-ch[t3];
		}
	}
	
	ar1 = 1.0;
	ai1 = 0.0;
	t1 = 0;
	t2 = ipp2*idl1;
	t3 = (ip-1)*idl1;
	for ( l=1; l<ipph; l++ ) {
		t1 += idl1;
		t2 -= idl1;
		ar1h = dcp*ar1-dsp*ai1;
		ai1 = dcp*ai1+dsp*ar1;
		ar1 = ar1h;
		t4 = t1;
		t5 = t2;
		t6 = t3;
		t7 = idl1;
		
		for ( ik=0; ik<idl1; ik++ ) {
			ch2[t4++] = c2[ik]+ar1*c2[t7++];
			ch2[t5++] = ai1*c2[t6++];
		}
		
		dc2 = ar1;
		ds2 = ai1;
		ar2 = ar1;
		ai2 = ai1;
		
		t4 = idl1;
		t5 = (ipp2-1)*idl1;
		for ( j=2; j<ipph; j++ ) {
			t4 += idl1;
			t5 -= idl1;
			
			ar2h = dc2*ar2-ds2*ai2;
			ai2 = dc2*ai2+ds2*ar2;
			ar2 = ar2h;
			
			t6 = t1;
			t7 = t2;
			t8 = t4;
			t9 = t5;
			for ( ik=0; ik<idl1; ik++ ) {
				ch2[t6++] += ar2*c2[t8++];
				ch2[t7++] += ai2*c2[t9++];
			}
		}
	}
	
	t1 = 0;
	for ( j=1; j<ipph; j++ ) {
		t1 += idl1;
		t2 = t1;
		for ( ik=0; ik<idl1; ik++ )
			ch2[ik] += c2[t2++];
	}
	
	if ( ido < l1 )
		goto L132;
	
	t1 = 0;
	t2 = 0;
	for ( k=0; k<l1; k++ ) {
		t3 = t1;
		t4 = t2;
		for ( i=0; i<ido; i++ )
			cc[t4++] = ch[t3++];
		t1 += ido;
		t2 += t10;
	}
	
	goto L135;
	
L132:
	for ( i=0; i<ido; i++ ) {
		t1 = i;
		t2 = i;
		for ( k=0; k<l1; k++ ) {
			cc[t2] = ch[t1];
			t1 += ido;
			t2 += t10;
		}
	}
	
L135:
	t1 = 0;
	t2 = ido<<1;
	t3 = 0;
	t4 = ipp2*t0;
	for ( j=1; j<ipph; j++ ) {
		
		t1 += t2;
		t3 += t0;
		t4 -= t0;
		
		t5 = t1;
		t6 = t3;
		t7 = t4;
		
		for ( k=0; k<l1; k++ ) {
			cc[t5-1] = ch[t6];
			cc[t5] = ch[t7];
			t5 += t10;
			t6 += ido;
			t7 += ido;
		}
	}
	
	if ( ido == 1 )
		return;
	if ( nbd < l1 )
		goto L141;
	
	t1 = -ido;
	t3 = 0;
	t4 = 0;
	t5 = ipp2*t0;
	for ( j=1; j<ipph; j++ ) {
		t1 += t2;
		t3 += t2;
		t4 += t0;
		t5 -= t0;
		t6 = t1;
		t7 = t3;
		t8 = t4;
		t9 = t5;
		for ( k=0; k<l1; k++ ) {
			for ( i=2; i<ido; i+=2 ) {
				ic = idp2-i;
				cc[i+t7-1] = ch[i+t8-1]+ch[i+t9-1];
				cc[ic+t6-1] = ch[i+t8-1]-ch[i+t9-1];
				cc[i+t7] = ch[i+t8]+ch[i+t9];
				cc[ic+t6] = ch[i+t9]-ch[i+t8];
			}
			t6 += t10;
			t7 += t10;
			t8 += ido;
			t9 += ido;
		}
	}
	return;
	
L141:
		
	t1 = -ido;
	t3 = 0;
	t4 = 0;
	t5 = ipp2*t0;
	for ( j=1; j<ipph; j++ ) {
		t1 += t2;
		t3 += t2;
		t4 += t0;
		t5 -= t0;
		for ( i=2; i<ido; i+=2 ) {
			t6 = idp2+t1-i;
			t7 = i+t3;
			t8 = i+t4;
			t9 = i+t5;
			for ( k=0; k<l1; k++ ) {
				cc[t7-1] = ch[t8-1]+ch[t9-1];
				cc[t6-1] = ch[t8-1]-ch[t9-1];
				cc[t7] = ch[t8]+ch[t9];
				cc[t6] = ch[t9]-ch[t8];
				t6 += t10;
				t7 += t10;
				t8 += ido;
				t9 += ido;
			}
		}
	}
}

/*
 * dradfg's arithmetic for ip = 5 in one sweep: every element gets the products and sums
 * dradfg computes for it, in the same order, with the same rotation factors (rot), reading
 * cc(ido, l1, 5) and writing ch(ido, 5, l1). ido must be odd.
 */
STIN void dradf5(int ido, int l1, double *cc, double *ch, double *wa, double *rot)
{
	int		i, ic, k, t0;
	double	*x0, *x1, *x2, *x3, *x4, *y, *wa1, *wa2, *wa3, *wa4;
	double	ar11, ai11, ar21, ai21, ar12, ai12, ar22, ai22;
	double	tr1, ti1, tr2, ti2, tr3, ti3, tr4, ti4;
	double	ar1, ar2, ar3, ar4, ai1, ai2, ai3, ai4;
	double	br0, br1, br2, br3, br4, bi0, bi1, bi2, bi3, bi4;

	ar11 = rot[0];
	ai11 = rot[1];
	ar21 = rot[2];
	ai21 = rot[3];
	ar12 = rot[4];
	ai12 = rot[5];
	ar22 = rot[6];
	ai22 = rot[7];
	wa1 = wa;
	wa2 = wa1+ido;
	wa3 = wa2+ido;
	wa4 = wa3+ido;
	t0 = l1*ido;

	for ( k=0; k<l1; k++ ) {
		x0 = cc+k*ido;
		x1 = x0+t0;
		x2 = x1+t0;
		x3 = x2+t0;
		x4 = x3+t0;
		y = ch+5*k*ido;

		ar1 = x1[0]+x4[0];
		ar4 = x4[0]-x1[0];
		ar2 = x2[0]+x3[0];
		ar3 = x3[0]-x2[0];
		br1 = x0[0]+ar11*ar1;
		br1 += ar21*ar2;
		br4 = ai11*ar4;
		br4 += ai21*ar3;
		br2 = x0[0]+ar12*ar1;
		br2 += ar22*ar2;
		br3 = ai12*ar4;
		br3 += ai22*ar3;
		br0 = x0[0]+ar1;
		br0 += ar2;
		y[0] = br0;
		y[2*ido-1] = br1;
		y[2*ido] = br4;
		y[4*ido-1] = br2;
		y[4*ido] = br3;

		for ( i=2; i<ido; i+=2 ) {
			ic = ido-i;
			tr1 = wa1[i-2]*x1[i-1]+wa1[i-1]*x1[i];
			ti1 = wa1[i-2]*x1[i]-wa1[i-1]*x1[i-1];
			tr2 = wa2[i-2]*x2[i-1]+wa2[i-1]*x2[i];
			ti2 = wa2[i-2]*x2[i]-wa2[i-1]*x2[i-1];
			tr3 = wa3[i-2]*x3[i-1]+wa3[i-1]*x3[i];
			ti3 = wa3[i-2]*x3[i]-wa3[i-1]*x3[i-1];
			tr4 = wa4[i-2]*x4[i-1]+wa4[i-1]*x4[i];
			ti4 = wa4[i-2]*x4[i]-wa4[i-1]*x4[i-1];

			ar1 = tr1+tr4;
			ar4 = ti1-ti4;
			ai1 = ti1+ti4;
			ai4 = tr4-tr1;
			ar2 = tr2+tr3;
			ar3 = ti2-ti3;
			ai2 = ti2+ti3;
			ai3 = tr3-tr2;

			br1 = x0[i-1]+ar11*ar1;
			br1 += ar21*ar2;
			br4 = ai11*ar4;
			br4 += ai21*ar3;
			br2 = x0[i-1]+ar12*ar1;
			br2 += ar22*ar2;
			br3 = ai12*ar4;
			br3 += ai22*ar3;
			br0 = x0[i-1]+ar1;
			br0 += ar2;
			bi1 = x0[i]+ar11*ai1;
			bi1 += ar21*ai2;
			bi4 = ai11*ai4;
			bi4 += ai21*ai3;
			bi2 = x0[i]+ar12*ai1;
			bi2 += ar22*ai2;
			bi3 = ai12*ai4;
			bi3 += ai22*ai3;
			bi0 = x0[i]+ai1;
			bi0 += ai2;

			y[i-1] = br0;
			y[i] = bi0;
			y[i-1+2*ido] = br1+br4;
			y[ic-1+ido] = br1-br4;
			y[i+2*ido] = bi1+bi4;
			y[ic+ido] = bi4-bi1;
			y[i-1+4*ido] = br2+br3;
			y[ic-1+3*ido] = br2-br3;
			y[i+4*ido] = bi2+bi3;
			y[ic+3*ido] = bi3-bi2;
		}
	}
}

/* The rotation factors dradfg derives for ip = 5, by its own recurrence */
STIN void radix5Rotations(double *rot)
{
	static double	tpi = 6.28318530717958647692528676655900577;
	double			arg, dcp, dsp, ar1, ai1, ar1h, ar2, ai2, ar2h, dc2, ds2;
	int				l;

	arg = tpi/(double)5;
	dcp = cos(arg);
	dsp = sin(arg);
	ar1 = 1.0;
	ai1 = 0.0;
	for ( l=1; l<3; l++ ) {
		ar1h = dcp*ar1-dsp*ai1;
		ai1 = dcp*ai1+dsp*ar1;
		ar1 = ar1h;
		dc2 = ar1;
		ds2 = ai1;
		ar2 = ar1;
		ai2 = ai1;
		ar2h = dc2*ar2-ds2*ai2;
		ai2 = dc2*ai2+ds2*ar2;
		ar2 = ar2h;
		rot[4*l-4] = ar1;
		rot[4*l-3] = ai1;
		rot[4*l-2] = ar2;
		rot[4*l-1] = ai2;
	}
}

/*
 * Two groups of a pass at once. When l1 is even, the data seen as pairs of doubles is the
 * same pass with l1/2 groups: lane 0 holds groups 0, 2, 4, ... and lane 1 groups 1, 3, 5, ...
 * Every lane does the scalar pass's operations in the scalar pass's order.
 */
typedef double pair __attribute__ ((vector_size (16), may_alias, aligned (8)));

/* dradf4 on pairs, ido odd or even */
STIN void dradf4p(int ido, int l1, pair *cc, pair *ch, double *wa1, double *wa2, double *wa3)
{
	static double	hsqt2 = .70710678118654752440084436210485;
	int				i, ic, k, t0;
	pair			*x0, *x1, *x2, *x3, *y;
	pair			ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;

	t0 = l1*ido;
	for ( k=0; k<l1; k++ ) {
		x0 = cc+k*ido;
		x1 = x0+t0;
		x2 = x1+t0;
		x3 = x2+t0;
		y = ch+4*k*ido;

		tr1 = x1[0]+x3[0];
		tr2 = x0[0]+x2[0];
		y[0] = tr1+tr2;
		y[4*ido-1] = tr2-tr1;
		y[2*ido-1] = x0[0]-x2[0];
		y[2*ido] = x3[0]-x1[0];

		for ( i=2; i<ido; i+=2 ) {
			ic = ido-i;
			cr2 = wa1[i-2]*x1[i-1]+wa1[i-1]*x1[i];
			ci2 = wa1[i-2]*x1[i]-wa1[i-1]*x1[i-1];
			cr3 = wa2[i-2]*x2[i-1]+wa2[i-1]*x2[i];
			ci3 = wa2[i-2]*x2[i]-wa2[i-1]*x2[i-1];
			cr4 = wa3[i-2]*x3[i-1]+wa3[i-1]*x3[i];
			ci4 = wa3[i-2]*x3[i]-wa3[i-1]*x3[i-1];

			tr1 = cr2+cr4;
			tr4 = cr4-cr2;
			ti1 = ci2+ci4;
			ti4 = ci2-ci4;
			ti2 = x0[i]+ci3;
			ti3 = x0[i]-ci3;
			tr2 = x0[i-1]+cr3;
			tr3 = x0[i-1]-cr3;

			y[i-1] = tr1+tr2;
			y[i] = ti1+ti2;
			y[ic-1+ido] = tr3-ti4;
			y[ic+ido] = tr4-ti3;
			y[i-1+2*ido] = ti4+tr3;
			y[i+2*ido] = tr4+ti3;
			y[ic-1+3*ido] = tr2-tr1;
			y[ic+3*ido] = ti1-ti2;
		}

		if ( ido%2 == 0 ) {
			ti1 = -hsqt2*(x1[ido-1]+x3[ido-1]);
			tr1 = hsqt2*(x1[ido-1]-x3[ido-1]);
			y[ido-1] = tr1+x0[ido-1];
			y[3*ido-1] = x0[ido-1]-tr1;
			y[ido] = ti1-x2[ido-1];
			y[3*ido] = ti1+x2[ido-1];
		}
	}
}

/* dradf5 on pairs; ido must be odd, as it is in the paired passes, where only 5s follow a 5
   in ifac */
STIN void dradf5p(int ido, int l1, pair *cc, pair *ch, double *wa, double *rot)
{
	int		i, ic, k, t0;
	pair	*x0, *x1, *x2, *x3, *x4, *y;
	double	*wa1, *wa2, *wa3, *wa4;
	double	ar11, ai11, ar21, ai21, ar12, ai12, ar22, ai22;
	pair	tr1, ti1, tr2, ti2, tr3, ti3, tr4, ti4;
	pair	ar1, ar2, ar3, ar4, ai1, ai2, ai3, ai4;
	pair	br0, br1, br2, br3, br4, bi0, bi1, bi2, bi3, bi4;

	ar11 = rot[0];
	ai11 = rot[1];
	ar21 = rot[2];
	ai21 = rot[3];
	ar12 = rot[4];
	ai12 = rot[5];
	ar22 = rot[6];
	ai22 = rot[7];
	wa1 = wa;
	wa2 = wa1+ido;
	wa3 = wa2+ido;
	wa4 = wa3+ido;
	t0 = l1*ido;

	for ( k=0; k<l1; k++ ) {
		x0 = cc+k*ido;
		x1 = x0+t0;
		x2 = x1+t0;
		x3 = x2+t0;
		x4 = x3+t0;
		y = ch+5*k*ido;

		ar1 = x1[0]+x4[0];
		ar4 = x4[0]-x1[0];
		ar2 = x2[0]+x3[0];
		ar3 = x3[0]-x2[0];
		br1 = x0[0]+ar11*ar1;
		br1 += ar21*ar2;
		br4 = ai11*ar4;
		br4 += ai21*ar3;
		br2 = x0[0]+ar12*ar1;
		br2 += ar22*ar2;
		br3 = ai12*ar4;
		br3 += ai22*ar3;
		br0 = x0[0]+ar1;
		br0 += ar2;
		y[0] = br0;
		y[2*ido-1] = br1;
		y[2*ido] = br4;
		y[4*ido-1] = br2;
		y[4*ido] = br3;

		for ( i=2; i<ido; i+=2 ) {
			ic = ido-i;
			tr1 = wa1[i-2]*x1[i-1]+wa1[i-1]*x1[i];
			ti1 = wa1[i-2]*x1[i]-wa1[i-1]*x1[i-1];
			tr2 = wa2[i-2]*x2[i-1]+wa2[i-1]*x2[i];
			ti2 = wa2[i-2]*x2[i]-wa2[i-1]*x2[i-1];
			tr3 = wa3[i-2]*x3[i-1]+wa3[i-1]*x3[i];
			ti3 = wa3[i-2]*x3[i]-wa3[i-1]*x3[i-1];
			tr4 = wa4[i-2]*x4[i-1]+wa4[i-1]*x4[i];
			ti4 = wa4[i-2]*x4[i]-wa4[i-1]*x4[i-1];

			ar1 = tr1+tr4;
			ar4 = ti1-ti4;
			ai1 = ti1+ti4;
			ai4 = tr4-tr1;
			ar2 = tr2+tr3;
			ar3 = ti2-ti3;
			ai2 = ti2+ti3;
			ai3 = tr3-tr2;

			br1 = x0[i-1]+ar11*ar1;
			br1 += ar21*ar2;
			br4 = ai11*ar4;
			br4 += ai21*ar3;
			br2 = x0[i-1]+ar12*ar1;
			br2 += ar22*ar2;
			br3 = ai12*ar4;
			br3 += ai22*ar3;
			br0 = x0[i-1]+ar1;
			br0 += ar2;
			bi1 = x0[i]+ar11*ai1;
			bi1 += ar21*ai2;
			bi4 = ai11*ai4;
			bi4 += ai21*ai3;
			bi2 = x0[i]+ar12*ai1;
			bi2 += ar22*ai2;
			bi3 = ai12*ai4;
			bi3 += ai22*ai3;
			bi0 = x0[i]+ai1;
			bi0 += ai2;

			y[i-1] = br0;
			y[i] = bi0;
			y[i-1+2*ido] = br1+br4;
			y[ic-1+ido] = br1-br4;
			y[i+2*ido] = bi1+bi4;
			y[ic+ido] = bi4-bi1;
			y[i-1+4*ido] = br2+br3;
			y[ic-1+3*ido] = br2-br3;
			y[i+4*ido] = bi2+bi3;
			y[ic+3*ido] = bi3-bi2;
		}
	}
}

/* dradf4 with l1 = 1 on data left in pairs by the passes before: row j, element i at
   cc[2*(ido*(j/2)+i) + j%2] */
STIN void dradf4last(int ido, double *cc, double *ch, double *wa1, double *wa2, double *wa3)
{
	static double	hsqt2 = .70710678118654752440084436210485;
	int				i, ic;
	double			*x0, *x1, *x2, *x3;
	double			ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;

	x0 = cc;
	x1 = cc+1;
	x2 = cc+2*ido;
	x3 = cc+2*ido+1;

	tr1 = x1[0]+x3[0];
	tr2 = x0[0]+x2[0];
	ch[0] = tr1+tr2;
	ch[4*ido-1] = tr2-tr1;
	ch[2*ido-1] = x0[0]-x2[0];
	ch[2*ido] = x3[0]-x1[0];

	for ( i=2; i<ido; i+=2 ) {
		ic = ido-i;
		cr2 = wa1[i-2]*x1[2*i-2]+wa1[i-1]*x1[2*i];
		ci2 = wa1[i-2]*x1[2*i]-wa1[i-1]*x1[2*i-2];
		cr3 = wa2[i-2]*x2[2*i-2]+wa2[i-1]*x2[2*i];
		ci3 = wa2[i-2]*x2[2*i]-wa2[i-1]*x2[2*i-2];
		cr4 = wa3[i-2]*x3[2*i-2]+wa3[i-1]*x3[2*i];
		ci4 = wa3[i-2]*x3[2*i]-wa3[i-1]*x3[2*i-2];

		tr1 = cr2+cr4;
		tr4 = cr4-cr2;
		ti1 = ci2+ci4;
		ti4 = ci2-ci4;
		ti2 = x0[2*i]+ci3;
		ti3 = x0[2*i]-ci3;
		tr2 = x0[2*i-2]+cr3;
		tr3 = x0[2*i-2]-cr3;

		ch[i-1] = tr1+tr2;
		ch[i] = ti1+ti2;
		ch[ic-1+ido] = tr3-ti4;
		ch[ic+ido] = tr4-ti3;
		ch[i-1+2*ido] = ti4+tr3;
		ch[i+2*ido] = tr4+ti3;
		ch[ic-1+3*ido] = tr2-tr1;
		ch[ic+3*ido] = ti1-ti2;
	}

	if ( ido%2 == 0 ) {
		ti1 = -hsqt2*(x1[2*ido-2]+x3[2*ido-2]);
		tr1 = hsqt2*(x1[2*ido-2]-x3[2*ido-2]);
		ch[ido-1] = tr1+x0[2*ido-2];
		ch[3*ido-1] = x0[2*ido-2]-tr1;
		ch[ido] = ti1-x2[2*ido-2];
		ch[3*ido] = ti1+x2[2*ido-2];
	}
}

/* All passes but the last on pairs: needs l1 even in each of them (n even, factors 4 and 5
   only) and a last factor of 4 */
STIN int pairedPasses(int *ifac)
{
	int		k1, nf = ifac[1];

	if ( (nf < 2) || (ifac[2] != 4) )
		return 0;
	for ( k1=1; k1<nf; k1++ )
		if ( (ifac[k1+2] != 4) && (ifac[k1+2] != 5) )
			return 0;
	return 1;
}

/* The first of the paired passes (ido = 1), its +1/-1 inputs read from the bits: the pair in
   row j of group q holds bits 2q+l1*j and 2q+1+l1*j */
STIN void firstPassBits(int ip, int l1, const uint64_t *bits, pair *ch, double *rot)
{
	static const pair	pm[4] = { { -1.0, -1.0 }, { -1.0, 1.0 }, { 1.0, -1.0 }, { 1.0, 1.0 } };
	long	q, p;
	int		j;
	pair	x[5], *y, ar1, ar2, ar3, ar4, br0, br1, br2, br3, br4, tr1, tr2;

	for ( q=0; q<l1/2; q++ ) {
		for ( j=0; j<ip; j++ ) {
			p = 2*q+(long)l1*j;
			x[j] = pm[(bits[p >> 6] >> (62 - (p & 63))) & 3];
		}
		y = ch+ip*q;
		if ( ip == 4 ) {
			tr1 = x[1]+x[3];
			tr2 = x[0]+x[2];
			y[0] = tr1+tr2;
			y[3] = tr2-tr1;
			y[1] = x[0]-x[2];
			y[2] = x[3]-x[1];
		}
		else {
			ar1 = x[1]+x[4];
			ar4 = x[4]-x[1];
			ar2 = x[2]+x[3];
			ar3 = x[3]-x[2];
			br1 = x[0]+rot[0]*ar1;
			br1 += rot[2]*ar2;
			br4 = rot[1]*ar4;
			br4 += rot[3]*ar3;
			br2 = x[0]+rot[4]*ar1;
			br2 += rot[6]*ar2;
			br3 = rot[5]*ar4;
			br3 += rot[7]*ar3;
			br0 = x[0]+ar1;
			br0 += ar2;
			y[0] = br0;
			y[1] = br1;
			y[2] = br4;
			y[3] = br2;
			y[4] = br3;
		}
	}
}

/*
 * The forward transform of the n values 2*bit-1 of the packed bits (bit i is bit 63-i%64 of
 * word i/64), with the arithmetic of __ogg_fdrfftf. c and wsave are the work arrays of
 * __ogg_fdrfftf; the result is left in c or in wsave[0..n-1], and that one is returned.
 * in holds the data before each pass; a pass writes out, except dradfg with ido > 1, which
 * works in place with out as scratch.
 */
double *__ogg_fdrfftfbits(int n, const uint64_t *bits, double *c, double *wsave, int *ifac)
{
	int		i, k1, l1, l2;
	int		nf, ip, iw, ido, idl1, paired;
	double	rot[8], *wa, *in, *out, *t;
	
	nf = ifac[1];
	l2 = n;
	iw = n;
	in = c;
	out = wsave;
	wa = wsave+n;
	paired = (n > 1) && pairedPasses(ifac);
	if ( !paired )
		for ( i=0; i<n; i++ )
			c[i] = 2*(int)((bits[i >> 6] >> (63 - (i & 63))) & 1) - 1;
	if ( n == 1 )
		return c;
	radix5Rotations(rot);
	
	for ( k1=0; k1<nf; k1++ ) {
		ip = ifac[nf-k1+1];
		l1 = l2/ip;
		ido = n/l2;
		idl1 = ido*l1;
		iw -= (ip-1)*ido;
		l2 = l1;
		
		if ( paired && (k1 == 0) )
			firstPassBits(ip, l1, bits, (pair *)out, rot);
		else if ( paired && (k1 < nf-1) && (ip == 4) )
			dradf4p(ido, l1/2, (pair *)in, (pair *)out, wa+iw-1, wa+iw+ido-1, wa+iw+2*ido-1);
		else if ( paired && (k1 < nf-1) )
			dradf5p(ido, l1/2, (pair *)in, (pair *)out, wa+iw-1, rot);
		else if ( paired )
			dradf4last(ido, in, out, wa+iw-1, wa+iw+ido-1, wa+iw+2*ido-1);
		else if ( ip == 4 )
			dradf4(ido, l1, in, out, wa+iw-1, wa+iw+ido-1, wa+iw+2*ido-1);
		else if ( ip == 2 )
			dradf2(ido, l1, in, out, wa+iw-1);
		else if ( (ip == 5) && (ido%2 == 1) )
			dradf5(ido, l1, in, out, wa+iw-1, rot);
		else if ( ido == 1 )
			dradfg(ido, ip, l1, idl1, out, out, out, in, in, wa+iw-1);
		else {
			dradfg(ido, ip, l1, idl1, in, in, in, out, out, wa+iw-1);
			continue;
		}
		t = in;
		in = out;
		out = t;
	}
	return in;
}
