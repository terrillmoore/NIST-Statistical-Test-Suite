#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/externs.h"
#include "../include/matrix.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
R A N K  A L G O R I T H M  R O U T I N E S
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * Ranks over GF(2) of two 32x32 matrices, each held in 16 words with two 32-bit rows per
 * word. The forward/backward elimination of the original byte-matrix code always ends with
 * as many nonzero rows as the rank: after the forward pass the matrix is upper triangular,
 * and every row the backward pass leaves without a pivot is already zero.
 *
 * Row k is used as column k of the transpose, which has the same rank; the bits of a
 * 32-bit half are the transpose's rows, so a pivot needs no row swap. The two matrices go
 * through the elimination together because each step waits on the one before it.
 */
void
computeRanks(const uint64_t *a, const uint64_t *b, int *rankA, int *rankB)
{
	uint64_t	x[16], y[16], mx, my, bx, by;
	uint32_t	usedX = 0, usedY = 0, cx, cy, rx, ry;
	int			col, u, sx, sy, kx = 0, ky = 0;

	memcpy(x, a, sizeof(x));
	memcpy(y, b, sizeof(y));
	for ( col=0; col<32; col++ ) {
		cx = (uint32_t)(x[col >> 1] >> ((~col & 1) << 5)) & ~usedX;
		cy = (uint32_t)(y[col >> 1] >> ((~col & 1) << 5)) & ~usedY;
		/* without a pivot the step below changes nothing */
		sx = cx ? __builtin_ctz(cx) : 0;
		sy = cy ? __builtin_ctz(cy) : 0;
		usedX |= cx & (0u - cx);
		usedY |= cy & (0u - cy);
		kx += cx != 0;
		ky += cy != 0;
		rx = cx & (cx - 1);
		ry = cy & (cy - 1);
		bx = ((uint64_t)rx << 32) | rx;
		by = ((uint64_t)ry << 32) | ry;
		/* columns before col are finished, so the word holding col-1 can change too */
		for ( u=col >> 1; u<16; u++ ) {
			mx = (x[u] >> sx) & 0x0000000100000001ULL;
			my = (y[u] >> sy) & 0x0000000100000001ULL;
			x[u] ^= bx & ((mx << 32) - mx);
			y[u] ^= by & ((my << 32) - my);
		}
	}
	*rankA = kx;
	*rankB = ky;
}
