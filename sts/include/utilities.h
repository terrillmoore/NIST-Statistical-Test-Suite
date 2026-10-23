/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
              U T I L I T Y  F U N C T I O N  P R O T O T Y P E S 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <stdint.h>

int		displayGeneratorOptions();
int		generatorOptions(char** streamFile);
void	chooseTests();
void	fixParameters();
void	fileBasedBitStreams(char *streamFile);
void	readBinaryDigitsInASCIIFormat(FILE *fp, char *streamFile);
void	readHexDigitsInBinaryFormat(FILE *fp);
int		convertToBits(BYTE *x, int xBitLength, int bitsNeeded, int *num_0s, int *num_1s, int *bitsRead);
void	openOutputStreams(int option);
void	invokeTestSuite(int option, char *streamFile);
void	nist_test_suite();

/* Partial sums S_k of the +1/-1 walk over the current stream, k = 1..n */
typedef struct _walkSummary {
	int		S;					/* S_n */
	int		sup, inf;			/* max and min of 0 and all S_k */
	int		zeros;				/* number of k with S_k = 0 */
	int		visits[19];			/* visits[x+9]: number of k with S_k = x */
	int		cycles[6][8];		/* cycles[j][i]: cycles (runs between zeros) with j visits */
								/* (5: five or more) to state i of -4..-1, 1..4 */
} WalkSummary;

int		sequenceOnes(void);
int		rangeOnes(const uint64_t *w, long start, long len);
int		popcountWords(const uint64_t *w, long count);
const WalkSummary	*walkSummary(void);
void	windowCounts(unsigned int *counts, long start, long count, int m);
void	cyclicWindowCounts(unsigned int *counts, int m, int n);

static inline int
popcount64(uint64_t x)
{
	x = x - ((x >> 1) & 0x5555555555555555ULL);
	x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
	x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
	return (int)((x * 0x0101010101010101ULL) >> 56);
}

/* The 64 bits of a packed stream that start at bit pos, the first one on top */
static inline uint64_t
bitWindow(const uint64_t *w, long pos)
{
	const uint64_t	*p = w + (pos >> 6);
	int				s = (int)(pos & 63);

	return (p[0] << s) | ((p[1] >> 1) >> (63 - s));
}

#endif
