/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
U T I L I T I E S
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/generators.h"
#include "../include/stat_fncs.h"

#define ASCII_BUFFER	(1 << 20)

static long				wordsAllocated;
static int				onesCount;
static unsigned long	sequenceId;

static unsigned char	*asciiBuffer;
static size_t			asciiStart, asciiEnd;
static int				asciiAtEnd;

static void	runTests(void);

int
displayGeneratorOptions()
{
	int		option = 0;

	printf("           G E N E R A T O R    S E L E C T I O N \n");
	printf("           ______________________________________\n\n");
	printf("    [0] Input File                 [1] Linear Congruential\n");
	printf("    [2] Quadratic Congruential I   [3] Quadratic Congruential II\n");
	printf("    [4] Cubic Congruential         [5] XOR\n");
	printf("    [6] Modular Exponentiation     [7] Blum-Blum-Shub\n");
	printf("    [8] Micali-Schnorr             [9] G Using SHA-1\n\n");
	printf("   Enter Choice: ");
	scanf("%d", &option);
	printf("\n\n");

	return option;
}


int
generatorOptions(char** streamFile)
{
	char	file[200];
	int		option = NUMOFGENERATORS+1;
	FILE	*fp;
	
	while ( (option < 0) || (option > NUMOFGENERATORS) ) {
		option = displayGeneratorOptions();
		switch( option ) {
			case 0:
				printf("\t\tUser Prescribed Input File: ");
				scanf("%s", file);
				*streamFile = (char*)calloc(200, sizeof(char));
				sprintf(*streamFile, "%s", file);
				printf("\n");
				if ( (fp = fopen(*streamFile, "r")) == NULL ) {
					printf("File Error:  file %s could not be opened.\n",  *streamFile);
					exit(-1);
				}
				else
					fclose(fp);
				break;
			case 1:
				*streamFile = "Linear-Congruential";
				break;
			case 2:
				*streamFile = "Quadratic-Congruential-1";
				break;
			case 3:
				*streamFile = "Quadratic-Congruential-2";
				break;
			case 4:
				*streamFile = "Cubic-Congruential";
				break;
			case 5:
				*streamFile = "XOR";
				break;
			case 6:
				*streamFile = "Modular-Exponentiation";
				break;
			case 7:
				*streamFile = "Blum-Blum-Shub";
				break;
			case 8:
				*streamFile = "Micali-Schnorr";
				break;
			case 9:
				*streamFile = "G using SHA-1";
				break;
				
			/* INTRODUCE NEW PRNG NAMES HERE */
			/*
			case 10:  *streamFile = "myNewPRNG";
				break;
			*/
			default:
				printf("Error:  Out of range - Try again!\n");
				break;
		}
	}
	return option;
}


void
chooseTests()
{
	int		i;
	
	printf("                S T A T I S T I C A L   T E S T S\n");
	printf("                _________________________________\n\n");
	printf("    [01] Frequency                       [02] Block Frequency\n");
	printf("    [03] Cumulative Sums                 [04] Runs\n");
	printf("    [05] Longest Run of Ones             [06] Rank\n");
	printf("    [07] Discrete Fourier Transform      [08] Nonperiodic Template Matchings\n");
	printf("    [09] Overlapping Template Matchings  [10] Universal Statistical\n");
	printf("    [11] Approximate Entropy             [12] Random Excursions\n");
	printf("    [13] Random Excursions Variant       [14] Serial\n");
	printf("    [15] Linear Complexity\n\n");
	printf("         INSTRUCTIONS\n");
	printf("            Enter 0 if you DO NOT want to apply all of the\n");
	printf("            statistical tests to each sequence and 1 if you DO.\n\n");
	printf("   Enter Choice: ");
	scanf("%d", &testVector[0]);
	printf("\n");
	if ( testVector[0] == 1 )
		for( i=1; i<=NUMOFTESTS; i++ )
			testVector[i] = 1;
	else {
		printf("         INSTRUCTIONS\n");
		printf("            Enter a 0 or 1 to indicate whether or not the numbered statistical\n");
		printf("            test should be applied to each sequence.\n\n");
		printf("      123456789111111\n");
		printf("               012345\n");
		printf("      ");
		for ( i=1; i<=NUMOFTESTS; i++ ) 
			scanf("%1d", &testVector[i]);
		printf("\n\n");
	}
}


void
fixParameters()
{
	int		counter, testid;
	
	//  Check to see if any parameterized tests are selected
	if ( (testVector[TEST_BLOCK_FREQUENCY] != 1) && (testVector[TEST_NONPERIODIC] != 1) && 
		 (testVector[TEST_OVERLAPPING] != 1) && (testVector[TEST_APEN] != 1) &&
		 (testVector[TEST_SERIAL] != 1) && (testVector[TEST_LINEARCOMPLEXITY] != 1) )
			return;
		
	do {
		counter = 1;
		printf("        P a r a m e t e r   A d j u s t m e n t s\n");
		printf("        -----------------------------------------\n");
		if ( testVector[TEST_BLOCK_FREQUENCY] == 1 )
			printf("    [%d] Block Frequency Test - block length(M):         %d\n", counter++, tp.blockFrequencyBlockLength);
		if ( testVector[TEST_NONPERIODIC] == 1 )
			printf("    [%d] NonOverlapping Template Test - block length(m): %d\n", counter++, tp.nonOverlappingTemplateBlockLength);
		if ( testVector[TEST_OVERLAPPING] == 1 )
			printf("    [%d] Overlapping Template Test - block length(m):    %d\n", counter++, tp.overlappingTemplateBlockLength);
		if ( testVector[TEST_APEN] == 1 )
			printf("    [%d] Approximate Entropy Test - block length(m):     %d\n", counter++, tp.approximateEntropyBlockLength);
		if ( testVector[TEST_SERIAL] == 1 )
			printf("    [%d] Serial Test - block length(m):                  %d\n", counter++, tp.serialBlockLength);
		if ( testVector[TEST_LINEARCOMPLEXITY] == 1 )
			printf("    [%d] Linear Complexity Test - block length(M):       %d\n", counter++, tp.linearComplexitySequenceLength);
		printf("\n");
		printf("   Select Test (0 to continue): ");
		scanf("%1d", &testid);
		printf("\n");
		
		counter = 0;
		if ( testVector[TEST_BLOCK_FREQUENCY] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Block Frequency Test block length: ");
				scanf("%d", &tp.blockFrequencyBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_NONPERIODIC] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter NonOverlapping Template Test block Length: ");
				scanf("%d", &tp.nonOverlappingTemplateBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_OVERLAPPING] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Overlapping Template Test block Length: ");
				scanf("%d", &tp.overlappingTemplateBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_APEN] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Approximate Entropy Test block Length: ");
				scanf("%d", &tp.approximateEntropyBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_SERIAL] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Serial Test block Length: ");
				scanf("%d", &tp.serialBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_LINEARCOMPLEXITY] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Linear Complexity Test block Length: ");
				scanf("%d", &tp.linearComplexitySequenceLength);
				printf("\n");
				continue;
			}
		}
	} while ( testid != 0 );
}


void
fileBasedBitStreams(char *streamFile)
{
	FILE	*fp;
	int		mode;
	
	printf("   Input File Format:\n");
	printf("    [0] ASCII - A sequence of ASCII 0's and 1's\n");
	printf("    [1] Binary - Each byte in data file contains 8 bits of data\n\n");
	printf("   Select input mode:  ");
	scanf("%1d", &mode);
	printf("\n");
	if ( mode == 0 ) {
		if ( (fp = fopen(streamFile, "r")) == NULL ) {
			printf("ERROR IN FUNCTION fileBasedBitStreams:  file %s could not be opened.\n",  streamFile);
			exit(-1);
		}
		readBinaryDigitsInASCIIFormat(fp, streamFile);
		fclose(fp);
	}
	else if ( mode == 1 ) {
		if ( (fp = fopen(streamFile, "rb")) == NULL ) {
			printf("ERROR IN FUNCTION fileBasedBitStreams:  file %s could not be opened.\n", streamFile);
			exit(-1);
		}
		readHexDigitsInBinaryFormat(fp);
		fclose(fp);
	}
}


static int
allocateSequence(void)
{
	long	words = ((long)tp.n + 63) / 64 + 2;

	if ( words > wordsAllocated ) {
		free(epsilonWords);
		if ( (epsilonWords = (uint64_t *) calloc(words, sizeof(uint64_t))) == NULL ) {
			wordsAllocated = 0;
			return 0;
		}
		wordsAllocated = words;
	}
	return 1;
}


/* Called once the packed words of a new stream are in place */
static void
sequenceLoaded(void)
{
	long	full = tp.n / 64, words = ((long)tp.n + 63) / 64;

	if ( tp.n % 64 )
		epsilonWords[full] &= ~0ULL << (64 - tp.n % 64);
	epsilonWords[words] = 0;
	epsilonWords[words+1] = 0;
	onesCount = popcountWords(epsilonWords, words);
	sequenceId++;
}


static inline uint64_t
load64be(const unsigned char *p)
{
	return ((uint64_t)p[0] << 56) | ((uint64_t)p[1] << 48) | ((uint64_t)p[2] << 40) | ((uint64_t)p[3] << 32) |
		   ((uint64_t)p[4] << 24) | ((uint64_t)p[5] << 16) | ((uint64_t)p[6] << 8) | (uint64_t)p[7];
}


static inline uint64_t
load64le(const unsigned char *p)
{
	return (uint64_t)p[0] | ((uint64_t)p[1] << 8) | ((uint64_t)p[2] << 16) | ((uint64_t)p[3] << 24) |
		   ((uint64_t)p[4] << 32) | ((uint64_t)p[5] << 40) | ((uint64_t)p[6] << 48) | ((uint64_t)p[7] << 56);
}


static void
fillAsciiBuffer(FILE *fp)
{
	size_t	left = asciiEnd - asciiStart;

	memmove(asciiBuffer, asciiBuffer + asciiStart, left);
	asciiStart = 0;
	asciiEnd = left;
	asciiEnd += fread(asciiBuffer + asciiEnd, 1, ASCII_BUFFER - asciiEnd, fp);
	if ( asciiEnd < ASCII_BUFFER )
		asciiAtEnd = 1;
	/* past the end of the data the buffer reads as blanks */
	memset(asciiBuffer + asciiEnd, ' ', 8);
}


/* The top bit of each byte of the eight characters in x that is not a 0 or a 1 */
static inline uint64_t
nonDigits(uint64_t x)
{
	uint64_t	y = (x ^ 0x3030303030303030ULL) & 0xFEFEFEFEFEFEFEFEULL;

	return (((y & 0x7F7F7F7F7F7F7F7FULL) + 0x7F7F7F7F7F7F7F7FULL) | y) & 0x8080808080808080ULL;
}


/* The low bits of the eight characters in x as one byte, the first character on top */
static inline uint64_t
digitByte(uint64_t x)
{
	return ((x & 0x0101010101010101ULL) * 0x8040201008040201ULL) >> 56;
}


/* Reads the next n digits into epsilonWords, skipping whitespace; returns the count read */
static int
readAsciiBits(FILE *fp, int n)
{
	uint64_t			acc = 0, x, nondigit, b, bad, v;
	int					got = 0, accBits = 0, lead, take, r, j;
	long				k = 0;
	const unsigned char	*p;

	while ( got < n ) {
		if ( asciiStart + 64 > asciiEnd && !asciiAtEnd )
			fillAsciiBuffer(fp);
		if ( asciiStart >= asciiEnd )
			break;
		p = asciiBuffer + asciiStart;
		x = load64le(p);
		nondigit = nonDigits(x);
		if ( !nondigit && (n - got >= 64) && (asciiStart + 64 <= asciiEnd) && !nonDigits(load64le(p + 56)) ) {
			/* 64 digits in a row, as in files with long lines */
			v = digitByte(x);
			bad = 0;
			for ( j=8; j<64; j+=8 ) {
				x = load64le(p + j);
				bad |= nonDigits(x);
				v = (v << 8) | digitByte(x);
			}
			if ( !bad ) {
				epsilonWords[k++] = accBits ? (acc << (64 - accBits)) | (v >> accBits) : v;
				acc = v & ((1ULL << accBits) - 1);
				got += 64;
				asciiStart += 64;
				continue;
			}
			x = load64le(p);
		}
		lead = nondigit ? __builtin_ctzll(nondigit) >> 3 : 8;
		take = lead < n - got ? lead : n - got;
		if ( take > 0 ) {
			b = digitByte(x) >> (8 - take);
			if ( accBits + take < 64 ) {
				acc = (acc << take) | b;
				accBits += take;
			}
			else {
				r = 64 - accBits;
				epsilonWords[k++] = (acc << r) | (b >> (take - r));
				accBits = take - r;
				acc = b & ((1ULL << accBits) - 1);
			}
			got += take;
			asciiStart += take;
		}
		if ( take == lead && lead < 8 )
			asciiStart++;
	}
	if ( accBits > 0 )
		epsilonWords[k] = acc << (64 - accBits);
	return got;
}


void
readBinaryDigitsInASCIIFormat(FILE *fp, char *streamFile)
{
	int		i, bitsRead;
	
	if ( !allocateSequence() ||
		 ((asciiBuffer == NULL) && ((asciiBuffer = (unsigned char *) malloc(ASCII_BUFFER + 8)) == NULL)) ) {
		printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
		printf("Statistical Testing Aborted!\n");
		return;
	}
	asciiStart = asciiEnd = 0;
	asciiAtEnd = 0;
	printf("     Statistical Testing In Progress.........\n\n");   
	for ( i=0; i<tp.numOfBitStreams; i++ ) {
		bitsRead = readAsciiBits(fp, tp.n);
		if ( bitsRead < tp.n ) {
			printf("ERROR:  Insufficient data in file %s.  %d bits were read.\n", streamFile, bitsRead);
			break;
		}
		sequenceLoaded();
		fprintf(freqfp, "\t\tBITSREAD = %d 0s = %d 1s = %d\n", tp.n, tp.n - onesCount, onesCount);
		runTests();
	}
}


void
readHexDigitsInBinaryFormat(FILE *fp)
{
	int				i;
	long			k, full = tp.n / 64;
	size_t			streamBytes = ((size_t)tp.n + 31) / 32 * 4;
	unsigned char	*buffer = NULL;
	
	if ( !allocateSequence() ||
		 ((buffer = (unsigned char *) calloc(streamBytes + 8, 1)) == NULL) ) {
		printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
		return;
	}

	printf("     Statistical Testing In Progress.........\n\n");   
	for ( i=0; i<tp.numOfBitStreams; i++ ) {
		if ( fread(buffer, 1, streamBytes, fp) != streamBytes ) {
			printf("READ ERROR:  Insufficient data in file.\n");
			break;
		}
		for ( k=0; k<full; k++ )
			epsilonWords[k] = load64be(buffer + 8*k);
		if ( tp.n % 64 )
			epsilonWords[full] = load64be(buffer + 8*full);
		sequenceLoaded();
		fprintf(freqfp, "\t\tBITSREAD = %d 0s = %d 1s = %d\n", tp.n, tp.n - onesCount, onesCount);
		runTests();
	}
	free(buffer);
}


int
convertToBits(BYTE *x, int xBitLength, int bitsNeeded, int *num_0s, int *num_1s, int *bitsRead)
{
	int		i, j, count, bit;
	BYTE	mask;
	int		zeros, ones;

	count = 0;
	zeros = ones = 0;
	for ( i=0; i<(xBitLength+7)/8; i++ ) {
		mask = 0x80;
		for ( j=0; j<8; j++ ) {
			if ( *(x+i) & mask ) {
				bit = 1;
				(*num_1s)++;
				ones++;
			}
			else {
				bit = 0;
				(*num_0s)++;
				zeros++;
			}
			mask >>= 1;
			epsilon[*bitsRead] = bit;
			(*bitsRead)++;
			if ( *bitsRead == bitsNeeded )
				return 1;
			if ( ++count == xBitLength )
				return 0;
		}
	}
	
	return 0;
}


void
openOutputStreams(int option)
{
	int		i, numOfBitStreams, numOfOpenFiles = 0;
	char	freqfn[200], summaryfn[200], statsDir[200], resultsDir[200];
	
	sprintf(freqfn, "experiments/%s/freq.txt", generatorDir[option]);
	if ( (freqfp = fopen(freqfn, "w")) == NULL ) {
		printf("\t\tMAIN:  Could not open freq file: <%s>", freqfn);
		exit(-1);
	}
	sprintf(summaryfn, "experiments/%s/finalAnalysisReport.txt", generatorDir[option]);
	if ( (summary = fopen(summaryfn, "w")) == NULL ) {
		printf("\t\tMAIN:  Could not open stats file: <%s>", summaryfn);
		exit(-1);
	}
	
	for( i=1; i<=NUMOFTESTS; i++ ) {
		if ( testVector[i] == 1 ) {
			sprintf(statsDir, "experiments/%s/%s/stats.txt", generatorDir[option], testNames[i]);
			sprintf(resultsDir, "experiments/%s/%s/results.txt", generatorDir[option], testNames[i]);
			if ( (stats[i] = fopen(statsDir, "w")) == NULL ) {	/* STATISTICS LOG */
				printf("ERROR: LOG FILES COULD NOT BE OPENED.\n");
				printf("       MAX # OF OPENED FILES HAS BEEN REACHED = %d\n", numOfOpenFiles);
				printf("-OR-   THE OUTPUT DIRECTORY DOES NOT EXIST.\n");
				exit(-1);
			}
			else
				numOfOpenFiles++;
			if ( (results[i] = fopen(resultsDir, "w")) == NULL ) {	/* P_VALUES LOG   */
				 printf("ERROR: LOG FILES COULD NOT BE OPENED.\n");
				 printf("       MAX # OF OPENED FILES HAS BEEN REACHED = %d\n", numOfOpenFiles);
				 printf("-OR-   THE OUTPUT DIRECTORY DOES NOT EXIST.\n");
				 exit(-1);
			}
			else
				numOfOpenFiles++;
		}
	}
	printf("   How many bitstreams? ");
	scanf("%d", &numOfBitStreams);
	tp.numOfBitStreams = numOfBitStreams;
	printf("\n");
}


void
invokeTestSuite(int option, char *streamFile)
{
	fprintf(freqfp, "________________________________________________________________________________\n\n");
	fprintf(freqfp, "\t\tFILE = %s\t\tALPHA = %6.4f\n", streamFile, ALPHA);
	fprintf(freqfp, "________________________________________________________________________________\n\n");
	if ( option != 0 )
		printf("     Statistical Testing In Progress.........\n\n");
	switch( option ) {
		case 0:
			fileBasedBitStreams(streamFile);
			break;
		case 1:
			lcg();
			break;
		case 2:
			quadRes1();
			break;
		case 3:
			quadRes2();
			break;
		case 4:
			cubicRes();
			break;
		case 5:
			exclusiveOR();
			break;
		case 6:
			modExp();
			break;
		case 7:
			bbs();
			break;
		case 8:
			micali_schnorr();
			break;
		case 9:
			SHA1();
			break;
			
		/* INTRODUCE NEW PSEUDO RANDOM NUMBER GENERATORS HERE */
			
		default:
			printf("Error in invokeTestSuite!\n");
			break;
	}
	printf("     Statistical Testing Complete!!!!!!!!!!!!\n\n");
}


static void
runTests(void)
{
	if ( (testVector[0] == 1) || (testVector[TEST_FREQUENCY] == 1) ) 
		Frequency(tp.n);
	
	if ( (testVector[0] == 1) || (testVector[TEST_BLOCK_FREQUENCY] == 1) ) 
		BlockFrequency(tp.blockFrequencyBlockLength, tp.n);
	
	if ( (testVector[0] == 1) || (testVector[TEST_CUSUM] == 1) )
		CumulativeSums(tp.n);
	
	if ( (testVector[0] == 1) || (testVector[TEST_RUNS] == 1) )
		Runs(tp.n); 
	
	if ( (testVector[0] == 1) || (testVector[TEST_LONGEST_RUN] == 1) )
		LongestRunOfOnes(tp.n);
	
	if ( (testVector[0] == 1) || (testVector[TEST_RANK] == 1) )
		Rank(tp.n);
	
	if ( (testVector[0] == 1) || (testVector[TEST_FFT] == 1) )
		DiscreteFourierTransform(tp.n);
	
	if ( (testVector[0] == 1) || (testVector[TEST_NONPERIODIC] == 1) )
		NonOverlappingTemplateMatchings(tp.nonOverlappingTemplateBlockLength, tp.n);
	
	if ( (testVector[0] == 1) || (testVector[TEST_OVERLAPPING] == 1) )
		OverlappingTemplateMatchings(tp.overlappingTemplateBlockLength, tp.n);
	
	if ( (testVector[0] == 1) || (testVector[TEST_UNIVERSAL] == 1) )
		Universal(tp.n);
	
	if ( (testVector[0] == 1) || (testVector[TEST_APEN] == 1) )
		ApproximateEntropy(tp.approximateEntropyBlockLength, tp.n);
	
	if ( (testVector[0] == 1) || (testVector[TEST_RND_EXCURSION] == 1) )
		RandomExcursions(tp.n);
	
	if ( (testVector[0] == 1) || (testVector[TEST_RND_EXCURSION_VAR] == 1) )
		RandomExcursionsVariant(tp.n);
	
	if ( (testVector[0] == 1) || (testVector[TEST_SERIAL] == 1) )
		Serial(tp.serialBlockLength,tp.n);
	
	if ( (testVector[0] == 1) || (testVector[TEST_LINEARCOMPLEXITY] == 1) )
		LinearComplexity(tp.linearComplexitySequenceLength, tp.n);
}


/* Entry point for the built-in generators, which fill only epsilon */
void
nist_test_suite()
{
	long		i, k, words = ((long)tp.n + 63) / 64;
	uint64_t	w;

	if ( !allocateSequence() ) {
		printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
		return;
	}
	for ( k=0; k<words; k++ ) {
		w = 0;
		for ( i=64*k; (i<64*k+64) && (i<tp.n); i++ )
			if ( epsilon[i] )
				w |= 1ULL << (63 - (i - 64*k));
		epsilonWords[k] = w;
	}
	sequenceLoaded();
	runTests();
}


int
sequenceOnes(void)
{
	return onesCount;
}


int
popcountWords(const uint64_t *w, long count)
{
	long		i, chunk;
	uint64_t	x, acc;
	int			total = 0;

	while ( count > 0 ) {
		/* byte counters overflow after 31 words */
		chunk = count < 31 ? count : 31;
		acc = 0;
		for ( i=0; i<chunk; i++ ) {
			x = w[i];
			x = x - ((x >> 1) & 0x5555555555555555ULL);
			x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
			acc += (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
		}
		acc = (acc & 0x00FF00FF00FF00FFULL) + ((acc >> 8) & 0x00FF00FF00FF00FFULL);
		total += (int)((acc * 0x0001000100010001ULL) >> 48);
		w += chunk;
		count -= chunk;
	}
	return total;
}


/* Ones among bits start .. start+len-1 of a packed array */
int
rangeOnes(const uint64_t *w, long start, long len)
{
	long	first, last, end = start + len;

	if ( len <= 0 )
		return 0;
	first = start >> 6;
	last = (end - 1) >> 6;
	if ( first == last )
		return popcount64((w[first] << (start & 63)) >> (64 - len));
	return popcount64(w[first] << (start & 63)) + popcountWords(w + first + 1, last - first - 1) +
		   popcount64(w[last] >> ((64 - (end & 63)) & 63));
}


#define SPAN	25		/* visits to the states -SPAN..SPAN are counted */

/* A cycle ends: file its visits to -4..-1, 1..4 by count and start counting again */
static void
closeCycle(WalkSummary *ws, const int *visits, int *start)
{
	int		i, x, c;

	for ( i=0; i<8; i++ ) {
		x = i < 4 ? i-4 : i-3;
		c = visits[SPAN+x] - start[i];
		ws->cycles[c < 5 ? c : 5][i]++;
		start[i] = visits[SPAN+x];
	}
}


/* Steps taken one at a time, from the top bit of bits down */
static void
walkSteps(WalkSummary *ws, int *visits, int *start, int *S, int *sup, int *inf, uint64_t bits, int count)
{
	int		s = *S;

	for ( ; count>0; count--, bits<<=1 ) {
		s += (bits >> 63) ? 1 : -1;
		*sup = s > *sup ? s : *sup;
		*inf = s < *inf ? s : *inf;
		if ( (s >= -SPAN) && (s <= SPAN) ) {
			visits[SPAN+s]++;
			if ( s == 0 )
				closeCycle(ws, visits, start);
		}
	}
	*S = s;
}


/* Shared by CumulativeSums, RandomExcursions and RandomExcursionsVariant; computed once per stream */
const WalkSummary *
walkSummary(void)
{
	static WalkSummary		ws;
	static unsigned long	id = 0;
	static signed char		sum8[256], max8[256], min8[256];
	static unsigned char	hits8[256][17];		/* hits8[b][d]: steps of byte b that end d-8 above its start */
	static int				ready = 0;
	static int				*atWord = NULL;		/* S before each whole word */
	static long				atWordSize = 0;
	int						*grown;
	uint64_t	x;
	long		k, full = tp.n / 64;
	int			S = 0, sup = 0, inf = 0, visits[2*SPAN+1], start[8];
	int			i, j, b, d, t, hi, lo, ends;

	if ( id == sequenceId )
		return &ws;
	if ( full + 1 > atWordSize ) {
		if ( (grown = (int *)realloc(atWord, (full + 1) * sizeof(int))) == NULL ) {
			printf("Partial sums:  Insufficient memory available.\n");
			exit(-1);
		}
		atWord = grown;
		atWordSize = full + 1;
	}
	if ( !ready ) {
		memset(hits8, 0, sizeof(hits8));
		for ( b=0; b<256; b++ ) {
			t = 0;
			hi = -8;
			lo = 8;
			for ( j=7; j>=0; j-- ) {
				t += ((b >> j) & 1) ? 1 : -1;
				hits8[b][t+8]++;
				if ( t > hi )
					hi = t;
				if ( t < lo )
					lo = t;
			}
			sum8[b] = t;
			max8[b] = hi;
			min8[b] = lo;
		}
		ready = 1;
	}
	memset(&ws, 0, sizeof(ws));
	memset(visits, 0, sizeof(visits));
	memset(start, 0, sizeof(start));

	/* S at the word boundaries first; the extremes of those are bounds the words must beat */
	for ( k=0; k<full; k++ ) {
		atWord[k] = S;
		S += 2*popcount64(epsilonWords[k]) - 64;
		sup = S > sup ? S : sup;
		inf = S < inf ? S : inf;
	}
	atWord[full] = S;

	/*
	 * Over 64 steps from atWord[k] to atWord[k+1] the walk stays within half the sum of the
	 * two, plus or minus 32. Words that cannot pass the extremes found so far or come near
	 * zero are skipped. Near zero a byte adds its visits from a table, unless it reaches zero
	 * itself and so may end a cycle.
	 */
	for ( k=0; k<full; k++ ) {
		S = atWord[k];
		ends = S + atWord[k+1];
		if ( ((ends > 82) || (ends < -82)) && (ends + 64 <= 2*sup) && (ends - 64 >= 2*inf) )
			continue;
		x = epsilonWords[k];
		for ( j=56; j>=0; j-=8 ) {
			b = (int)(x >> j) & 0xFF;
			t = S + max8[b];
			sup = t > sup ? t : sup;
			t = S + min8[b];
			inf = t < inf ? t : inf;
			if ( (S >= -8) && (S <= 8) && hits8[b][8-S] )
				walkSteps(&ws, visits, start, &S, &sup, &inf, (uint64_t)b << 56, 8);
			else {
				if ( (S >= -17) && (S <= 17) )
					for ( d=0; d<17; d++ )
						visits[SPAN+S-8+d] += hits8[b][d];
				S += sum8[b];
			}
		}
	}
	S = atWord[full];
	if ( tp.n % 64 )
		walkSteps(&ws, visits, start, &S, &sup, &inf, epsilonWords[full], tp.n % 64);
	if ( S != 0 )
		closeCycle(&ws, visits, start);
	ws.S = S;
	ws.sup = sup;
	ws.inf = inf;
	ws.zeros = visits[SPAN];
	for ( i=-9; i<=9; i++ )
		ws.visits[i+9] = visits[SPAN+i];
	id = sequenceId;

	return &ws;
}
