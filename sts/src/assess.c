/* --------------------------------------------------------------------------
   Title       :  The NIST Statistical Test Suite

   Date        :  December 1999

   Programmer  :  Juan Soto

   Summary     :  For use in the evaluation of the randomness of bitstreams
                  produced by cryptographic random number generators.

   Package     :  Version 1.0

   Copyright   :  (c) 1999 by the National Institute Of Standards & Technology

   History     :  Version 1.0 by J. Soto, October 1999
                  Revised by J. Soto, November 1999
                  Revised by Larry Bassham, March 2008

   Keywords    :  Pseudorandom Number Generator (PRNG), Randomness, Statistical 
                  Tests, Complementary Error functions, Incomplete Gamma 
                  Function, Random Walks, Rank, Fast Fourier Transform, 
                  Template, Cryptographically Secure PRNG (CSPRNG),
                  Approximate Entropy (ApEn), Secure Hash Algorithm (SHA-1), 
                  Blum-Blum-Shub (BBS) CSPRNG, Micali-Schnorr (MS) CSPRNG, 

   Source      :  David Banks, Elaine Barker, James Dray, Allen Heckert, 
                  Stefan Leigh, Mark Levenson, James Nechvatal, Andrew Rukhin, 
                  Miles Smid, Juan Soto, Mark Vangel, and San Vo.

   Technical
   Assistance  :  Larry Bassham, Ron Boisvert, James Filliben, Daniel Lozier,
                  and Bert Rust.

   Warning     :  Portability Issues.

   Limitation  :  Amount of memory allocated for workspace.

   Restrictions:  Permission to use, copy, and modify this software without 
                  fee is hereby granted, provided that this entire notice is 
                  included in all copies of any software which is or includes
                  a copy or modification of this software and in all copies 
                  of the supporting documentation for such software.
   -------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/decls.h"
#include "../include/cephes.h"  
#include "../include/utilities.h"

void	partitionResultFile(int numOfFiles, int numOfSequences, int option, int testNameID);
void	postProcessResults(int option);
int		computeMetrics(const float *values, int test);

/* p-values of each test's data<k>.txt files, numOfBitStreams per file */
static float	*partitioned[NUMOFTESTS+1];

int
main(int argc, char *argv[])
{
	int		i;
	int		option;			/* TEMPLATE LENGTH/STREAM LENGTH/GENERATOR*/
	char	*streamFile;	/* STREAM FILENAME     */
	

	if ( argc != 2 ) {
		printf("Usage: %s <stream length>\n", argv[0]);
		printf("   <stream length> is the length of the individual bit stream(s) to be processed\n");
		return 0;
	}

	tp.n = atoi(argv[1]);
	tp.blockFrequencyBlockLength = 128;
	tp.nonOverlappingTemplateBlockLength = 9;
	tp.overlappingTemplateBlockLength = 9;
	tp.approximateEntropyBlockLength = 10;
	tp.serialBlockLength = 16;
	tp.linearComplexitySequenceLength = 500;
	tp.numOfBitStreams = 1;
	option = generatorOptions(&streamFile);
	chooseTests();
	fixParameters();
	openOutputStreams(option);
	invokeTestSuite(option, streamFile);
	fclose(freqfp);
	for( i=1; i<=NUMOFTESTS; i++ ) {
		if ( stats[i] != NULL )
			fclose(stats[i]);
		if ( results[i] != NULL )
			fclose(results[i]);
	}
	if ( (testVector[0] == 1) || (testVector[TEST_CUSUM] == 1) ) 
		partitionResultFile(2, tp.numOfBitStreams, option, TEST_CUSUM);
	if ( (testVector[0] == 1) || (testVector[TEST_NONPERIODIC] == 1) ) 
		partitionResultFile(MAXNUMOFTEMPLATES, tp.numOfBitStreams, option, TEST_NONPERIODIC);
	if ( (testVector[0] == 1) || (testVector[TEST_RND_EXCURSION] == 1) )
		partitionResultFile(8, tp.numOfBitStreams, option, TEST_RND_EXCURSION);
	if ( (testVector[0] == 1) || (testVector[TEST_RND_EXCURSION_VAR] == 1) )
		partitionResultFile(18, tp.numOfBitStreams, option, TEST_RND_EXCURSION_VAR);
	if ( (testVector[0] == 1) || (testVector[TEST_SERIAL] == 1) )
		partitionResultFile(2, tp.numOfBitStreams, option, TEST_SERIAL);
	fprintf(summary, "------------------------------------------------------------------------------\n");
	fprintf(summary, "RESULTS FOR THE UNIFORMITY OF P-VALUES AND THE PROPORTION OF PASSING SEQUENCES\n");
	fprintf(summary, "------------------------------------------------------------------------------\n");
	fprintf(summary, "   generator is <%s>\n", streamFile);
	fprintf(summary, "------------------------------------------------------------------------------\n");
	fprintf(summary, " C1  C2  C3  C4  C5  C6  C7  C8  C9 C10  P-VALUE  PROPORTION  STATISTICAL TEST\n");
	fprintf(summary, "------------------------------------------------------------------------------\n");
	postProcessResults(option);
	fclose(summary);

	return 1;
}

static void
fileNotFound(const char *name)
{
	printf("%s", name);
	printf(" -- file not found. Exiting program.\n");
	exit(-1);
}

static char *
readTextFile(const char *name)
{
	FILE	*fp;
	char	*text = NULL, *grown;
	size_t	size = 0, used = 0, got;

	if ( (fp = fopen(name, "r")) == NULL )
		return NULL;
	do {
		if ( size - used < 65536 ) {
			size = 2*size + 65536;
			if ( (grown = (char *)realloc(text, size + 1)) == NULL ) {
				free(text);
				fclose(fp);
				return NULL;
			}
			text = grown;
		}
		got = fread(text + used, 1, size - used, fp);
		used += got;
	} while ( got > 0 );
	fclose(fp);
	text[used] = '\0';

	return text;
}

/* One fscanf(fp, "%f", c) on the text at *pos; when nothing can be read c keeps its value */
static void
scanValue(char **pos, float *c)
{
	char	*end;
	float	value = strtof(*pos, &end);

	if ( end != *pos ) {
		*c = value;
		*pos = end;
	}
}

void
partitionResultFile(int numOfFiles, int numOfSequences, int option, int testNameID)
{ 
	int		i, num, len;
	long	k, total = (long)numOfFiles * numOfSequences;
	size_t	used;
	float	c = 0, *read, *values;
	char	resultsDir[200], s[200], *text, *pos, *buffer;
	FILE	*fp;
	
	sprintf(resultsDir, "experiments/%s/%s/results.txt", generatorDir[option], testNames[testNameID]);
	if ( (text = readTextFile(resultsDir)) == NULL )
		fileNotFound(resultsDir);
	if ( ((read = (float *)malloc((total+1)*sizeof(float))) == NULL) ||
		 ((values = (float *)malloc((total+1)*sizeof(float))) == NULL) ||
		 ((buffer = (char *)malloc((size_t)numOfSequences*64 + 1)) == NULL) ) {
		printf("Final Analysis Report aborted due to insufficient workspace\n");
		exit(-1);
	}
	
	/* the values go out in turn to data1 .. data<numOfFiles>, one round per sequence */
	pos = text;
	for ( k=0; k<total; k++ ) {
		scanValue(&pos, &c);
		read[k] = c;
	}
	for ( i=0; i<numOfFiles; i++ ) {
		used = 0;
		for ( num=0; num<numOfSequences; num++ ) {
			len = sprintf(buffer + used, "%f\n", read[(long)num*numOfFiles + i]);
			values[(long)i*numOfSequences + num] = strtof(buffer + used, NULL);
			used += len;
		}
		sprintf(s, "experiments/%s/%s/data%d.txt", generatorDir[option], testNames[testNameID], i+1);
		if ( (fp = fopen(s, "w")) == NULL )
			fileNotFound(s);
		fwrite(buffer, 1, used, fp);
		fclose(fp);
	}
	free(text);
	free(read);
	free(buffer);
	free(partitioned[testNameID]);
	partitioned[testNameID] = values;
}

void
postProcessResults(int option)
{
	int		i, j, k, randomExcursionSampleSize, generalSampleSize;
	int		passRate, case1, case2, numOfFiles = 2;
	double	p_hat;
	char	s[200], *text, *pos;
	float	c, *values;
	
	for ( i=1; i<=NUMOFTESTS; i++ ) {		// FOR EACH TEST
		if ( testVector[i] ) {
			// SPECIAL CASES -- HANDLING MULTIPLE FILES FOR A SINGLE TEST
			if ( ((i == TEST_CUSUM) && testVector[TEST_CUSUM] ) ||
				 ((i == TEST_NONPERIODIC) && testVector[TEST_NONPERIODIC] ) ||
				 ((i == TEST_RND_EXCURSION) && testVector[TEST_RND_EXCURSION]) ||
				 ((i == TEST_RND_EXCURSION_VAR) && testVector[TEST_RND_EXCURSION_VAR]) || 
				 ((i == TEST_SERIAL) && testVector[TEST_SERIAL]) ) {
				
				if ( (i == TEST_NONPERIODIC) && testVector[TEST_NONPERIODIC] )  
					numOfFiles = MAXNUMOFTEMPLATES;
				else if ( (i == TEST_RND_EXCURSION) && testVector[TEST_RND_EXCURSION] ) 
					numOfFiles = 8;
				else if ( (i == TEST_RND_EXCURSION_VAR) && testVector[TEST_RND_EXCURSION_VAR] ) 
					numOfFiles = 18;
				else
					numOfFiles = 2;
				for ( k=0; k<numOfFiles; k++ ) {
					if ( partitioned[i] == NULL ) {
						sprintf(s, "experiments/%s/%s/data%d.txt", generatorDir[option], testNames[i], k+1);
						fileNotFound(s);
					}
					values = partitioned[i] + (long)k*tp.numOfBitStreams;
					if ( (i == TEST_RND_EXCURSION) || (i == TEST_RND_EXCURSION_VAR) )
						randomExcursionSampleSize = computeMetrics(values, i);
					else
						generalSampleSize = computeMetrics(values, i);
				}
			}
			else {
				sprintf(s, "experiments/%s/%s/results.txt", generatorDir[option], testNames[i]);
				if ( ((text = readTextFile(s)) == NULL) ||
					 ((values = (float *)malloc((tp.numOfBitStreams+1)*sizeof(float))) == NULL) )
					fileNotFound(s);
				c = 0;
				pos = text;
				for ( j=0; j<tp.numOfBitStreams; j++ ) {
					scanValue(&pos, &c);
					values[j] = c;
				}
				generalSampleSize = computeMetrics(values, i);
				free(values);
				free(text);
			}
		}
	}

	fprintf(summary, "\n\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
	case1 = 0;
	case2 = 0;
	if ( testVector[TEST_RND_EXCURSION] || testVector[TEST_RND_EXCURSION_VAR] ) 
		case2 = 1;
	for ( i=1; i<=NUMOFTESTS; i++ ) {
		if ( testVector[i] && (i != TEST_RND_EXCURSION) && (i != TEST_RND_EXCURSION_VAR) ) {
			case1 = 1;
			break;
		}
	}
	p_hat = 1.0 - ALPHA;
	if ( case1 ) {
		if ( generalSampleSize == 0 ) {
			fprintf(summary, "The minimum pass rate for each statistical test with the exception of the\n");
			fprintf(summary, "random excursion (variant) test is undefined.\n\n");
		}
		else {
			passRate = (p_hat - 3.0 * sqrt((p_hat*ALPHA)/generalSampleSize)) * generalSampleSize;
			fprintf(summary, "The minimum pass rate for each statistical test with the exception of the\n");
			fprintf(summary, "random excursion (variant) test is approximately = %d for a\n", generalSampleSize ? passRate : 0);
			fprintf(summary, "sample size = %d binary sequences.\n\n", generalSampleSize);
		}
	}
	if ( case2 ) {
		if ( randomExcursionSampleSize == 0 )
			fprintf(summary, "The minimum pass rate for the random excursion (variant) test is undefined.\n\n");
		else {
			passRate = (p_hat - 3.0 * sqrt((p_hat*ALPHA)/randomExcursionSampleSize)) * randomExcursionSampleSize;
			fprintf(summary, "The minimum pass rate for the random excursion (variant) test\n");
			fprintf(summary, "is approximately = %d for a sample size = %d binary sequences.\n\n", passRate, randomExcursionSampleSize);
		}
	}
	fprintf(summary, "For further guidelines construct a probability table using the MAPLE program\n");
	fprintf(summary, "provided in the addendum section of the documentation.\n");
	fprintf(summary, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
}

/* values: the numOfBitStreams p-values of one results.txt or data<k>.txt file */
int
computeMetrics(const float *values, int test)
{
	int		j, pos, count, passCount, sampleSize, expCount, proportion_threshold_min, proportion_threshold_max;
	int		freqPerBin[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double	*A, chi2, uniformity, p_hat;
	float	c;
	
	if ( (A = (double *)calloc(tp.numOfBitStreams+1, sizeof(double))) == NULL ) {
		printf("Final Analysis Report aborted due to insufficient workspace\n");
		return 0;
	}
	
	/* Compute Metric 1: Proportion of Passing Sequences */
	
	count = 0; 		
	sampleSize = tp.numOfBitStreams;
	
	if ( (test == TEST_RND_EXCURSION) || (test == TEST_RND_EXCURSION_VAR) ) { /* Special Case: Random Excursion Tests */
		for ( j=0; j<sampleSize; j++ ) {
			c = values[j];
			if ( c > 0.000000 )
				A[count++] = c;
		}
		
		sampleSize = count;
		count = 0;
		for ( j=0; j<sampleSize; j++ )
			if ( A[j] < ALPHA )
				count++;
	}
	else {
		for ( j=0; j<sampleSize; j++ ) {
			c = values[j];
			if ( c < ALPHA )
				count++;
			A[j] = c;
		}
	}
	if ( sampleSize == 0 )
		passCount = 0;
	else
		passCount = sampleSize - count;
	
	p_hat = 1.0 - ALPHA;
	proportion_threshold_max = (p_hat + 3.0 * sqrt((p_hat*ALPHA)/sampleSize)) * sampleSize;
	proportion_threshold_min = (p_hat - 3.0 * sqrt((p_hat*ALPHA)/sampleSize)) * sampleSize;
	
	/* Compute Metric 2: Histogram (the bin counts do not depend on the order) */
	
	for ( j=0; j<sampleSize; j++ ) {
		pos = (int)floor(A[j]*10);
		if ( pos == 10 )
			pos--;
		freqPerBin[pos]++;
	}
	chi2 = 0.0;
	expCount = sampleSize/10;
	if ( expCount == 0 )
		uniformity = 0.0;
	else {
		for ( j=0; j<10; j++ )
			chi2 += pow(freqPerBin[j]-expCount, 2)/expCount;
		uniformity = cephes_igamc(9.0/2.0, chi2/2.0);
	}
	
	for ( j=0; j<10; j++ )			/* DISPLAY RESULTS */
		fprintf(summary, "%3d ", freqPerBin[j]);
	
	if ( expCount == 0 )
		fprintf(summary, "    ----    ");
	else if ( uniformity < 0.0001 )
		fprintf(summary, " %8.6f * ", uniformity);
	else
		fprintf(summary, " %8.6f   ", uniformity);
	
	if ( sampleSize == 0 )
		fprintf(summary, " ------     %s\n", testNames[test]);
	//	else if ( proportion < 0.96 )
	else if ( (passCount < proportion_threshold_min) || (passCount > proportion_threshold_max))
		fprintf(summary, "%4d/%-4d *  %s\n", passCount, sampleSize, testNames[test]);
	else
		fprintf(summary, "%4d/%-4d    %s\n", passCount, sampleSize, testNames[test]);
	
	free(A);
	
	return sampleSize;
}
