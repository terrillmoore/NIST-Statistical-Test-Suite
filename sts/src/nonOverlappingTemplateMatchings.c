#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"  

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          N O N O V E R L A P P I N G  T E M P L A T E  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * The templates of a length are read from their file once and kept for later streams; the
 * reads are the original's (fscanf per bit, fseek past SKIP-1 lines), so the same templates
 * are used. A template that cannot overlap itself is never skipped by the greedy scan, so its
 * count in a block is the number of windows equal to it, and one histogram of the m-bit
 * windows serves all such templates. Any other template takes the greedy scan.
 */
static int			loadedLength = -1, loadedCount, histogramOnly;
static BitSequence	*loaded = NULL;			/* loadedCount templates of loadedLength values */
static unsigned int	*loadedValue = NULL;	/* each template as an m-bit number, first bit on top */
static char			*labels = NULL;			/* each template as printed, then a space */
static int			*labelStart = NULL;

static int
loadTemplates(int m, const char *directory)
{
	int				numOfTemplates[100] = {0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116,
						2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152};
	/*----------------------------------------------------------------------------
	NOTE:  Should additional templates lengths beyond 21 be desired, they must 
	first be constructed, saved into files and then the corresponding 
	number of nonperiodic templates for that file be stored in the m-th 
	position in the numOfTemplates variable.
	----------------------------------------------------------------------------*/
	unsigned int	bit = 0, value;
	FILE			*fp;
	int				jj, k, SKIP, count, total, used;

	if ( m == loadedLength )
		return 1;
	if ( (fp = fopen(directory, "r")) == NULL )
		return 0;
	total = ((m >= 0) && (m < 100)) ? numOfTemplates[m] : 0;
	if ( total < MAXNUMOFTEMPLATES )
		SKIP = 1;
	else
		SKIP = (int)(total/MAXNUMOFTEMPLATES);
	count = MIN(MAXNUMOFTEMPLATES, (int)total/SKIP);
	free(loaded);
	free(loadedValue);
	free(labels);
	free(labelStart);
	loaded = (BitSequence *) malloc((size_t)count*m + 1);
	loadedValue = (unsigned int *) malloc((count + 1)*sizeof(unsigned int));
	labels = (char *) malloc((size_t)count*(4*m + 2) + 1);
	labelStart = (int *) malloc((count + 1)*sizeof(int));
	loadedLength = -1;
	if ( (loaded == NULL) || (loadedValue == NULL) || (labels == NULL) || (labelStart == NULL) ) {
		fclose(fp);
		return 0;
	}
	histogramOnly = (m >= 1) && (m <= 20);
	used = 0;
	for ( jj=0; jj<count; jj++ ) {
		labelStart[jj] = used;
		value = 0;
		for ( k=0; k<m; k++ ) {
			fscanf(fp, "%d", &bit);
			loaded[jj*m+k] = bit;
			used += sprintf(labels + used, "%d", loaded[jj*m+k]);
			if ( loaded[jj*m+k] > 1 )
				histogramOnly = 0;
			value = (value << 1) | (loaded[jj*m+k] & 1);
		}
		used += sprintf(labels + used, " ");
		loadedValue[jj] = value;
		/* a template that overlaps itself: some proper prefix equals the suffix */
		for ( k=1; (k<m) && histogramOnly; k++ )
			if ( (value >> k) == (value & ((1U << (m-k)) - 1)) )
				histogramOnly = 0;
		if ( SKIP > 1 )
			fseek(fp, (long)(SKIP-1)*2*m, SEEK_CUR);
	}
	labelStart[count] = used;
	fclose(fp);
	loadedLength = m;
	loadedCount = count;
	return 1;
}

/* The original greedy scan of one template over block i */
static unsigned int
greedyCount(const BitSequence *sequence, int m, int M, int i)
{
	unsigned int	W_obs = 0;
	long			pos;
	int				j, k, match;

	for ( j=0; j<M-m+1; j++ ) {
		match = 1;
		for ( k=0; k<m; k++ ) {
			pos = (long)i*M + j + k;
			if ( (int)sequence[k] != (int)((epsilonWords[pos >> 6] >> (63 - (pos & 63))) & 1) ) {
				match = 0;
				break;
			}
		}
		if ( match == 1 ) {
			W_obs++;
			j += m-1;
		}
	}
	return W_obs;
}

void
NonOverlappingTemplateMatchings(int m, int n)
{
	static unsigned int	*hist = NULL, *Wj = NULL;
	static double		*term = NULL;
	static long			histSize = 0, WjSize = 0, termSize = 0;
	static int			termLength = -1, termBlock = -1;
	unsigned int	W_obs, *W;
	double			chi2, p_value, lambda, varWj;
	int				i, jj, k, M, N;
	char			directory[100];

	N = 8;
	M = n/N;

	lambda = (M-m+1)/pow(2, m);
	varWj = M*(1.0/pow(2.0, m) - (2.0*m-1.0)/pow(2.0, 2.0*m));
	sprintf(directory, "templates/template%d", m);

	if ( ((isNegative(lambda)) || (isZero(lambda))) || !loadTemplates(m, directory) ) {
		fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
		fprintf(stats[TEST_NONPERIODIC], "\tLambda (%f) not being positive!\n", lambda);
		fprintf(stats[TEST_NONPERIODIC], "\tTemplate file <%s> not existing\n", directory);
		fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
	}
	else {
		fprintf(stats[TEST_NONPERIODIC], "\t\t  NONPERIODIC TEMPLATES TEST\n");
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		fprintf(stats[TEST_NONPERIODIC], "\t\t  COMPUTATIONAL INFORMATION\n");
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		fprintf(stats[TEST_NONPERIODIC], "\tLAMBDA = %f\tM = %d\tN = %d\tm = %d\tn = %d\n", lambda, M, N, m, n);
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		fprintf(stats[TEST_NONPERIODIC], "\t\tF R E Q U E N C Y\n");
		fprintf(stats[TEST_NONPERIODIC], "Template   W_1  W_2  W_3  W_4  W_5  W_6  W_7  W_8    Chi^2   P_value Assignment Index\n");
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");

		if ( (long)loadedCount*N > WjSize ) {
			free(Wj);
			WjSize = (long)loadedCount*N;
			if ( (Wj = (unsigned int *) malloc(WjSize*sizeof(unsigned int))) == NULL ) {
				WjSize = 0;
				return;
			}
		}
		if ( histogramOnly && ((1L << m) > histSize) ) {
			free(hist);
			histSize = 1L << m;
			if ( (hist = (unsigned int *) malloc(histSize*sizeof(unsigned int))) == NULL ) {
				histSize = 0;
				return;
			}
		}
		/* W[jj*N+i]: matches of template jj in block i */
		W = Wj;
		if ( histogramOnly ) {
			memset(hist, 0, (1L << m)*sizeof(unsigned int));
			for ( i=0; i<N; i++ ) {
				windowCounts(hist, (long)i*M, M-m+1, m);
				for ( jj=0; jj<loadedCount; jj++ ) {
					W_obs = hist[loadedValue[jj]];
					for ( k=0; k<i; k++ )
						W_obs -= W[jj*N+k];
					W[jj*N+i] = W_obs;
				}
			}
		}
		else {
			for ( jj=0; jj<loadedCount; jj++ )
				for ( i=0; i<N; i++ )
					W[jj*N+i] = greedyCount(loaded + (long)jj*m, m, M, i);
		}
		
		/* chi-square terms, memoized per count (negative: not yet computed); a count in a
		   block is at most M/m */
		if ( (termLength != m) || (termBlock != M) ) {
			free(term);
			termSize = M/(m > 0 ? m : 1) + 2;
			if ( (term = (double *) malloc(termSize*sizeof(double))) == NULL ) {
				termSize = 0;
				termLength = -1;
				return;
			}
			for ( k=0; k<termSize; k++ )
				term[k] = -1;
			termLength = m;
			termBlock = M;
		}
		for ( jj=0; jj<loadedCount; jj++ ) {
			fwrite(labels + labelStart[jj], 1, labelStart[jj+1] - labelStart[jj], stats[TEST_NONPERIODIC]);
			chi2 = 0.0;                                   /* Compute Chi Square */
			for ( i=0; i<N; i++ ) {
				W_obs = W[jj*N+i];
				if ( m == 10 )
					fprintf(stats[TEST_NONPERIODIC], "%3d  ", W_obs);
				else
					fprintf(stats[TEST_NONPERIODIC], "%4d ", W_obs);
				if ( term[W_obs] < 0 )
					term[W_obs] = pow(((double)W_obs - lambda)/pow(varWj, 0.5), 2);
				chi2 += term[W_obs];
			}
			p_value = cephes_igamc(N/2.0, chi2/2.0);
		
			if ( isNegative(p_value) || isGreaterThanOne(p_value) )
				fprintf(stats[TEST_NONPERIODIC], "\t\tWARNING:  P_VALUE IS OUT OF RANGE.\n");

			fprintf(stats[TEST_NONPERIODIC], "%9.6f %f %s %3d\n", chi2, p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS", jj);
			fprintf(results[TEST_NONPERIODIC], "%f\n", p_value);
		}
	}
	
	fprintf(stats[TEST_NONPERIODIC], "\n");
}
