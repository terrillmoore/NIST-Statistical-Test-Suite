#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"  

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                     R A N D O M  E X C U R S I O N S  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
RandomExcursions(int n)
{
	int		i, k, J, x;
	int		stateX[8] = { -4, -3, -2, -1, 1, 2, 3, 4 };
	double	p_value, sum, constraint, nu[6][8];
	double	pi[5][6] = { {0.0000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.0000000000}, 
						 {0.5000000000, 0.25000000000, 0.12500000000, 0.06250000000, 0.03125000000, 0.0312500000},
						 {0.7500000000, 0.06250000000, 0.04687500000, 0.03515625000, 0.02636718750, 0.0791015625},
						 {0.8333333333, 0.02777777778, 0.02314814815, 0.01929012346, 0.01607510288, 0.0803755143},
						 {0.8750000000, 0.01562500000, 0.01367187500, 0.01196289063, 0.01046752930, 0.0732727051} };
	const WalkSummary	*walk = walkSummary();
	
	if ( walk->zeros > MAX(1000, n/100) ) {
		printf("ERROR IN FUNCTION randomExcursions:  EXCEEDING THE MAX NUMBER OF CYCLES EXPECTED\n.");
		return;
	}
	J = walk->zeros;
	if ( walk->S != 0 )
		J++;

	fprintf(stats[TEST_RND_EXCURSION], "\t\t\t  RANDOM EXCURSIONS TEST\n");
	fprintf(stats[TEST_RND_EXCURSION], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_RND_EXCURSION], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_RND_EXCURSION], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_RND_EXCURSION], "\t\t(a) Number Of Cycles (J) = %04d\n", J);
	fprintf(stats[TEST_RND_EXCURSION], "\t\t(b) Sequence Length (n)  = %d\n", n);

	constraint = MAX(0.005*pow(n, 0.5), 500);
	if (J < constraint) {
		fprintf(stats[TEST_RND_EXCURSION], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_RND_EXCURSION], "\t\tWARNING:  TEST NOT APPLICABLE.  THERE ARE AN\n");
		fprintf(stats[TEST_RND_EXCURSION], "\t\t\t  INSUFFICIENT NUMBER OF CYCLES.\n");
		fprintf(stats[TEST_RND_EXCURSION], "\t\t---------------------------------------------\n");
		for(i = 0; i < 8; i++)
			fprintf(results[TEST_RND_EXCURSION], "%f\n", 0.0);
	}
	else {
		fprintf(stats[TEST_RND_EXCURSION], "\t\t(c) Rejection Constraint = %f\n", constraint);
		fprintf(stats[TEST_RND_EXCURSION], "\t\t-------------------------------------------\n");

		for ( k=0; k<6; k++ )
			for ( i=0; i<8; i++ )
				nu[k][i] = walk->cycles[k][i];
		
		for ( i=0; i<8; i++ ) {
			x = stateX[i];
			sum = 0.;
			for ( k=0; k<6; k++ )
				sum += pow(nu[k][i] - J*pi[(int)fabs(x)][k], 2) / (J*pi[(int)fabs(x)][k]);
			p_value = cephes_igamc(2.5, sum/2.0);
			
			if ( isNegative(p_value) || isGreaterThanOne(p_value) )
				fprintf(stats[TEST_RND_EXCURSION], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

			fprintf(stats[TEST_RND_EXCURSION], "%s\t\tx = %2d chi^2 = %9.6f p_value = %f\n",
					p_value < ALPHA ? "FAILURE" : "SUCCESS", x, sum, p_value);
			fprintf(results[TEST_RND_EXCURSION], "%f\n", p_value);
		}
	} 
	fprintf(stats[TEST_RND_EXCURSION], "\n");
}
