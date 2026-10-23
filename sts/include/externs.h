
#include "../include/defs.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                   G L O B A L   D A T A  S T R U C T U R E S 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

extern BitSequence	*epsilon;				// BIT STREAM

/* epsilonWords is epsilon packed 64 bits per word: bit i of the stream is bit 63-i%64 of
   word i/64. Bits from tp.n on are zero, and two more zero words follow the last one. */
extern uint64_t		*epsilonWords;

extern TP			tp;						// TEST PARAMETER STRUCTURE
extern FILE			*stats[NUMOFTESTS+1];	// FILE OUTPUT STREAM
extern FILE			*results[NUMOFTESTS+1];	// FILE OUTPUT STREAM
extern FILE			*freqfp;				// FILE OUTPUT STREAM
extern FILE			*summary;				// FILE OUTPUT STREAM
extern int			testVector[NUMOFTESTS+1];

extern char	generatorDir[NUMOFGENERATORS][20];
extern char	testNames[NUMOFTESTS+1][32];
