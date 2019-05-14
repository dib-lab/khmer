/*
 * =====================================================================================
 *
 *       Filename:  main_release.c
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  2017-02-04 03:40:58 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (ppandey@cs.stonybrook.edu)
 *                  Rob Johnson (rob@cs.stonybrook.edu)
 *   Organization:  Stony Brook University
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <openssl/rand.h>

#include "gqf.h"

int main(int argc, char **argv)
{
	QF cf;
	QFi cfi;
	uint64_t qbits = 15;
	uint64_t nhashbits = qbits + 8;
	uint64_t nslots = (1ULL << qbits);

	/* Initialise the CQF q=15 r= 8 */
	printf("Constructing CQF q=15 and r=8\n");
	qf_init(&cf, nslots, nhashbits, 0);
        printf("Constructing CQF with q=15 and r=8 was done successfully\n");

	nhashbits=qbits+3;
	/* Initialise the CQF q=15 r= 3 */
	printf("Constructing CQF q=15 and r=3\n");
	qf_init(&cf, nslots, nhashbits, 0);
        printf("Constructing CQF with q=15 and r=3 was done successfully\n");

}
