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
        QF cf,cf1,cf2;
	QFi cfi;
	uint64_t qbits = 18;
	uint64_t small_qbits=qbits;
	qbits+=2;
	uint64_t nhashbits = qbits + 8;
	uint64_t small_nhashbits=small_qbits+8;
	uint64_t nslots = (1ULL << qbits);
	uint64_t small_nslots=(1ULL << small_qbits);
	uint64_t nvals = 250*nslots/1000;
	uint64_t *vals;

	/* Initialise the CQF */
	printf("Initialize first cqf size =%d, hashbits=%d\n", nslots, nhashbits);
	qf_init(&cf, nslots, nhashbits, 0);
	printf("Initialize Second cqf size =%d, hashbits=%d\n",small_nslots,small_nhashbits);
	qf_init(&cf1, small_nslots,small_nhashbits, 0);
	printf("Initialize Third cqf size =%d, hashbits=%d\n",small_nslots,small_nhashbits);
	qf_init(&cf2,small_nslots,small_nhashbits, 0);
	/* Generate random values */
	vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
	RAND_pseudo_bytes((unsigned char *)vals, sizeof(*vals) * nvals);
	for (uint64_t i = 0; i < nvals; i++) {
		vals[i] = (1 * vals[i]);
	}
	printf("Inserting\n");
	/* Insert vals in the CQF */
	for (uint64_t i = 0; i < (nvals*2)/3; i++) {
	  vals[i]=vals[i]%cf1.range;
	  if(i%2==1){
	    qf_insert(&cf2, vals[i], 0, 50);
	  }
	  else{
	    qf_insert(&cf1, vals[i], 0, 50);
	  }

	}
	printf("Merging\n");

	qf_merge(&cf1,&cf2,&cf);
	printf("Inserting again\n");
	for (uint64_t i = (nvals*2)/3; i <nvals; i++) {
	  vals[i]=vals[i]%cf.range;
	    qf_insert(&cf, vals[i], 0, 50);
	  }

	for (uint64_t i = 0; i < nvals; i++) {
		uint64_t count = qf_count_key_value(&cf, vals[i], 0);
		if (count < 50) {
			fprintf(stderr, "failed lookup after insertion for %lx %ld.\n", vals[i],
							count);
			abort();
		}
	}

	/* Initialize an iterator */
	qf_iterator(&cf, &cfi, 0);
	do {
		uint64_t key, value, count;
		qfi_get(&cfi, &key, &value, &count);
		if (qf_count_key_value(&cf, key, 0) < 50) {
			fprintf(stderr, "Failed lookup from A for: %ld. Returned count: %ld\n",
							key, qf_count_key_value(&cf, key, 0));
			abort();
		}
	} while(!qfi_next(&cfi));

	fprintf(stdout, "Validated the CQF.\n");
}
