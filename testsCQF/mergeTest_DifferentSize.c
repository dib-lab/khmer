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
	uint64_t qbits = 17;
	uint64_t qbits1=qbits;
	uint64_t qbits2=qbits+1;
	uint64_t large_qbits=qbits+2;
	uint64_t nhashbits = large_qbits + 8;
	uint64_t nhashbits1 = qbits1 + 8;
	uint64_t nhashbits2 = qbits2 + 8;

	uint64_t nslots = (1ULL << large_qbits);
	uint64_t nslots1=(1ULL << qbits1);
	uint64_t nslots2=(1ULL << qbits2);
	uint64_t nvals = 250*nslots1/1000;
	uint64_t *vals;

	/* Initialise the CQF */
	printf("Initialize first cqf size =%d, hashbits=%d\n", nslots, nhashbits);
	qf_init(&cf, nslots, nhashbits, 0);
	printf("Initialize Second cqf size =%d, hashbits=%d\n",nslots1,nhashbits1);
	qf_init(&cf1, nslots1,nhashbits1, 0);
	printf("Initialize Third cqf size =%d, hashbits=%d\n",nslots2,nhashbits2);
	qf_init(&cf2,nslots2,nhashbits2, 0);
	/* Generate random values */
	vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
	RAND_pseudo_bytes((unsigned char *)vals, sizeof(*vals) * nvals);
	for (uint64_t i = 0; i < nvals; i++) {
		vals[i] = (1 * vals[i]);
	}
	vals[0]=131074;
	printf("Inserting\n");
	/* Insert vals in the CQF */
	for (uint64_t i = 0; i < (nvals*2)/3; i++) {

	  if(i%2==1){
//printf("%d\n",vals[i] );
	    qf_insert(&cf2, vals[i]%cf2.range, 0, 50);
	  }
	  else{
	    qf_insert(&cf1, vals[i]%cf1.range, 0, 50);
	  }

	}

	printf("Merging\n");
	qf_merge(&cf1,&cf2,&cf);

	printf("Inserting again into Big one\n");
	for (uint64_t i = (nvals*2)/3; i <nvals; i++) {
	    qf_insert(&cf, vals[i]%cf.range, 0, 50);
	  }

	for (uint64_t i = 0; i < nvals; i++) {
		uint64_t count = qf_count_key_value(&cf, vals[i]%cf.range, 0);
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
