#include "gqf.h"
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>
#include<iostream>
using namespace std;

 int main(int argc, char const *argv[]) {
  QF qf;
  int counter_size=1;
  srand (1);
  uint64_t qbits=16;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size,0, true, "", 2038074761);

  uint64_t nvals = (1ULL<<qbits);
  //uint64_t nvals = 3;
  uint64_t *vals;
  uint64_t *nRepetitions;
  vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
  nRepetitions= (uint64_t*)malloc(nvals*sizeof(nRepetitions[0]));
  uint64_t count;

  for(int i=0;i<nvals;i++)
  {
    vals[i]=rand();
    vals[i]=(vals[i]<<32)|rand();
    vals[i]=vals[i]%(qf.metadata->range);

    nRepetitions[i]=(rand()%257)+1;
  }

  double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  uint64_t insertedItems=0;
  while(insertedItems<nvals){

    try{
    qf_insert(&qf,vals[insertedItems],nRepetitions[insertedItems],false,false);
    }
    catch (std::exception& e)
  {
    std::cerr << "exception caught: " << e.what() << '\n';
    std::cout<< loadFactor<<endl;
    break;
  }
    //qf_dump(&qf);
    count = qf_count_key(&qf, vals[insertedItems]);
    insertedItems++;
    loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;

  }


  for(int i=0;i<insertedItems;i++)
  {
    count = qf_count_key(&qf, vals[i]);
    if(count < nRepetitions[i] ){
      cout<<"Item  "<<vals[i]<<" is not counted correctly."<<endl;
      return -1;
    }
  }
  cout<<"All kmers are counted successfully"<<endl;
  qf_destroy(&qf);
  return 0;
}
