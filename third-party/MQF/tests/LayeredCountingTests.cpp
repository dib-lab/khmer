#include "../gqf.h"
#include "../LayeredMQF.h"
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>
#include<iostream>
#include "../catch.hpp"
using namespace std;





TEST_CASE( "simple counting test(layered)","[layered]" ) {
  //except first item is inserted 5 times to full test _insert1
  layeredMQF qf;
  int counter_size=2;
  uint64_t qbits=5;
  uint64_t singleQbits=6;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  uint64_t count,fixed_counter;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  layeredMQF_init(&qf,(1ULL<<singleQbits) ,(1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);



  for(uint64_t i=0;i<=10;i++){
    layeredMQF_insert(&qf,100,1,false,false);
    count = layeredMQF_count_key(&qf, 100);
    //fixed_counter=qf_get_fixed_counter(&qf,100);
    INFO("Counter = "<<count<<" fixed counter = "<<fixed_counter)
    CHECK(count == (1+i));
  }


  layeredMQF_insert(&qf,1500,50,false,false);

  count = layeredMQF_count_key(&qf, 1500);
  //  fixed_counter=qf_get_fixed_counter(&qf,1500);
  INFO("Counter = "<<count<<" fixed counter = "<<fixed_counter)
  CHECK(count == (50));

  layeredMQF_insert(&qf,1600,60,false,false);
  count = layeredMQF_count_key(&qf, 1600);
  //  fixed_counter=qf_get_fixed_counter(&qf,1600);
  INFO("Counter = "<<count<<" fixed counter = "<<fixed_counter)
  CHECK(count == (60));


  layeredMQF_insert(&qf,2000,4000,false,false);
  count = layeredMQF_count_key(&qf, 2000);
  //  fixed_counter=qf_get_fixed_counter(&qf,2000);
  INFO("Counter = "<<count<<" fixed counter = "<<fixed_counter)
  CHECK(count == (4000));

}

TEST_CASE( "Maximum count(layered)","[layered]" ) {
  layeredMQF qf;
  int counter_size=4;
  srand (1);
  uint64_t qbits=5;
  uint64_t singleQbits=6;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  layeredMQF_init(&qf,(1ULL<<singleQbits), (1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);
  qf.secondLayer->metadata->maximum_count=10;
  layeredMQF_insert(&qf,100,100000,false,false);
  uint64_t count = layeredMQF_count_key(&qf, 100);
  CHECK(count==10);

  layeredMQF_insert(&qf,150,8,false,false);
  layeredMQF_insert(&qf,150,8,false,false);
  count = layeredMQF_count_key(&qf, 150);
  CHECK(count==10);

}
//
TEST_CASE( "Big count(layered)","[layered]" ) {
  layeredMQF qf;
  int counter_size=4;
  srand (1);
  uint64_t qbits=5;
  uint64_t singleQbits=6;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  layeredMQF_init(&qf,(1ULL<<singleQbits), (1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);
  layeredMQF_insert(&qf,100,100000,false,false);
  uint64_t count = layeredMQF_count_key(&qf, 100);

  CHECK(count==100000);

}

TEST_CASE( "Inserting items( repeated 1 time) in cqf(90% load factor )(layered)" ,"[layered]") {
  //except first item is inserted 5 times to full test _insert1
  layeredMQF qf;
  int counter_size=2;
  uint64_t qbits=15;
  uint64_t singleQbits=16;
  uint64_t num_hash_bits=qbits+9;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  layeredMQF_init(&qf, (1ULL<<singleQbits) ,(1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);

  uint64_t nvals = (1ULL<<qbits)*2;
  uint64_t *vals;
  vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
  for(uint64_t i=0;i<nvals;i++)
  {
    vals[i]=rand();
    vals[i]=(vals[i]<<32)|rand();
    vals[i]=vals[i]%(qf.firstLayer_singletons->metadata->range);
  }

  double loadFactor=(double)qf.firstLayer_singletons->metadata->noccupied_slots/(double)qf.firstLayer_singletons->metadata->nslots;
  uint64_t insertedItems=0;

    layeredMQF_insert(&qf,vals[0],1,false,false);
    layeredMQF_insert(&qf,vals[0],1,false,false);
    layeredMQF_insert(&qf,vals[0],1,false,false);
    layeredMQF_insert(&qf,vals[0],1,false,false);
  // for(uint64_t i=0;i<32;i++)
  // {
  //   cout<<get_fixed_counter(&qf,i)<<"-";
  // }
  //cout<<endl;
  while(loadFactor<0.9){

    layeredMQF_insert(&qf,vals[insertedItems],1,false,false);
    // for(uint64_t i=0;i<32;i++)
    // {
    //   cout<<get_fixed_counter(&qf,i)<<"-";
    // }
    // cout<<endl;
    insertedItems++;
    loadFactor=(double)qf.firstLayer_singletons->metadata->noccupied_slots/(double)qf.firstLayer_singletons->metadata->nslots;
  }
  INFO("Inserted Items = "<<insertedItems);

  uint64_t count;
  //INFO("Fixed counter = "<<qf_get_fixed_counter(&qf,vals[0]));
  count = layeredMQF_count_key(&qf, vals[0]);
  CHECK(count >= 5);

  for(uint64_t i=1;i<insertedItems;i++)
  {
    count = layeredMQF_count_key(&qf, vals[i]);
    CHECK(count >= 1);
  }
  // layeredMQFIterator qfi;
  // layeredMQF_iterator(&qf, &qfi, 0);
  // do {
  //   uint64_t key, value, count;
  //   qfi_get(&qfi, &key, &value, &count);
  //   count=layeredMQF_count_key(&qf, key);
  //   if(key==vals[0]){
  //     CHECK(count >= 5);
  //   }
  //   else{
  //     CHECK(count >= 1);
  //   }
  //
  // } while(!qfi_next(&qfi));

  layeredMQF_destroy(&qf);

}
//
//
// TEST_CASE( "Inserting items( repeated 50 times) in cqf(90% load factor )" ) {
//   layeredMQF qf;
//   int counter_size=4;
//   uint64_t qbits=15;
//   uint64_t num_hash_bits=qbits+8;
//   uint64_t maximum_count=(1ULL<<counter_size)-1;
//   INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
//   layeredMQF_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);
//
//   uint64_t nvals = (1ULL<<qbits);
//   uint64_t *vals;
//   vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
//   for(uint64_t i=0;i<nvals;i++)
//   {
//     vals[i]=rand();
//     vals[i]=(vals[i]<<32)|rand();
//     vals[i]=vals[i]%(qf.metadata->range);
//   }
//   double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
//   uint64_t insertedItems=0;
//   uint64_t count;
//   while(loadFactor<0.9){
//     layeredMQF_insert(&qf,vals[insertedItems],50,false,false);
//
//     insertedItems++;
//     loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
//   }
//
//   for(uint64_t i=0;i<insertedItems;i++)
//   {
//     count = layeredMQF_count_key(&qf, vals[i]);
//     CHECK(count >= 50);
//   }
//   layeredMQFIterator qfi;
//   layeredMQF_iterator(&qf, &qfi, 0);
//   do {
//     uint64_t key, value, count;
//     qfi_get(&qfi, &key, &value, &count);
//     count=layeredMQF_count_key(&qf, key);
//     CHECK(count >= 50);
//   } while(!qfi_next(&qfi));
//
//   layeredMQF_destroy(&qf);
//
// }
//
//
TEST_CASE( "Inserting items( repeated 1-1000 times) in cqf(90% load factor )(layered)","[layered]" ) {
  layeredMQF qf;
  int counter_size=4;
  srand (1);
  uint64_t qbits=15;
  uint64_t singleQbits=16;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  layeredMQF_init(&qf,(1ULL<<singleQbits) ,(1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);

  uint64_t nvals = (1ULL<<singleQbits);
  //uint64_t nvals = 3;
  uint64_t *vals;
  uint64_t *nRepetitions;
  vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
  nRepetitions= (uint64_t*)malloc(nvals*sizeof(nRepetitions[0]));
  uint64_t count;

  for(uint64_t i=0;i<nvals;i++)
  {
    uint64_t newvalue=0;
    while(newvalue==0){
      newvalue=rand();
      newvalue=(newvalue<<32)|rand();
      newvalue=newvalue%(qf.secondLayer->metadata->range);
      for(uint64_t j=0;j<i;j++)
      {
        if(vals[j]==newvalue)
        {
          newvalue=0;
          break;
        }
      }
    }
    vals[i]=newvalue;


    nRepetitions[i]=(rand()%1000)+1;
  }
  int loadFactor=layeredMQF_space(&qf);
  uint64_t insertedItems=0;
  while(insertedItems<nvals && loadFactor<90){
  //  printf("inserting %lu count = %lu\n",vals[insertedItems],nRepetitions[insertedItems] );
    INFO("Inserting "<< vals[insertedItems] << " Repeated "<<nRepetitions[insertedItems]);
    layeredMQF_insert(&qf,vals[insertedItems],nRepetitions[insertedItems],false,false);
    //qf_dump(&qf);
    INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);
    count = layeredMQF_count_key(&qf, vals[insertedItems]);
    CHECK(count == nRepetitions[insertedItems]);
    insertedItems++;
    loadFactor=layeredMQF_space(&qf);

  }
  INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);

  for(uint64_t i=0;i<insertedItems;i++)
  {
    count = layeredMQF_count_key(&qf, vals[i]);
    INFO("value = "<<vals[i]<<" Repeated " <<nRepetitions[i]);
    CHECK(count == nRepetitions[i]);
  }

  layeredMQF_destroy(&qf);

}
//
// TEST_CASE( "Migrate" ) {
//   layeredMQF qf,qf2;
//   int counter_size=4;
//   srand (1);
//   uint64_t qbits=16;
//   uint64_t num_hash_bits=qbits+8;
//   uint64_t maximum_count=(1ULL<<counter_size)-1;
//   INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
//   layeredMQF_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);
//   layeredMQF_init(&qf2, (1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);
//
//   uint64_t nvals = (1ULL<<qbits);
//   //uint64_t nvals = 3;
//   uint64_t *vals;
//   uint64_t *nRepetitions;
//   vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
//   nRepetitions= (uint64_t*)malloc(nvals*sizeof(nRepetitions[0]));
//   uint64_t count;
//
//   for(uint64_t i=0;i<nvals;i++)
//   {
//     uint64_t newvalue=0;
//     while(newvalue==0){
//       newvalue=rand();
//       newvalue=(newvalue<<32)|rand();
//       newvalue=newvalue%(qf.metadata->range);
//       for(uint64_t j=0;j<i;j++)
//       {
//         if(vals[j]==newvalue)
//         {
//           newvalue=0;
//           break;
//         }
//       }
//     }
//     vals[i]=newvalue;
//
//
//     nRepetitions[i]=(rand()%1000)+1;
//   }
//   double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
//   uint64_t insertedItems=0;
//   while(insertedItems<nvals && loadFactor<0.9){
//   //  printf("inserting %lu count = %lu\n",vals[insertedItems],nRepetitions[insertedItems] );
//     INFO("Inserting "<< vals[insertedItems] << " Repeated "<<nRepetitions[insertedItems]);
//     layeredMQF_insert(&qf2,vals[insertedItems],nRepetitions[insertedItems],false,false);
//     //qf_dump(&qf);
//     INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);
//     count = layeredMQF_count_key(&qf2, vals[insertedItems]);
//     CHECK(count == nRepetitions[insertedItems]);
//     insertedItems++;
//     loadFactor=(double)qf2.metadata->noccupied_slots/(double)qf.metadata->nslots;
//
//   }
//   INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);
//   layeredMQF_migrate(&qf2,&qf);
//   for(uint64_t i=0;i<insertedItems;i++)
//   {
//     count = layeredMQF_count_key(&qf, vals[i]);
//     INFO("value = "<<vals[i]<<" Repeated " <<nRepetitions[i]);
//     CHECK(count == nRepetitions[i]);
//   }
//
//   layeredMQF_destroy(&qf);
//
// }
//
// TEST_CASE( "Counting Big counters" ){
//   layeredMQF qf;
//   int counter_size=2;
//   srand (1);
//   uint64_t qbits=16;
//   uint64_t num_hash_bits=qbits+8;
//   uint64_t maximum_count=(1ULL<<counter_size)-1;
//   INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
//   layeredMQF_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);
//
//   uint64_t nvals = (1ULL<<qbits);
//   //uint64_t nvals = 3;
//   uint64_t *vals;
//   uint64_t *nRepetitions;
//   vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
//   nRepetitions= (uint64_t*)malloc(nvals*sizeof(nRepetitions[0]));
//   uint64_t count;
//
//   for(uint64_t i=0;i<nvals;i++)
//   {
//     vals[i]=rand();
//     vals[i]=(vals[i]<<32)|rand();
//     vals[i]=vals[i]%(qf.metadata->range);
//
//     nRepetitions[i]=(rand())+1;
//   }
//   double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
//   uint64_t insertedItems=0;
//   while(insertedItems<nvals && loadFactor<0.9){
//   //  printf("inserting %lu count = %lu\n",vals[insertedItems],nRepetitions[insertedItems] );
//     INFO("Inserting "<< vals[insertedItems] << " Repeated "<<nRepetitions[insertedItems]);
//     layeredMQF_insert(&qf,vals[insertedItems],nRepetitions[insertedItems],false,false);
//     //layeredMQF_dump(&qf);
//     INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);
//     count = layeredMQF_count_key(&qf, vals[insertedItems]);
//     CHECK(count >= nRepetitions[insertedItems]);
//     insertedItems++;
//     loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
//
//   }
//   INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);
//
//   for(uint64_t i=0;i<insertedItems;i++)
//   {
//     count = layeredMQF_count_key(&qf, vals[i]);
//     INFO("value = "<<vals[i]<<" Repeated " <<nRepetitions[i]);
//     CHECK(count >= nRepetitions[i]);
//   }
//
//   layeredMQF_destroy(&qf);
//
//
// }
//
// TEST_CASE( "Removing items from cqf(90% load factor )") {
//   layeredMQF qf;
//   int counter_size=2;
//   uint64_t qbits=16;
//   uint64_t num_hash_bits=qbits+8;
//   uint64_t maximum_count=(1ULL<<counter_size)-1;
//   layeredMQF_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);
//
//   uint64_t nvals = (1ULL<<qbits);
//   uint64_t *vals;
//   vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
//   for(uint64_t i=0;i<nvals;i++)
//   {
//     uint64_t newvalue=0;
//     while(newvalue==0){
//       newvalue=rand();
//       newvalue=(newvalue<<32)|rand();
//       newvalue=newvalue%(qf.metadata->range);
//       for(uint64_t j=0;j<i;j++)
//       {
//         if(vals[j]==newvalue)
//         {
//           newvalue=0;
//           break;
//         }
//       }
//     }
//     vals[i]=newvalue;
//   }
//   double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
//   uint64_t insertedItems=0;
//   while(loadFactor<0.9){
//
//     layeredMQF_insert(&qf,vals[insertedItems],50,false,false);
//     insertedItems++;
//     loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
//   }
//   uint64_t count;
//   for(uint64_t i=0;i<insertedItems;i++)
//   {
//     if(i%2==0){
//       count = layeredMQF_count_key(&qf, vals[i]);
//       if(count==100){
//         printf("coubn ==100\n" );
//       }
//     layeredMQF_remove(&qf,vals[i],50,false,false);
//     count = layeredMQF_count_key(&qf, vals[i]);
//     CHECK(count ==0);
//     }
//   }
//   for(uint64_t i=0;i<insertedItems;i++)
//   {
//     count = layeredMQF_count_key(&qf, vals[i]);
//     if(i%2==1){
//     CHECK(count >= 50);
//     }
//     else{
//       if(count!=0){
//         INFO("ERROR "<<vals[i]<<" Not deleted index= "<<i)
//         //printf("%lu not delete at index %lu\n", vals[i],i);
//       }
//       CHECK(count ==0);
//
//     }
//   }
//   layeredMQFIterator qfi;
//   layeredMQF_iterator(&qf, &qfi, 0);
//   do {
//     uint64_t key, value, count;
//     qfi_get(&qfi, &key, &value, &count);
//     count=layeredMQF_count_key(&qf, key);
//     CHECK(count >= 50);
//   } while(!qfi_next(&qfi));
//
//   layeredMQF_destroy(&qf);
//
// }
