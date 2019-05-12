#include "../gqf.h"
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>
#include<iostream>
#include <unordered_map>
#include "../catch.hpp"
using namespace std;





TEST_CASE( "simple counting test" ) {
  //except first item is inserted 5 times to full test _insert1
  QF qf;
  int counter_size=2;
  uint64_t qbits=5;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  uint64_t count,fixed_counter;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size,0, true, "", 2038074761);



  for(uint64_t i=0;i<=10;i++){
    qf_insert(&qf,100,1,false,false);
    count = qf_count_key(&qf, 100);
    //fixed_counter=qf_get_fixed_counter(&qf,100);
    INFO("Counter = "<<count<<" fixed counter = "<<fixed_counter)
    CHECK(count == (1+i));
  }


  qf_insert(&qf,1500,50,false,false);

  count = qf_count_key(&qf, 1500);
  //  fixed_counter=qf_get_fixed_counter(&qf,1500);
  INFO("Counter = "<<count<<" fixed counter = "<<fixed_counter)
  CHECK(count == (50));

  qf_insert(&qf,1600,60,false,false);
  count = qf_count_key(&qf, 1600);
  //  fixed_counter=qf_get_fixed_counter(&qf,1600);
  INFO("Counter = "<<count<<" fixed counter = "<<fixed_counter)
  CHECK(count == (60));


  qf_insert(&qf,2000,4000,false,false);
  count = qf_count_key(&qf, 2000);
  //  fixed_counter=qf_get_fixed_counter(&qf,2000);
  INFO("Counter = "<<count<<" fixed counter = "<<fixed_counter)
  CHECK(count == (4000));

}

TEST_CASE( "Maximum count" ) {
  QF qf;
  int counter_size=4;
  srand (1);
  uint64_t qbits=5;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size,0, true, "", 2038074761);
  qf.metadata->maximum_count=10;
  qf_insert(&qf,100,100000,false,false);
  uint64_t count = qf_count_key(&qf, 100);
  CHECK(count==10);

  qf_insert(&qf,150,8,false,false);
  qf_insert(&qf,150,8,false,false);
  count = qf_count_key(&qf, 150);
  CHECK(count==10);

}

TEST_CASE( "Big count" ) {
  QF qf;
  int counter_size=4;
  srand (1);
  uint64_t qbits=5;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size,0, true, "", 2038074761);
  qf_insert(&qf,100,100000,false,false);
  uint64_t count = qf_count_key(&qf, 100);

  CHECK(count==100000);

}
TEST_CASE( "Inserting items( repeated 1 time) in cqf(90% load factor )" ) {
  //except first item is inserted 5 times to full test _insert1
  QF qf;
  int counter_size=2;
  uint64_t qbits=15;
  uint64_t num_hash_bits=qbits+9;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size,0, true, "", 2038074761);

  uint64_t nvals = (1ULL<<qbits)*2;
  uint64_t *vals;
  vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
  for(uint64_t i=0;i<nvals;i++)
  {
    vals[i]=rand();
    vals[i]=(vals[i]<<32)|rand();
    vals[i]=vals[i]%(qf.metadata->range);
  }

  double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  uint64_t insertedItems=0;

    qf_insert(&qf,vals[0],1,false,false);
    qf_insert(&qf,vals[0],1,false,false);
    qf_insert(&qf,vals[0],1,false,false);
    qf_insert(&qf,vals[0],1,false,false);
  // for(uint64_t i=0;i<32;i++)
  // {
  //   cout<<get_fixed_counter(&qf,i)<<"-";
  // }
  //cout<<endl;
  while(loadFactor<0.9){

    qf_insert(&qf,vals[insertedItems],1,false,false);
    // for(uint64_t i=0;i<32;i++)
    // {
    //   cout<<get_fixed_counter(&qf,i)<<"-";
    // }
    // cout<<endl;
    insertedItems++;
    loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  }
  INFO("Inserted Items = "<<insertedItems);

  uint64_t count;
  //INFO("Fixed counter = "<<qf_get_fixed_counter(&qf,vals[0]));
  count = qf_count_key(&qf, vals[0]);
  CHECK(count >= 5);

  for(uint64_t i=1;i<insertedItems;i++)
  {
    count = qf_count_key(&qf, vals[i]);
    CHECK(count >= 1);
  }
  QFi qfi;
  qf_iterator(&qf, &qfi, 0);
  do {
    uint64_t key, value, count;
    qfi_get(&qfi, &key, &value, &count);
    count=qf_count_key(&qf, key);
    if(key==vals[0]){
      CHECK(count >= 5);
    }
    else{
      CHECK(count >= 1);
    }

  } while(!qfi_next(&qfi));

  qf_destroy(&qf);

}


TEST_CASE( "Inserting items( repeated 50 times) in cqf(90% load factor )" ) {
  QF qf;
  int counter_size=4;
  uint64_t qbits=15;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size,0, true, "", 2038074761);

  uint64_t nvals = (1ULL<<qbits);
  uint64_t *vals;
  vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
  for(uint64_t i=0;i<nvals;i++)
  {
    vals[i]=rand();
    vals[i]=(vals[i]<<32)|rand();
    vals[i]=vals[i]%(qf.metadata->range);
  }
  double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  uint64_t insertedItems=0;
  uint64_t count;
  while(loadFactor<0.9){
    qf_insert(&qf,vals[insertedItems],50,false,false);

    insertedItems++;
    loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  }

  for(uint64_t i=0;i<insertedItems;i++)
  {
    count = qf_count_key(&qf, vals[i]);
    CHECK(count >= 50);
  }
  QFi qfi;
  qf_iterator(&qf, &qfi, 0);
  do {
    uint64_t key, value, count;
    qfi_get(&qfi, &key, &value, &count);
    count=qf_count_key(&qf, key);
    CHECK(count >= 50);
  } while(!qfi_next(&qfi));

  qf_destroy(&qf);

}


TEST_CASE( "Inserting items( repeated 1-1000 times) in cqf(90% load factor )" ) {
  QF qf;
  int counter_size=4;
  srand (1);
  uint64_t qbits=16;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size,0, true, "", 2038074761);

  uint64_t nvals = (1ULL<<qbits);
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
      newvalue=newvalue%(qf.metadata->range);
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
  double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  uint64_t insertedItems=0;
  while(insertedItems<nvals && loadFactor<0.9){
  //  printf("inserting %lu count = %lu\n",vals[insertedItems],nRepetitions[insertedItems] );
    INFO("Inserting "<< vals[insertedItems] << " Repeated "<<nRepetitions[insertedItems]);
    qf_insert(&qf,vals[insertedItems],nRepetitions[insertedItems],false,false);
    //qf_dump(&qf);
    INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);
    count = qf_count_key(&qf, vals[insertedItems]);
    CHECK(count == nRepetitions[insertedItems]);
    insertedItems++;
    loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;

  }
  INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);

  for(uint64_t i=0;i<insertedItems;i++)
  {
    count = qf_count_key(&qf, vals[i]);
    INFO("value = "<<vals[i]<<" Repeated " <<nRepetitions[i]);
    CHECK(count == nRepetitions[i]);
  }

  qf_destroy(&qf);

}

TEST_CASE( "test kmers order" ) {
  QF qf;
  int counter_size=4;
  srand (1);
  uint64_t qbits=16;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size,0, true, "", 2038074761);

  uint64_t nvals = (1ULL<<qbits);
  //nvals = 10000;
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
      newvalue=newvalue%(qf.metadata->range);
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
  double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  uint64_t insertedItems=0;
  while(insertedItems<nvals && loadFactor<0.9){
  //  printf("inserting %lu count = %lu\n",vals[insertedItems],nRepetitions[insertedItems] );
    INFO("Inserting "<< vals[insertedItems] << " Repeated "<<nRepetitions[insertedItems]);
    qf_insert(&qf,vals[insertedItems],nRepetitions[insertedItems],false,false);
    //qf_dump(&qf);
    INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);
    count = qf_count_key(&qf, vals[insertedItems]);
  //  CHECK(count == nRepetitions[insertedItems]);
    insertedItems++;
    loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->xnslots;

  }
  INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);

  QFi it,itBlock;
  qf_iterator(&qf,&it,0);
  uint64_t key,value,newkey;
  uint64_t currBlock=0;
  uint64_t currCount=0;
  unordered_map<uint64_t,uint64_t> kmerOrders;
  while(!qfi_end(&it))
  {
    qfi_get(&it,&key,&value,&count);
    if(it.current/64 == currBlock)
    {
      kmerOrders[key]=currCount;
      currCount++;
    }
    else{
      currCount=0;
      currBlock=it.current/64;
      kmerOrders[key]=currCount;
      currCount++;
    }
    qfi_next(&it);
  }
  for(uint64_t i=0;i<insertedItems;i++)
  {

    bool res = qfi_find(&qf,&it, vals[i]);
    qfi_get(&it,&key,&value,&count);
    INFO("value = "<<vals[i]<<" Repeated " <<nRepetitions[i]);
    CHECK(count == nRepetitions[i]);

    uint64_t blockId=(it.current/64)*64;
    qfi_firstInBlock(&qf,&it,&itBlock);
    qfi_get(&itBlock,&newkey,&value,&count);
//    cout<<"search for "<<key<<" at "<<it.current<<"\n";
    uint64_t order=0;
    while(newkey!=vals[i])
    {
      qfi_next(&itBlock);
      qfi_get(&itBlock,&newkey,&value,&count);
      //cout<<newkey<<"\n";
      order++;
      if(order>64){
        cout<<"Block "<<blockId<<endl;
        break;
      }
    }
    if(kmerOrders[vals[i]]!=order)
    {
      cout<<it.current<<endl;
    }
    CHECK(kmerOrders[vals[i]]==order);
  }

  qf_destroy(&qf);

}

TEST_CASE( "test get_iterator" ) {
  QF qf;
  int counter_size=4;
  srand (1);
  uint64_t qbits=16;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size,0, true, "", 2038074761);

  uint64_t nvals = (1ULL<<qbits);
  //nvals = 5000;
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
      newvalue=newvalue%(qf.metadata->range);
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
  double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  uint64_t insertedItems=0;
  while(insertedItems<nvals && loadFactor<0.9){
  //  printf("inserting %lu count = %lu\n",vals[insertedItems],nRepetitions[insertedItems] );
    INFO("Inserting "<< vals[insertedItems] << " Repeated "<<nRepetitions[insertedItems]);
    qf_insert(&qf,vals[insertedItems],nRepetitions[insertedItems],false,false);
    //qf_dump(&qf);
    INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);
    count = qf_count_key(&qf, vals[insertedItems]);
  //  CHECK(count == nRepetitions[insertedItems]);
    insertedItems++;
    loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->xnslots;

  }
  INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);

  for(uint64_t i=0;i<insertedItems;i++)
  {
    QFi it,itBlock;
    bool res = qfi_find(&qf,&it, vals[i]);
    uint64_t key,value,count,newkey;
    qfi_get(&it,&key,&value,&count);


    INFO("value = "<<vals[i]<<" Repeated " <<nRepetitions[i]);
    CHECK(count == nRepetitions[i]);

    uint64_t blockId=(it.current/64)*64;
    qfi_firstInBlock(&qf,&it,&itBlock);
    qfi_get(&itBlock,&newkey,&value,&count);
//    cout<<"search for "<<key<<" at "<<it.current<<"\n";
    uint64_t tmp=0;
    while(newkey!=key)
    {
      qfi_next(&itBlock);
      qfi_get(&itBlock,&newkey,&value,&count);
      //cout<<newkey<<"\n";
      tmp++;
      if(tmp>64){
        cout<<"Block "<<blockId<<endl;
        break;
      }
    }
    CHECK(tmp<64);
  }

  qf_destroy(&qf);

}


TEST_CASE( "Migrate" ) {
  QF qf,qf2;
  int counter_size=4;
  srand (1);
  uint64_t qbits=16;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size,0, true, "", 2038074761);
  qf_init(&qf2, (1ULL<<qbits), num_hash_bits, 0,counter_size,0, true, "", 2038074761);

  uint64_t nvals = (1ULL<<qbits);
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
      newvalue=newvalue%(qf.metadata->range);
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
  double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  uint64_t insertedItems=0;
  while(insertedItems<nvals && loadFactor<0.9){
  //  printf("inserting %lu count = %lu\n",vals[insertedItems],nRepetitions[insertedItems] );
    INFO("Inserting "<< vals[insertedItems] << " Repeated "<<nRepetitions[insertedItems]);
    qf_insert(&qf2,vals[insertedItems],nRepetitions[insertedItems],false,false);
    //qf_dump(&qf);
    INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);
    count = qf_count_key(&qf2, vals[insertedItems]);
    CHECK(count == nRepetitions[insertedItems]);
    insertedItems++;
    loadFactor=(double)qf2.metadata->noccupied_slots/(double)qf.metadata->nslots;

  }
  INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);
  qf_migrate(&qf2,&qf);
  for(uint64_t i=0;i<insertedItems;i++)
  {
    count = qf_count_key(&qf, vals[i]);
    INFO("value = "<<vals[i]<<" Repeated " <<nRepetitions[i]);
    CHECK(count == nRepetitions[i]);
  }

  qf_destroy(&qf);

}

TEST_CASE( "Counting Big counters" ){
  QF qf;
  int counter_size=2;
  srand (1);
  uint64_t qbits=16;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size,0, true, "", 2038074761);

  uint64_t nvals = (1ULL<<qbits);
  //uint64_t nvals = 3;
  uint64_t *vals;
  uint64_t *nRepetitions;
  vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
  nRepetitions= (uint64_t*)malloc(nvals*sizeof(nRepetitions[0]));
  uint64_t count;

  for(uint64_t i=0;i<nvals;i++)
  {
    vals[i]=rand();
    vals[i]=(vals[i]<<32)|rand();
    vals[i]=vals[i]%(qf.metadata->range);

    nRepetitions[i]=(rand())+1;
  }
  double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  uint64_t insertedItems=0;
  while(insertedItems<nvals && loadFactor<0.9){
  //  printf("inserting %lu count = %lu\n",vals[insertedItems],nRepetitions[insertedItems] );
    INFO("Inserting "<< vals[insertedItems] << " Repeated "<<nRepetitions[insertedItems]);
    qf_insert(&qf,vals[insertedItems],nRepetitions[insertedItems],false,false);
    //qf_dump(&qf);
    INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);
    count = qf_count_key(&qf, vals[insertedItems]);
    CHECK(count >= nRepetitions[insertedItems]);
    insertedItems++;
    loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;

  }
  INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);

  for(uint64_t i=0;i<insertedItems;i++)
  {
    count = qf_count_key(&qf, vals[i]);
    INFO("value = "<<vals[i]<<" Repeated " <<nRepetitions[i]);
    CHECK(count >= nRepetitions[i]);
  }

  qf_destroy(&qf);


}

TEST_CASE( "Removing items from cqf(90% load factor )") {
  QF qf;
  int counter_size=2;
  uint64_t qbits=16;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size,0, true, "", 2038074761);

  uint64_t nvals = (1ULL<<qbits);
  uint64_t *vals;
  vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
  for(uint64_t i=0;i<nvals;i++)
  {
    uint64_t newvalue=0;
    while(newvalue==0){
      newvalue=rand();
      newvalue=(newvalue<<32)|rand();
      newvalue=newvalue%(qf.metadata->range);
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
  }
  double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  uint64_t insertedItems=0;
  while(loadFactor<0.9){

    qf_insert(&qf,vals[insertedItems],50,false,false);
    insertedItems++;
    loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  }
  uint64_t count;
  for(uint64_t i=0;i<insertedItems;i++)
  {
    if(i%2==0){
      count = qf_count_key(&qf, vals[i]);
      if(count==100){
        printf("coubn ==100\n" );
      }
    qf_remove(&qf,vals[i],50,false,false);
    count = qf_count_key(&qf, vals[i]);
    CHECK(count ==0);
    }
  }
  for(uint64_t i=0;i<insertedItems;i++)
  {
    count = qf_count_key(&qf, vals[i]);
    if(i%2==1){
    CHECK(count >= 50);
    }
    else{
      if(count!=0){
        INFO("ERROR "<<vals[i]<<" Not deleted index= "<<i)
        //printf("%lu not delete at index %lu\n", vals[i],i);
      }
      CHECK(count ==0);

    }
  }
  QFi qfi;
  qf_iterator(&qf, &qfi, 0);
  do {
    uint64_t key, value, count;
    qfi_get(&qfi, &key, &value, &count);
    count=qf_count_key(&qf, key);
    CHECK(count >= 50);
  } while(!qfi_next(&qfi));

  qf_destroy(&qf);

}
