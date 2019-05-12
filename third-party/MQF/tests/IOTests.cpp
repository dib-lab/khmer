#include "../gqf.h"
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>
#include<iostream>
#include "../catch.hpp"
#include <unordered_map>
using namespace std;

TEST_CASE( "Writing and Reading to/from Disk") {
  QF qf;
  int counter_size=2;
  uint64_t qbits=16;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 3,counter_size,0, true, "", 2038074761);

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

    nRepetitions[i]=(rand()%257)+1;
  }
  double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  uint64_t insertedItems=0;
  while(insertedItems<nvals && loadFactor<0.9){
    qf_insert(&qf,vals[insertedItems],nRepetitions[insertedItems],false,false);
    qf_add_tag(&qf,vals[insertedItems],insertedItems%8);
    count = qf_count_key(&qf, vals[insertedItems]);
    CHECK(count >= nRepetitions[insertedItems]);
    insertedItems++;
    loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  }
  INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);
  INFO("nslots ="<<qf.metadata->nslots);
  qf_serialize(&qf,"tmp.ser");
  qf_destroy(&qf);

  SECTION("Reading using qf_read(mmap)"){
    QF qf2;
    qf_read(&qf2,"tmp.ser");
    INFO("nslots ="<<qf2.metadata->nslots);
    for(uint64_t i=0;i<insertedItems;i++)
    {
      count = qf_count_key(&qf2, vals[i]);
      INFO("value = "<<vals[i]<<" Repeated " <<nRepetitions[i]);
      CHECK(count >= nRepetitions[i]);
      CHECK(qf_get_tag(&qf2,vals[i])== i%8);
    }

    qf_destroy(&qf2);
  }

  SECTION("Reading using deserialize "){
    qf_deserialize(&qf,"tmp.ser");

    for(uint64_t i=0;i<insertedItems;i++)
    {
      count = qf_count_key(&qf, vals[i]);
      INFO("value = "<<vals[i]<<" Repeated " <<nRepetitions[i]);
      CHECK(count >= nRepetitions[i]);
      CHECK(qf_get_tag(&qf,vals[i])== i%8);
    }

    qf_destroy(&qf);
  }



}
TEST_CASE("BLOCK Tag after loading")
{
  QF qf;
  int tag_size=8;
  uint64_t qbits=15;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<tag_size)-1;
  INFO("Counter size = "<<tag_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,3,tag_size, true, "", 2038074761);

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
  unordered_map<uint64_t,uint64_t> tags_Map;
  uint64_t currTag=0;
  while(loadFactor<0.9){
    qf_insert(&qf,vals[insertedItems],50,false,false);
    insertedItems++;
    loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  }

//  cout<<endl;
  QFi qfi;
  qf_iterator(&qf, &qfi, 0);
  do {
    uint64_t key, value, count;
    qfi_get(&qfi, &key, &value, &count);

    char* tag;
    bool flag=qf_getBlockTag_pointer_byItem(&qf,key,tag);
    REQUIRE(flag==true);
    uint64_t block_index = qfi.current/64;
    auto it=tags_Map.find(block_index);
    if(it==tags_Map.end())
    {
      tags_Map[block_index]=currTag++;
    }
    *((uint64_t*)tag)=tags_Map[block_index];

  } while(!qfi_next(&qfi));

  qf_serialize(&qf,"tmp.ser");
  qf_destroy(&qf);

  QF qf2;
  qf_deserialize(&qf2,"tmp.ser");
  qf_iterator(&qf2, &qfi, 0);
  do {
    uint64_t key, value, count;
    qfi_get(&qfi, &key, &value, &count);
    count=qf_count_key(&qf2, key);
    CHECK(count >= 50);
    char* tag;
    bool flag=qf_getBlockTag_pointer_byItem(&qf2,key,tag);
    REQUIRE(flag==true);
    uint64_t block_index = qfi.current/64;
    CHECK(tags_Map[block_index]==*((uint64_t*)tag));
    // count = qf_get_tag(&qf,key);2
    // CHECK(count == key%(maximum_count+1));
  } while(!qfi_next(&qfi));




  qf_destroy(&qf2);

}

TEST_CASE( "MMap test") {
  QF qf;
  int counter_size=2;
  uint64_t qbits=16;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 3,counter_size,0, false, "tmp.ser", 2038074761);

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

    nRepetitions[i]=(rand()%257)+1;
  }
  double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  uint64_t insertedItems=0;
  while(insertedItems<nvals && loadFactor<0.9){
    qf_insert(&qf,vals[insertedItems],nRepetitions[insertedItems],false,false);
    qf_add_tag(&qf,vals[insertedItems],insertedItems%8);
    count = qf_count_key(&qf, vals[insertedItems]);
    CHECK(count >= nRepetitions[insertedItems]);
    insertedItems++;
    loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  }
  INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems<<" "<<vals[0]);

  for(uint64_t i=0;i<insertedItems;i++)
  {
    INFO("Check = "<<vals[i]);
    count = qf_count_key(&qf, vals[i]);
    INFO("value = "<<vals[i]<<" Repeated " <<nRepetitions[i]);
    CHECK(count >= nRepetitions[i]);
    CHECK(qf_get_tag(&qf,vals[i])== i%8);
  }

  qf_destroy(&qf);

}
