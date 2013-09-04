//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
//

#ifndef KHMER_HH
#   define KHMER_HH

// C standard integer types are used almost ubiquitously.
#   if (__cplusplus >= 201103L)
#	include <cstdint>
#   else
extern "C"
{
#	include <stdint.h>
}
#   endif
// Provide standard limits to all.
#include <climits>
#ifndef SSIZE_MAX
#   define SSIZE_MAX	((ssize_t)(SIZE_MAX / 2))
#endif

// C++ standard exceptions are subclassed almost ubiquitously.
#include <exception>
#include <set>
#include <map>
#include <queue>

#   define VERSION "0.4"

#   define MAX_COUNT 255
#   define MAX_BIGCOUNT 65535
#   define DEFAULT_TAG_DENSITY 40   // must be even

#   define MAX_CIRCUM 3		// @CTB remove
#   define CIRCUM_RADIUS 2	// @CTB remove
#   define CIRCUM_MAX_VOL 200	// @CTB remove

#   define SAVED_FORMAT_VERSION 3
#   define SAVED_COUNTING_HT 1
#   define SAVED_HASHBITS 2
#   define SAVED_TAGS 3
#   define SAVED_STOPTAGS 4
#   define SAVED_SUBSET 5

#   define VERBOSE_REPARTITION 0

#   define MIN( a, b )	(((a) > (b)) ? (b) : (a))
#   define MAX( a, b )	(((a) < (b)) ? (b) : (a))

namespace khmer {
  // largest number we can count up to, exactly. (8 bytes)
  typedef unsigned long long int ExactCounterType;

  // largest number we're going to hash into. (8 bytes/64 bits/32 nt)
  typedef unsigned long long int HashIntoType;

  // largest size 'k' value for k-mer calculations.  (1 byte/255)
  typedef unsigned char WordLength;

  typedef unsigned short int BoundedCounterType;

  // A single-byte type.
  typedef unsigned char Byte;

  typedef void (*CallbackFn)(const char * info, void * callback_data,
			     unsigned long long n_reads,
			     unsigned long long other);

  struct InvalidStreamBuffer : public std:: exception
  { };


  typedef unsigned int PartitionID;
  typedef std::set<HashIntoType> SeenSet;
  typedef std::set<PartitionID> PartitionSet;
  typedef std::map<HashIntoType, PartitionID*> PartitionMap;
  typedef std::map<PartitionID, PartitionID*> PartitionPtrMap;
  typedef std::map<PartitionID, SeenSet*> PartitionsToTagsMap;
  typedef std::set<PartitionID *> PartitionPtrSet;
  typedef std::map<PartitionID, PartitionPtrSet*> ReversePartitionMap;
  typedef std::queue<HashIntoType> NodeQueue;
  typedef std::map<PartitionID, PartitionID*> PartitionToPartitionPMap;
  typedef std::map<HashIntoType, unsigned int> TagCountMap;
  typedef std::map<PartitionID, unsigned int> PartitionCountMap;
  typedef std::map<unsigned long long, unsigned long long> PartitionCountDistribution;
}

#endif // KHMER_HH
