/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2010-2015, Michigan State University.
Copyright (C) 2015-2016, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/
#ifndef OXLI_HH
#   define OXLI_HH

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

/* The checker automatically defines this preprocessor name when creating
   the custom attribute: */
#if defined(WITH_CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF_ATTRIBUTE)
#define CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF(typename) \
__attribute__((cpychecker_type_object_for_typedef(typename)))
#else
/* This handles the case where we're compiling with a "vanilla"
   compiler that doesn't supply this attribute: */
#define CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF(typename)
#endif

#define NONCOPYABLE(className)\
private:\
    className(const className&);\
    const className& operator=(const className&)

#include <set>
#include <map>
#include <unordered_map>
#include <queue>
#include <list>
#include <functional>

#include "oxli_exception.hh"

#   define MAX_KCOUNT 255
#   define MAX_BIGCOUNT 65535
#   define DEFAULT_TAG_DENSITY 40   // must be even

#   define MAX_CIRCUM 3		// @CTB remove
#   define CIRCUM_RADIUS 2	// @CTB remove
#   define CIRCUM_MAX_VOL 200	// @CTB remove

#   define SAVED_SIGNATURE "OXLI"
#   define SAVED_FORMAT_VERSION 4
#   define SAVED_COUNTING_HT 1
#   define SAVED_HASHBITS 2
#   define SAVED_TAGS 3
#   define SAVED_STOPTAGS 4
#   define SAVED_SUBSET 5
#   define SAVED_LABELSET 6
#   define SAVED_SMALLCOUNT 7
#   define SAVED_QFCOUNT 8

#   define TRAVERSAL_LEFT 0
#   define TRAVERSAL_RIGHT 1

#   define VERBOSE_REPARTITION 0

#   define MIN( a, b )	(((a) > (b)) ? (b) : (a))
#   define MAX( a, b )	(((a) < (b)) ? (b) : (a))

namespace oxli
{

// largest number we can count up to, exactly. (8 bytes)
typedef unsigned long long int ExactCounterType;

// largest number we're going to hash into. (8 bytes/64 bits/32 nt)
typedef unsigned long long int HashIntoType;
const unsigned char KSIZE_MAX = sizeof(HashIntoType)*4;

// largest size 'k' value for k-mer calculations.  (1 byte/255)
typedef unsigned char WordLength;

typedef unsigned short int BoundedCounterType;

// A single-byte type.
typedef unsigned char Byte;

typedef void (*CallbackFn)(const char * info, void * callback_data,
                           unsigned long long n_reads,
                           unsigned long long other);

typedef unsigned int PartitionID;
typedef std::set<HashIntoType> SeenSet;
typedef std::set<PartitionID> PartitionSet;
typedef std::unordered_map<HashIntoType, PartitionID*> PartitionMap;
typedef std::unordered_map<PartitionID, PartitionID*> PartitionPtrMap;
typedef std::unordered_map<PartitionID, SeenSet*> PartitionsToTagsMap;
typedef std::set<PartitionID *> PartitionPtrSet;
typedef std::unordered_map<PartitionID, PartitionPtrSet*> ReversePartitionMap;
typedef std::queue<HashIntoType> NodeQueue;
typedef std::unordered_map<PartitionID, PartitionID*> PartitionToPartitionPMap;
typedef std::unordered_map<HashIntoType, unsigned int> TagCountMap;
typedef std::unordered_map<PartitionID, unsigned int> PartitionCountMap;
typedef std::map<unsigned long long, unsigned long long>
PartitionCountDistribution;

// types used in @camillescott's sparse labeling extension
typedef unsigned long long int Label;
typedef std::unordered_multimap<HashIntoType, Label> TagLabelMap;
typedef std::unordered_multimap<Label, HashIntoType> LabelTagMap;
typedef std::pair<HashIntoType, Label> TagLabelPair;
typedef std::pair<Label, HashIntoType> LabelTagPair;
typedef std::set<Label> LabelSet;
typedef std::set<HashIntoType> TagSet;


template <typename T>
void deallocate_ptr_set(T& s)
{
    for (typename T::iterator i = s.begin(); i != s.end(); ++i) {
        delete *i;
    }
}

class Kmer;
typedef std::queue<Kmer> KmerQueue;
typedef std::set<Kmer> KmerSet;

// A function which takes a Kmer and returns true if it
// is to be filtered / ignored
typedef std::function<bool (const Kmer&)> KmerFilter;
typedef std::list<KmerFilter> KmerFilterList;
typedef std::vector<std::string> StringVector;
}

#endif // OXLI_HH
