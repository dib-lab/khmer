#define VERSION "0.3"

#define MAX_COUNT 255
#define DEFAULT_TAG_DENSITY 40		// must be even

#define MAX_CIRCUM 3		// @CTB remove
#define CIRCUM_RADIUS 2		// @CTB remove
#define CIRCUM_MAX_VOL 200	// @CTB remove

#define SAVED_FORMAT_VERSION 1
#define SAVED_COUNTING_HT 1
#define SAVED_HASHBITS 2

#define VERBOSE_REPARTITION 0

namespace khmer {
  // largest number we can count up to, exactly. (8 bytes)
  typedef unsigned long long int ExactCounterType;

  // largest number we're going to hash into. (8 bytes/64 bits/32 nt)
  typedef unsigned long long int HashIntoType;

  // largest size 'k' value for k-mer calculations.  (1 byte/255)
  typedef unsigned char WordLength;

  // largest number we can count up to, approximately. (8 bytes/127).
  // see MAX_COUNT, above.
  typedef unsigned char BoundedCounterType;

  // A single-byte type.
  typedef unsigned char Byte;

  typedef void (*CallbackFn)(const char * info, void * callback_data,
			     unsigned int n_reads, unsigned long long other);

};
