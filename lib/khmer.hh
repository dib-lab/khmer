#define VERSION "0.3"

#define MAX_COUNT 255
#define TAG_DENSITY 90
#define CONNECTED_THRESHOLD 0

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

  typedef void (*CallbackFn)(const char * info, void * callback_data,
			     unsigned int n_reads, unsigned long long other);
};
