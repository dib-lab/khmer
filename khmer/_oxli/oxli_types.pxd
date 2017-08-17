from libcpp cimport bool
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.string cimport string


cdef extern from "oxli/oxli.hh":
    cdef int MAX_BIGCOUNT

cdef extern from "oxli/oxli.hh" namespace "oxli":
    ctypedef unsigned long long int HashIntoType
    ctypedef set[HashIntoType] HashIntoTypeSet

    ctypedef unsigned int PartitionID
    ctypedef set[PartitionID] PartitionSet
    ctypedef map[HashIntoType, PartitionID*] PartitionMap
    ctypedef map[PartitionID, PartitionID*] PartitionPtrMap
    ctypedef map[PartitionID, HashIntoTypeSet*] PartitionToTagsMap
    ctypedef map[PartitionID, unsigned int] PartitionCountMap
    ctypedef map[unsigned long long, unsigned long long] PartitionCountDistribution
    ctypedef map[HashIntoType, unsigned int] TagCountMap
    ctypedef unsigned char WordLength
    ctypedef unsigned short int BoundedCounterType

    ctypedef unsigned long long int Label
    ctypedef set[Label] LabelSet
    ctypedef set[HashIntoType] TagSet

    ctypedef void (*CallbackFn)(const char *, void *, uint64_t, uint64_t)
