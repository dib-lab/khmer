namespace khmer
{

template <typename HashFunctor>
class BaseFilter
{
protected:
    HashFunctor hash_function;

    explicit BaseFilter(HashFunctor hash_function)
        : hash_function(hash_function)
    {

    }

    explicit BaseFilter(const BaseFilter&);
    BaseFilter& operator=(const BaseFilter&);

    virtual ~BaseFilter( ) {}

    virtual void _allocate_counters();

public:

    virtual void count(const char * kmer) = 0;
    virtual void count(Kmer khash) = 0;

    // get the count for the given k-mer.
    virtual const BoundedCounterType get_count(const char * kmer) const = 0;
    virtual const BoundedCounterType get_count(Kmer khash) const = 0;

    virtual void save(std::string) = 0;
    virtual void load(std::string) = 0;

    // number of unique k-mers
    virtual const HashIntoType n_unique_kmers() const = 0;

    // count number of occupied bins
    virtual const HashIntoType n_occupied() const = 0;

};

}
