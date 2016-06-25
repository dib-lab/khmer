#include "hashtable.hh"

namespace khmer {
    class CompositeGraph : public Hashtable {
    protected:
        Hashtable * _primary;
    public:
        explicit CompositeGraph(Hashtable * p) : Hashtable(p->ksize()), _primary(p) { }
        void count(const char * kmer) { }; // do nothing
        void count(HashIntoType khash) { }; // do nothing

        virtual const BoundedCounterType get_count(const char * kmer) const;
        virtual const BoundedCounterType get_count(HashIntoType khash) const;
        
        void save(std::string) { };
        void load(std::string) { };
        const HashIntoType n_unique_kmers() const { return 0; }
        const HashIntoType n_occupied() const { return 0; }
        void save_tagset(std::string) { }
        void load_tagset(std::string, bool clear_tags=true) { };

        BoundedCounterType test_and_set_bits(const char * kmer) { };
        BoundedCounterType test_and_set_bits(HashIntoType khash) { };

        std::vector<HashIntoType> get_tablesizes() const;
        const size_t n_tables() const;

        void print_tagset(std::string) { };
        void print_stop_tags(std::string) { };
        void save_stop_tags(std::string) { };
    };
}
