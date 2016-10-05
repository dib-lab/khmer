#include "hashtable.hh"

namespace khmer {
    class CompositeGraph : public Hashtable {
    public:
        explicit CompositeGraph(WordLength ksize) : Hashtable(ksize) { }
        virtual const BoundedCounterType get_count(const char * kmer)
            const {
            HashIntoType hash = _hash(kmer, _ksize);
            return get_count(hash);
        }
        virtual const BoundedCounterType get_count(HashIntoType khash)
            const = 0;

        // these graphs are read-only; all below functions =/=.
        
        void count(const char * kmer) { }; // do nothing
        void count(HashIntoType khash) { }; // do nothing

        void save(std::string) { };
        void load(std::string) { };
        const HashIntoType n_unique_kmers() const { return 0; }
        const HashIntoType n_occupied() const { return 0; }
        void save_tagset(std::string) { }
        void load_tagset(std::string, bool clear_tags=true) { };

        BoundedCounterType test_and_set_bits(const char * kmer) { };
        BoundedCounterType test_and_set_bits(HashIntoType khash) { };

        std::vector<HashIntoType> get_tablesizes() const {
            std::vector<HashIntoType> ts;
            return ts;
        }
        const size_t n_tables() const {
            return 0;
        }

        void print_tagset(std::string) { };
        void print_stop_tags(std::string) { };
        void save_stop_tags(std::string) { };
    };

    class CompositeGraph_AndNot : public CompositeGraph {
    protected:
        Hashtable * primary, * secondary;
    public:
        explicit CompositeGraph_AndNot(Hashtable *p, Hashtable *s) :
            CompositeGraph(p->ksize()), primary(p), secondary(s) { };

        const BoundedCounterType get_count(HashIntoType khash) const {
            if (!secondary->get_count(khash)) {
                return primary->get_count(khash);
            }
            return 0;
        }
    };
}
