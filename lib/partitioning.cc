
void StreamingPartitioner::consume_sequence(const std::string& seq,
                                            uint64_t& n_consumed)
{
    
    KmerIterator kmers(seq.c_str(), graph->ksize());
    HashIntoType kmer;
    bool kmer_tagged;
    unsigned int since = graph->_tag_density / 2 + 1;

    std::set<HashIntoType> tags;
    bool known_component = false;

    // First check if we overlap any tags
    kmer = kmers.next();
    tags->insert(kmer); //always tag the first k-mer

    while(!kmers.done()) {
        kmer = kmers.next();
        bool is_new_kmer;
        is_new_kmer = graph->test_and_set_bits(kmer);

        if (is_new_kmer) {
            ++since;
        } else {
            kmer_tagged = tag_component_map.contains(kmer);
            if (kmer_tagged) {
                since = 1;
                tags->insert(kmer);
                known_component = true;
            } else {
                ++since;
            }
        }

        if (since >= graph->_tag_density) {
            tags.insert(kmer);
            since = 1;
        }
    }
    tags.insert(kmer);	// always tag the last k-mer

    if (known_component) {
        
    } else {

    }
}
