//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#include <iostream>
#include "hashtable.hh"

using namespace std;

int main()
{
#if 0
    khmer::Hashtable ht(5, 1024);

    ht.count("AAAAA");
    ht.count("AAAAT");
    ht.count("AAATG");
    ht.count("AATGG");

    ht.count("CCCCC");

    std::string k = "AATGG";
    ht.mark_connected_graph(k);

    std::cout << "a\n";

    ht.dump_kmers_and_counts();

    std::cout << "b\n";

    ht.empty_bins(false);
    ht.dump_kmers_and_counts();

    std::cout << "c\n";

    ht.count("CCCCC");
    ht.empty_bins(true);
    ht.dump_kmers_and_counts();

    std::cout << "d\n";
    ht.count("AAAAA");
    ht.count("AAAAT");
    ht.count("AAATG");
    ht.count("AATGG");
    ht.count("CCCCC");

    std::cout << ht.calc_connected_graph_size("AAAAA") << "\n";
    std::cout << ht.calc_connected_graph_size("AAAAA") << "\n";
    ht.clear_marks();
    std::cout << ht.calc_connected_graph_size("AAAAA") << "\n";
    ht.clear_marks();
    std::cout << ht.calc_connected_graph_size("AAAAT") << "\n";
    ht.clear_marks();
    std::cout << ht.calc_connected_graph_size("AAATG") << "\n";
    ht.clear_marks();
    std::cout << ht.calc_connected_graph_size("AATGG") << "\n";
    ht.clear_marks();
    std::cout << ht.calc_connected_graph_size("CCCCC") << "\n";

    std::cout << "e\n";
    ht.clear_marks();
    std::cout << ht.calc_connected_graph_size("AAAAT") << "\n";
    std::cout << ht.calc_connected_graph_size("CCCCC") << "\n";
    ht.clear_marks();
    ht.trim_graphs(2);

#endif

    std::string filename = "../foo.fa";
    unsigned int total_reads;
    unsigned long long n_consumed;
    khmer::Hashtable ht2(16, 4294967296);
    ht2.consume_fasta(filename, total_reads, n_consumed);
    ht2.trim_graphs(100);

    return 0;
}
