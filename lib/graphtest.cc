#include <iostream>
#include "hashtable.hh"

using namespace std;

int main()
{
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

  return 0;
}
