#include <iostream>
#include <chrono>
#include <ctime>
#include <memory>
#include <vector>
#include <cstdlib>

#include <oxli/kmer_hash.hh>
#include <oxli/read_parsers.hh>
#include <oxli/gmap.hh>

using namespace oxli;
using namespace oxli::read_parsers;

#define K 21

unsigned long long llrand() {
    unsigned long long r = 0;

    for (int i = 0; i < 5; ++i) {
        r = (r << 15) | (rand() & 0x7FFF);
    }

    return r & 0xFFFFFFFFFFFFFFFFULL;
}

FastxParserPtr get_test_reads() {
    FastxParserPtr parser = get_parser<FastxReader>("../../tests/test-data/test-reads.fa");
    return parser;
}

vector<HashIntoType> * get_test_kmers(int num_hashes=5000000) {
    vector<HashIntoType> * hashes = new vector<HashIntoType>();
    while(num_hashes > 0) {
        hashes->push_back(llrand());
        num_hashes--;
    }
    return hashes;
}


void fill_gmap(GuardedHashMap<int, false>& _map, vector<HashIntoType> * hashes) {
    for(auto hash: *hashes) {
        _map.set(hash, rand());
    }
}


void fill_uomap(std::unordered_map<HashIntoType, int>& _map, vector<HashIntoType> * hashes) {
    for (auto hash: *hashes) {
        _map[hash] = rand();
    }
}

void fill_map(std::map<HashIntoType, int>& _map, vector<HashIntoType> * hashes) {
    for (auto hash: *hashes) {
        _map[hash] = rand();
    }
}

void test_gmap(vector<HashIntoType> * hashes) {

    std::cout << "=== GMAP ===" << std::endl;

    vector<double> get_full_times;
    vector<double> get_empty_times;
    vector<double> get_bad_times;
    GuardedHashMap<int, false> _map(K, 4, 1000000);
    std::chrono::time_point<std::chrono::system_clock> start, end;

    fill_gmap(_map, hashes);
    for (auto hash: *hashes) {
        start = std::chrono::system_clock::now();
        int result = _map.get(hash);
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        get_full_times.push_back(elapsed_seconds.count());
    }
    double avg_get_full_time = std::accumulate(get_full_times.begin(), 
                                               get_full_times.end(), 0.0) / get_full_times.size();
    std::cout << "Avg full get time: " << avg_get_full_time << std::endl;


    vector<HashIntoType> * newhashes = get_test_kmers();
    for (auto hash: *newhashes) {
        start = std::chrono::system_clock::now();
        int result = _map.get(hash);
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        get_bad_times.push_back(elapsed_seconds.count());
    }
    double avg_get_bad_time = std::accumulate(get_bad_times.begin(), 
                                               get_bad_times.end(), 0.0) / get_bad_times.size();
    std::cout << "Avg bad get time: " << avg_get_bad_time << std::endl;
    delete newhashes;


    _map = GuardedHashMap<int, false>(K, 4, 1000000);

    for (auto hash: *hashes) {
        start = std::chrono::system_clock::now();
        int result = _map.get(hash);
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;

        get_empty_times.push_back(elapsed_seconds.count());
    }

    double avg_get_empty_time = std::accumulate(get_empty_times.begin(), 
                                               get_empty_times.end(), 0.0) / get_empty_times.size();
    std::cout << "Avg empty get time: " << avg_get_empty_time << std::endl;
}

void test_uomap(vector<HashIntoType> * hashes) {

    std::cout << "=== UOMAP ===" << std::endl;

    vector<double> get_full_times;
    vector<double> get_empty_times;
    vector<double> get_bad_times;
    std::unordered_map<HashIntoType, int> _map;
    std::chrono::time_point<std::chrono::system_clock> start, end;

    fill_uomap(_map, hashes);

    for (auto hash: *hashes) {
        start = std::chrono::system_clock::now();
        int result;
        auto search = _map.find(hash);
        if (search != _map.end()) {
            result = search->second;
        }
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;

        get_full_times.push_back(elapsed_seconds.count());
    }
    double avg_get_full_time = std::accumulate(get_full_times.begin(), 
                                               get_full_times.end(), 0.0) / get_full_times.size();
    std::cout << "Avg full get time: " << avg_get_full_time << std::endl;


    vector<HashIntoType> * newhashes = get_test_kmers();
    for (auto hash: *newhashes) {
        start = std::chrono::system_clock::now();
        int result;
        auto search = _map.find(hash);
        if (search != _map.end()) {
            result = search->second;
        }
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        get_bad_times.push_back(elapsed_seconds.count());
    }
    double avg_get_bad_time = std::accumulate(get_bad_times.begin(), 
                                               get_bad_times.end(), 0.0) / get_bad_times.size();
    std::cout << "Avg bad get time: " << avg_get_bad_time << std::endl;
    delete newhashes;


    _map = std::unordered_map<HashIntoType, int>();
    for (auto hash: *hashes) {
        start = std::chrono::system_clock::now();
        int result;
        auto search = _map.find(hash);
        if (search != _map.end()) {
            result = search->second;
        }
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;

        get_empty_times.push_back(elapsed_seconds.count());
    }

    double avg_get_empty_time = std::accumulate(get_empty_times.begin(), 
                                               get_empty_times.end(), 0.0) / get_empty_times.size();
    std::cout << "Avg empty get time: " << avg_get_empty_time << std::endl;
}

void test_map(vector<HashIntoType> * hashes) {

    std::cout << "=== MAP ===" << std::endl;

    vector<double> get_full_times;
    vector<double> get_empty_times;
    vector<double> get_bad_times;
    std::map<HashIntoType, int> _map;
    std::chrono::time_point<std::chrono::system_clock> start, end;

    fill_map(_map, hashes);
    for (auto hash: *hashes) {
        start = std::chrono::system_clock::now();
        int result;
        auto search = _map.find(hash);
        if (search != _map.end()) {
            result = search->second;
        }
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;

        get_full_times.push_back(elapsed_seconds.count());
    }
    double avg_get_full_time = std::accumulate(get_full_times.begin(), 
                                               get_full_times.end(), 0.0) / get_full_times.size();
    std::cout << "Avg full get time: " << avg_get_full_time << std::endl;

    vector<HashIntoType> * newhashes = get_test_kmers();
    for (auto hash: *newhashes) {
        start = std::chrono::system_clock::now();
        int result;
        auto search = _map.find(hash);
        if (search != _map.end()) {
            result = search->second;
        }
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        get_bad_times.push_back(elapsed_seconds.count());
    }
    double avg_get_bad_time = std::accumulate(get_bad_times.begin(), 
                                               get_bad_times.end(), 0.0) / get_bad_times.size();
    std::cout << "Avg bad get time: " << avg_get_bad_time << std::endl;
    delete newhashes;

    _map = std::map<HashIntoType, int>();
    for (auto hash: *hashes) {
        start = std::chrono::system_clock::now();
        int result;
        auto search = _map.find(hash);
        if (search != _map.end()) {
            result = search->second;
        }
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;

        get_empty_times.push_back(elapsed_seconds.count());
    }

    double avg_get_empty_time = std::accumulate(get_empty_times.begin(), 
                                               get_empty_times.end(), 0.0) / get_empty_times.size();
    std::cout << "Avg empty get time: " << avg_get_empty_time << std::endl;
}


int main() {
    vector<HashIntoType> * hashes = get_test_kmers();
    test_gmap(hashes);
    test_uomap(hashes);
    test_map(hashes);
}
