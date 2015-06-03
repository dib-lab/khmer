#include <iostream>

#include <seqan/sequence.h>
#include <seqan/misc/name_store_cache.h>

using namespace seqan;

int main()
{
    typedef StringSet<CharString>      TNameStore;
    typedef NameStoreCache<TNameStore> TNameStoreCache;

    // Used below.
    unsigned idx = 0;
    bool success = false;

    // Define the name store and cache.
    TNameStore      nameStore;

    // Append two names to the name store without the cache.
    appendValue(nameStore, "I");
    appendValue(nameStore, "II");

    // Construct the name store.  The name store knows about both "I" and "II".
    TNameStoreCache nameStoreCache(nameStore);
    success = getIdByName(idx, nameStoreCache, "I");
    std::cout << "lookup I   => success=" << success << ", idx=" << idx << "\n";
    success = getIdByName(idx, nameStoreCache, "II");
    std::cout << "lookup II  => success=" << success << ", idx=" << idx << "\n";

    // Append value through cache.
    appendName(nameStoreCache, "III");
    success = getIdByName(idx, nameStoreCache, "III");
    std::cout << "lookup III => success=" << success << ", idx=" << idx << "\n";

    // Append value and query at the same time.
    idx = nameToId(nameStoreCache, "IV");
    std::cout << "lookup IV  => idx=" << idx << ", length(nameStore)="
              << length(nameStore) << "\n";

    return 0;
}
