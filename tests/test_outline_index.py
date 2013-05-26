import khmer
import khmer_tst_utils as utils

def test_simple_retrieval():
    data_filename = utils.get_test_data('test-outline-simple.fa')

    # index the reads
    bin_filename = khmer.convert_fasta_to_indexed_bin(data_filename)

    # this tag is the first 32 bases of the only sequence in
    # 'test-outline-simple.fa'
    tag = 'CGCAGGCTGGATTCTAGAGGCAGAGGTGAGCT'
    tags = [tag]

    ht = khmer.new_hashbits(32, 1, 1)
    ht.add_tag(tags)

    # return the sorted/numerical tags (well, tag, in this case).
    # @CTB could we do this with 'sort' if they weere
    sorted_tags = ht.get_sorted_tags(sorted_tags)

    # build the index connect tags <=> reads
    index_filename = khmer.build_outline_index(data_filename, sorted_tags)

    # OK, done building indices, etc.

    # retrieve list of IDs
    read_id_list = khmer.outline_retrieve_read_ids_by_tag(index_filename,
                                                          tags)

    # should only be one sequence!
    assert len(read_id_list) == 1

    # retrieve read
    read = khmer.outline_retrieve_read_by_id(bin_filename, read_id_list[0])
    assert read == \
           'CGCAGGCTGGATTCTAGAGGCAGAGGTGAGCTATAAGATATTGCATACGTTGAGCCAGC',
           read
