import khmer
import khmer_tst_utils as utils
import os.path

def test_bin_index_and_retrieve():
    data_filename = utils.get_test_data('test-outline-simple.fa')

    # index the reads
    bin_filename = khmer.convert_fasta_to_indexed_bin(data_filename)
    print 'bin_filename is', bin_filename
    assert os.path.exists(bin_filename)

    # retrieve read
    read = khmer.outline_retrieve_read_by_id(bin_filename, 1)
    print 'READ IS', read
    assert read == \
           'CGCAGGCTGGATTCTAGAGGCAGAGGTGAGCTATAAGATATTGCATACGTTGAGCCAGC', \
           read


def test_simple_retrieval():
    data_filename = utils.get_test_data('test-outline-simple.fa')

    # index the reads
    bin_filename = khmer.convert_fasta_to_indexed_bin(data_filename)

    # this tag is the first 32 bases of the only sequence in
    # 'test-outline-simple.fa'
    tag = 'CGCAGGCTGGATTCTAGAGGCAGAGGTGAGCT'
    tags = [khmer.forward_hash(tag, 32)]

    ht = khmer.new_hashbits(32, 1, 1)
    ht.add_tag(tag)

    # build the index connect tags <=> reads
    index_filename = ht.build_outline_index(bin_filename)

    # OK, done building indices, etc.

    # retrieve list of IDs
    read_id_list = khmer.outline_retrieve_read_ids_by_taglist(index_filename,
                                                              tags)

    # should only be one sequence!
    assert len(read_id_list) == 1

    # retrieve read
    read = khmer.outline_retrieve_read_by_id(bin_filename, read_id_list[0])
    print 'READ IS', (read,), read_id_list
    assert read == \
           'CGCAGGCTGGATTCTAGAGGCAGAGGTGAGCTATAAGATATTGCATACGTTGAGCCAGC', \
           read

def test_2seq_retrieval():
    data_filename = utils.get_test_data('test-outline-simple2.fa')

    # index the reads
    bin_filename = khmer.convert_fasta_to_indexed_bin(data_filename)

    # this tag is the first 32 bases of the first sequence in
    # 'test-outline-simple2.fa', and the last 32 bases of the second
    # sequence.
    tag = 'CGCAGGCTGGATTCTAGAGGCAGAGGTGAGCT'
    tags = [khmer.forward_hash(tag, 32)]

    ht = khmer.new_hashbits(32, 1, 1)
    ht.add_tag(tag)

    # build the index connect tags <=> reads
    index_filename = ht.build_outline_index(bin_filename)

    print data_filename, bin_filename, index_filename

    # OK, done building indices, etc.

    # retrieve list of IDs
    read_id_list = khmer.outline_retrieve_read_ids_by_taglist(index_filename,
                                                              tags)

    # should be two sequences
    assert len(read_id_list) == 2

    # retrieve read 1
    read = khmer.outline_retrieve_read_by_id(bin_filename, read_id_list[0])
    print 'READ IS', (read,), read_id_list
    assert read == \
           'CGCAGGCTGGATTCTAGAGGCAGAGGTGAGCTATAAGATATTGCATACGTTGAGCCAGC', \
           read

    # retrieve read 2
    read = khmer.outline_retrieve_read_by_id(bin_filename, read_id_list[1])
    print 'READ IS', (read,), read_id_list
    assert read == \
           'TAGAGCCGATGAGATGCAGAGTAGAGACGCAGGCTGGATTCTAGAGGCAGAGGTGAGCT', \
           read

def test_kmer_traversal_retrieval():
    data_filename = utils.get_test_data('test-outline-simple.fa')

    # index the reads
    bin_filename = khmer.convert_fasta_to_indexed_bin(data_filename)

    # this tag is the first 32 bases of the first sequence in
    # 'test-outline-simple.fa'.

    ht = khmer.new_hashbits(32, 1000, 2)
    tag = 'CGCAGGCTGGATTCTAGAGGCAGAGGTGAGCT'
    ht.add_tag(tag)

    # build the graph
    ht.consume_fasta(data_filename)

    # build the index connect tags <=> reads
    index_filename = ht.build_outline_index(bin_filename)

    # OK, done building indices, etc.

    # do find_all_tags
    kmer = 'GAGCTATAAGATATTGCATACGTTGAGCCAGC'
    
    # @CTB should just use find_all_tags...
    tags = ht.find_all_tags_to_taglist(kmer)

    print tags

    # retrieve list of IDs
    read_id_list = khmer.outline_retrieve_read_ids_by_taglist(index_filename,
                                                              tags)

    print read_id_list

    # should be two sequences
    assert len(read_id_list) == 1

    # retrieve read 1
    read = khmer.outline_retrieve_read_by_id(bin_filename, read_id_list[0])
    print 'READ IS', (read,), read_id_list
    assert read == \
           'CGCAGGCTGGATTCTAGAGGCAGAGGTGAGCTATAAGATATTGCATACGTTGAGCCAGC', \
           read
