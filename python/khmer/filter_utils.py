def filter_fasta_file(ht, filename, total_reads, outname, threshold):
    minmax = ht.fasta_file_to_minmax(filename, total_reads)
    readmask = ht.filter_fasta_file_max(filename, minmax, threshold)

    n_kept = ht.output_filtered_fasta_file(filename, outname, readmask)

    return total_reads, n_kept
