def filter_fasta_file_any(ht, filename, total_reads, outname, threshold):
    minmax = ht.fasta_file_to_minmax(filename, total_reads)
    readmask = ht.filter_fasta_file_any(filename, minmax, threshold)

    n_kept = readmask.filter_fasta_file(filename, outname)

    return total_reads, n_kept

def filter_fasta_file_all(ht, filename, total_reads, outname, threshold):
    minmax = ht.fasta_file_to_minmax(filename, total_reads)
    readmask = ht.filter_fasta_file_all(filename, minmax, threshold)

    n_kept = readmask.filter_fasta_file(filename, outname)

    return total_reads, n_kept
