import sys
from screed.fasta import fasta_iter

MAX_SIZE=5000

def read_partition_file(fp):
    for n, line in enumerate(fp):
        if n % 2 == 0:
            surrendered = False
            name, partition_id = line[1:].strip().rsplit('\t', 1)

            if '*' in partition_id:
                partition_id = int(partition_id[:-1])
                surrendered = True
            else:
                partition_id = int(partition_id)
        else:
            sequence = line.strip()

            yield name, partition_id, surrendered, sequence


select_pid = int(sys.argv[2])
count = 0
for n, (name, pid, _, seq) in enumerate(read_partition_file(open(sys.argv[1]))):
    if pid == select_pid:
        print '>%s\t%d\n%s' % (name, pid, seq)
        count += 1

    if n % 10000 == 0:
        sys.stderr.write('...%d\n' % (n,))

sys.stderr.write('found %d total in partition %d\n' % (count, pid))
