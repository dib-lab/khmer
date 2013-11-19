import khmer
import sys
import screed
import graph_tool.all as gt


input_fasta = sys.argv[3]
K = sys.argv[1]
x = sys.argv[2]


ht = khmer.new_hashbits(K, x, 4)

sparse_graph = gt.Graph()

for n, record in enumerate(screed.open(input_fasta)):
    
