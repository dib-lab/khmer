import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys


def help():
    print("Usage: python3 plotGoldHist.py <outputPrefix>\n"
          "Plot a histogram for the true count\n"
          "Image will be saved on <outputPrefix>.Hist.png"
    )

if len(sys.argv)!=2 or  sys.argv[1] in ['-h','--h', '--help' ]:
    help()
    exit()


outputPrefix=sys.argv[1]
inputFile=open(outputPrefix+".gold")


counts=[int(l.split("\t")[1]) for l in inputFile]

n, bins, patches = plt.hist(counts, 50, normed=1, facecolor='g', alpha=0.75)
plt.xlabel('Kmers count')
plt.ylabel('Frequency')
plt.title("Kmers Frequency Distribution")
plt.savefig(outputPrefix+".goldHist.png")
