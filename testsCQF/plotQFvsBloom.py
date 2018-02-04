import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys


def help():
    print("Usage: python3 plotQFvsBloom.py results.tsv outputImg.jpg\n"
          "results.tsv produced using commands at run.sh under Load factor test\n"
          "Outputs one image: \n"
          "     1) Line Graph for bloom filter and CQF accuracy \n"
    )

if len(sys.argv)!=3 or  sys.argv[1] in ['-h','--h', '--help' ]:
    help()
    exit()


inputTsv=open(sys.argv[1],'r')
outImg=sys.argv[2]



inputTsv=[l.strip().split("\t") for l in inputTsv]
inputTsv=inputTsv[1:]
x=np.array([int(l[0]) for l in inputTsv])
cqf=np.array([int(l[1]) for l in inputTsv])
bloom=np.array([int(l[2]) for l in inputTsv])

fig, ax = plt.subplots()

fig.set_size_inches(9, 6)

line1, = ax.plot(x,cqf,color='blue', linewidth=2,label='CQF')
line2, = ax.plot(x,bloom, color='red',linewidth=2,label='Bloom Filter')
ax.legend(loc='lower right')
plt.xlabel('No Kmers')
plt.ylabel('False Positives')
plt.title('Quotient Filter and Bloom Filter of size~(524KB)')
plt.grid(True)
plt.savefig(outImg)
plt.close()
