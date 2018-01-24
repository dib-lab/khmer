import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys


def help():
    print("Usage: python3 plotLoadingFactor.py results.tsv (normal|log) outputImg.jpg\n"
          "results.tsv produced using commands at run.sh under Load factor test\n"
          "normal|log describes the type of scale of the xaxis\n"
          "Outputs one images: \n"
          "     1) Line Graph for the maximum loading \n"
    )

if len(sys.argv)!=4 or  sys.argv[1] in ['-h','--h', '--help' ]:
    help()
    exit()


inputTsv=open(sys.argv[1],'r')
outImg=sys.argv[3]
if sys.argv[2]=='log':
    logScale=True
elif sys.argv[2]=='normal':
    logScale=False
else:
    help()
    exit()




inputTsv=[l.strip().split("\t") for l in inputTsv]
inputTsv=inputTsv[1:]
x=np.array([int(l[0]) for l in inputTsv])
y=np.array([int(l[1]) for l in inputTsv])

fig, ax = plt.subplots()

fig.set_size_inches(9, 6)
if logScale:
    line1, = ax.loglog(x, y,basex=2, linewidth=2,
                 label='CQF')
else:
    line1, = ax.plot(x, y, linewidth=2,
                 label='CQF')
ax.legend(loc='upper right')
plt.xlabel('M')
plt.ylabel('Maximum Number of Unique Kmers')
plt.title('Maximum Loading')
plt.grid(True)
for i,j in zip(x,y):
    ax.annotate(str(j),xy=(i,j))
plt.savefig(outImg)
plt.close()
