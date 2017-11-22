import sys
from itertools import product
import random
import numpy as np

# N number of uniq kmers
# k is kmer length
# M number of kmers in the nonexisted dataset
# s is the zipifan dataset coeffieint
def GenerateKmers(N,k,outPrefix,M,s=1.5):
    Gold=open(outPrefix+".gold",'w')
    dataset=open(outPrefix+".dat",'w')
    NonExisted=open(outPrefix+".none.dat",'w')
    skipMax=random.randint(1,10000)
    skip=0
    rank=1
    i=0
    j=0
    s= np.random.zipf(s,N)
    for candidate in product('ACGT',repeat=k):
        skip+=1
        if skip <skipMax:
            continue
        skipMax=random.randint(1,1000)
        skip=0
        candidate="".join(candidate)
        if i<N:
            numberOfOccurences=min([s[i],1000])
            Gold.write("%s\t%d\n"%(candidate,numberOfOccurences))
            for j in range(0,int(numberOfOccurences)):
                dataset.write("%s\n"%candidate)
            i+=1
        elif j<M:
            NonExisted.write("%s\n"%candidate)
            j+=1
        else:
            break
    Gold.close()
    dataset.close()
    NonExisted.close()
        



    
def help():
    print("python3 generateSeq.py <NumberOfUniqeKmers> <K> <NumberOfUnseenKmers> <outputPrefix> ")
    print("Script generate 3 files:\n"
          "1) Gold: Kmers with the true count\n"
          "2) Dataset: Kmers to be  counted\n"
          "3) Unseen: kmers that are not in the dataset"
    )
    
if __name__ =='__main__':
    if sys.argv[1] in ['-h','--h', '--help' ]:
        help()
        exit(0)
    NoUniqueKmers=int(sys.argv[1])
    K=int(sys.argv[2])
    NumberOfUnseenKmers=int(sys.argv[3])
    outPrefix=sys.argv[4]
    GenerateKmers(NoUniqueKmers,K,outPrefix,NumberOfUnseenKmers)

