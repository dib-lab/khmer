import sys,math
from khmer import QFCounttable, Counttable, Nodetable


def help():
    print("Usage : python3 testSketchesAccuracy.py <dataset> <unseen dataset>")
    print("   <dataset> contains unique kmers to test the accuracy.")
    print("   <unseen dataset>")
    print("Results will be printed in case of sucess to the standard output")




def bloomEquivleantSize(nslots):
    nslots+=10*math.sqrt(nslots)
    nblocks=int(float(nslots+64-1)/(64.0))
    nslots=nblocks*64
    return int(nslots*8 + 2.125*nslots)

if len(sys.argv)!=3 or  sys.argv[1] in ['-h','--h', '--help' ]:
    help()
    exit(0)

dataset=open(sys.argv[1])
testingDataset=open(sys.argv[2])
kmers=[x.strip() for x in dataset.readlines()]
NotSeenKmers=[x.strip() for x in testingDataset.readlines()]



cqfSize=2**int(math.log2(len(kmers)))

accuracy=0.1
k=20
bitsPerElement=-1.44*math.log2(accuracy)
bloomSize=int(bloomEquivleantSize(cqfSize)/7)



def createQF(size):
    res=QFCounttable(k,size)
    return res


def createBloom(size):
    res=Nodetable(k,size,7)
    return res

counterQF=createQF(cqfSize)
counterBloom=createBloom(bloomSize)

print("No kmers\tCQF FP\tBloom FP")

countedKmers=0
for kmer in kmers:
    counterQF.count(kmer)
    counterBloom.count(kmer)
    countedKmers+=1
    completionPercent=float(countedKmers)/float(len(kmers))*100
    if completionPercent>20 and completionPercent%2==0:
        qfError=0
        bloomError=0
        for newkmer in NotSeenKmers:
            if counterQF.get(newkmer):
                qfError+=1
            if counterBloom.get(newkmer):
                bloomError+=1
        print("%d\t%d\t%d"%(countedKmers,(qfError),(bloomError)))
        sys.stdout.flush()
