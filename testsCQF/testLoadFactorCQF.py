import sys,math
from khmer import QFCounttable, Counttable

def help():
    print("Usage : python3 testloadFactorCQF.py <dataset> <SketchSize> <No Repeat>")
    print("   <dataset> contains unique kmers to test the load factor. No Kmers in the dataset must be slightly more than the no slots ")
    print("   <SketchSize> No of slots")
    print("   <No Repeat> each kmer will repeated X times")
    print("Results will be printed in case of sucess to the standard output")



def noSlotsIncrease(nslots):
    nslots+=10*math.sqrt(nslots)
    nblocks=int(float(nslots+64-1)/(64.0))
    return nblocks*64



if len(sys.argv)!=4 or  sys.argv[1] in ['-h','--h', '--help' ]:
    help()
    exit(0)

dataset=open(sys.argv[1])
size=int(sys.argv[2]) #cqf
nrepeat=int(sys.argv[3])
kmers=[x.strip() for x in dataset.readlines()]
noSlots=size
trueSize=noSlotsIncrease(noSlots)

counter=QFCounttable(20,noSlots)

countedKmers=0
for kmer in kmers:
    if counter.get(kmer)==0:
        for j in range(nrepeat):
            counter.count(kmer)
        countedKmers+=1
        print("%d\t%d"%(nrepeat,countedKmers))
        sys.stdout.flush()
            
    
