import sys,math
from khmer import QFCounttable, Counttable

def help():
    print("Usage : python3 testloadFactor.py <dataset> <SketchSize>")
    print("   <dataset> contains unique kmers to test the load factor. No Kmers in the dataset must be slightly more than the no slots ")
    print("   <SketchSize> No of slots")
    print("Results will be printed in case of sucess to the standard output")

if len(sys.argv)!=3 or  sys.argv[1] in ['-h','--h', '--help' ]:
    help()
    exit(0)


def noSlotsIncrease(nslots):
    nslots+=10*math.sqrt(nslots)
    nblocks=int(float(nslots+64-1)/(64.0))
    return nblocks*64



dataset=open(sys.argv[1])
size=int(sys.argv[2]) #cqf
counterName='Unknown'
kmers=[x.strip() for x in dataset.readlines()]
noSlots=size
trueSize=noSlotsIncrease(noSlots)
print("True size",trueSize)
start=int(len(kmers)*0.80)

step=int(len(kmers)*0.01)

print("CQF slots=%d and size=%.2fM"%(noSlots, ((noSlots)*1.3)/(1000000)))
for i in range(start,len(kmers),step) :
    counter=QFCounttable(20,noSlots)
    counterName='CQF'
    print("Testing %i kmers load factor= %d"%(i,(float(i)/float(trueSize))*100.0))
    try:
        
        for kmer in kmers[:i]:
            counter.count(kmer)
    except:
        print(size,"FAILED")

    collisions=0
    for kmer in kmers[:i]:
        c=counter.get(kmer)
        if(c>1):
            collisions+=c
    newLoad=(float(i-collisions)/float(trueSize))*100
    print("Success and No Collisions=%d which makes the load factor=%d"%(collisions,newLoad))

            
    
