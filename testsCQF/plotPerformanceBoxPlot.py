import matplotlib.pyplot as plt
import numpy as np
import sys

def help():
    print("Usage: python3 plotPerformanceBoxPlot.py outputPrefix\n"
          "outputPrefix used through the script\n"
          "Outputs two images: \n"
          "     1) BoxPlots of the diff values between count in the sketch and true count from Gold\n"
          "     2) The same ase the first image, but New Kmers are used"
    )

if len(sys.argv)!=2 or  sys.argv[1] in ['-h','--h', '--help' ]:
    help()
    exit()



outputPrefix=sys.argv[1]

performanceFile=open(outputPrefix+'.result')
outImagePrefix=outputPrefix

results=[]

line=performanceFile.readline()
while line!='':
    tmp=[]
    counterName,size,err1,err2,NonExist=line.strip().split(" ")
    tmp.append(int(size))
    tmp.append(counterName)
    line=performanceFile.readline()
    tmp.append([int(s) for s in line.strip().split(',')])
    line=performanceFile.readline()
    tmp.append([int(s) for s in line.strip().split(',')])    

    results.append(tmp)

    line=performanceFile.readline().strip()

results.sort(key=lambda x:x[0])


plt.figure(figsize=(len(results)+3,4))
bp=plt.boxplot([r[2] for r in results],patch_artist=True)
boxnames=["%s%.2fM"%(x[1],x[0]/1000000.0) for x in results]

for i,box in enumerate(bp['boxes']):
    if results[i][1]=='CM':
        box.set(color='red', linewidth=2)
    else:
        box.set(color='blue', linewidth=2)
    box.set(hatch = '/')
    
plt.xticks(list(range(1,len(boxnames)+1)),boxnames)
plt.title("Kmers exists in the sketch")
plt.savefig(outImagePrefix+'.Exist.png')

plt.close()


plt.figure(figsize=(len(results)+3,4))
bp=plt.boxplot([r[3] for r in results],patch_artist=True)
boxnames=["%s%.2fM"%(x[1],x[0]/1000000.0) for x in results]

for i,box in enumerate(bp['boxes']):
    # change outline color
    if results[i][1]=='CM':
        box.set(color='red', linewidth=2)
    else:
        box.set(color='blue', linewidth=2)
    # change fill color
    # change hatch
    box.set(hatch = '/')
    
plt.xticks(list(range(1,len(boxnames)+1)),boxnames)
plt.title("Kmers doesnt exist in the sketch")
plt.savefig(outImagePrefix+'.NonExist.png')


