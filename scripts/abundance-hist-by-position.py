import sys

first = True
n = 0
countSum = []

fd = open(sys.argv[1])
for n, line in enumerate(fd):
   if n % 100000 == 0:
      print >>sys.stderr, '...', n

   tok = line.split()

   if first:
      countSum = [0]*len(tok)
      first = False

   for i in range(len(tok)):
      countSum[i] += int(tok[i])

y = [0.0]*len(countSum)

for i in range(len(countSum)):
   y[i] = float(countSum[i]) / n

for n, i in enumerate(y):
   print n, i
