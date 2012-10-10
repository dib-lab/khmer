import math

for p in [0.001, 0.01, 0.05, 0.1, 0.15, 0.20]:
   print p, 1 / math.log(2, 2) * math.log(1 / p, 2)
