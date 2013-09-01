#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
myargs <- commandArgs(trailingOnly=T)
filename <- myargs[1]
data = read.csv(filename, header=TRUE)

plot(data$FPR, data$AVG, xlab="False Positive Rate", xlim=c(0,0.175), ylim=c(0, 10), ylab="Average Component Size", col="black", pch=20)

#plot(1, 1, xlim=c(0, 0.18), ylim=c(0, 20), xlab="False Positive Rate", ylab="Cluster Size", col="black", pch=20)

arrows(data$FPR, data$UPPER, data$FPR, data$LOWER, angle=90, code=3, length=0.04)

warnings()
