myargs <- commandArgs(trailingOnly=T)
filename <- myargs[1]
data = read.csv(filename, header=TRUE)

plot(data$FPR, data$AVG, xlab="False Positive Rate", xlim=c(0, 0.31), ylim=c(0, 26), ylab="Length of Longest Shortest Path", col="black", pch=20)

arrows(data$FPR, data$UPPER, data$FPR, data$LOWER, angle=90, code=3, length=0.04)

arrows(0.175, 21, 0.185, 24.5)
