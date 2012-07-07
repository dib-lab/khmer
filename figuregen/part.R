myargs <- commandArgs(trailingOnly=T)
filename <- myargs[1]
filenametwo <- myargs[2]
data = read.csv(filename, header=FALSE)

plot(data$V1, data$V2, xlab="False Positive Rate", xlim=c(0, 0.16), ylim=c(0, 1010), ylab="Number of Partitions", col="black")

if (file.exists(filenametwo)) {

   datatwo = read.csv(filenametwo, header=FALSE)

   arrows(datatwo$V1, datatwo$V4, datatwo$V1, datatwo$V2, col="black", angle=90, code=3, length=0.04)

   points(datatwo$V1, datatwo$V3, pch=20, col="black")

}

