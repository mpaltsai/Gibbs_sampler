#setwd('/home/ismini/Documents/bioinfo-grad/Algorithms_Bioinformatics/nikolaou/2nd_assingment/')

library(ggplot2)
library(ggseqlogo)

data1 = read.csv('array.csv', sep=',', header=F)

data2 = as.matrix(data1)
data2[data2<0] <- 0
dimnames(data2) = list(c('A', 'C', 'G', 'T'))
ggseqlogo(data2, method='custom', seq_type='dna')

png('ic_seqlogo.png')
ggseqlogo(data2, method='custom', seq_type='dna')
dev.off()

