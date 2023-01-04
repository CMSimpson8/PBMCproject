library(ggplot2)
library(ShortRead)

#uploads the first FASTQ file into the environment
FQ_R1R1 = readFastq("/home/genomics/fastq_data/Post-R1_R1.fastq", header = TRUE, StringAsFactor = FALSE, row.names=2())

#Summary give information about e length of the DNA sequence
#summary(FQ_R1R1)
#head gives the number of reads
#head(FQ_R1R1)
# want to look at the DNA sequences 
reads =  sread(FQ_R1R1)
head(reads)
widths = as.data.frame(reads@ranges@width)

ggplot(widths) +
  geom_histogram(aes(x=reads@ranges@width), bins = 5)
#we want to graph the quality

quals = quality(FQ_R1R1)
head(quals)
numscores  = as(quals, 'matrix')
avgscore = rowMeans(numscores)