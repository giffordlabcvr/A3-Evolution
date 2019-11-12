#!/usr/bin/env R

###############Comments###############
##Description
#Calculate the strand bias of G-to-A mutation rates (i.e., mutation rate ratio between G-to-A and C-to-T in the positive strand)
#
##Usage
#R --vanilla --slave --args \
#  <family-level mutation info file> \
#  <output file name> \
#  < calc_mut_rate_family-level.R
#
##Output
#1, class: class
#2, family: family
#3, a.site.q: the number of A sites in the consensus sequence (total value among all alignments)
#4, t.site.q: the number of T sites in the consensus sequence (total value among all alignments)
#5, g.site.q: the number of G sites in the consensus sequence (total value among all alignments)
#6, c.site.q: the number of C sites in the consensus sequence (total value among all alignments)
#7, ac: the number of A-to-C mutations (total value among all alignments)
#8, ag: the number of A-to-G mutations (total value among all alignments)
#9, at: the number of A-to-T mutations (total value among all alignments)
#10, ca: the number of C-to-A mutations (total value among all alignments)
#11, cg: the number of C-to-G mutations (total value among all alignments)
#12, ct: the number of C-to-T mutations (total value among all alignments)
#13, ga: the number of G-to-A mutations (total value among all alignments)
#14, gc: the number of G-to-C mutations (total value among all alignments)
#15, gt: the number of G-to-T mutations (total value among all alignments)
#16, ta: the number of T-to-A mutations (total value among all alignments)
#17, tc: the number of T-to-C mutations (total value among all alignments)
#18, tg: the number of T-to-G mutations (total value among all alignments)
#19, Label: Label (basically equal to class, but the ambiguity is eliminated)
#20, total.mut: the total number of mutations among all alignments 
#21, g2a.rate: overall G-to-A mutation rate (across all alignments)
#22, c2t.rate: overall G-to-A mutation rate (across all alignments)
#23, g2a.c2t: mutation rate ratio between G-to-A and C-to-T in the positive strand (i.e., the strand bias of G-to-A mutation rates)
#24, g2a.c2t.log: log2-transformed mutation rate ratio
#25, pval: P value calculated by Fisher exact test
#26, FER: Family-wise error rate calculated by Bonferroni method
#
######################################


args <- commandArgs(trailingOnly = T)

fisherTest <- function(x){
  mt <- matrix(x, nrow=2, byrow=T)
  test.fisher <- fisher.test(mt,alternative="t")
  pval <- test.fisher$p.value
  return(pval)
}


data.name <- args[1]
out.name <- args[2]

data <- read.table(data.name,header=T)

Labels.v <- c()
for(i in 1:nrow(data)){
  Label <- 'Others'
  int3 <- substr(as.character(data$class[i]),1,3)
  if(int3 == 'LTR'){
    Label <- 'LTR'
  }
  if(int3 == 'LIN'){
    Label <- 'LINE'
  }
  if(int3 == 'SIN'){
    Label <- 'SINE'
  }
  if(int3 == 'DNA'){
    Label <- 'DNA'
  }
  Labels.v <- c(Labels.v,Label)
}

data$Label <- Labels.v

data$total.mut <- apply(data[,7:18],1,sum)

data <- data[data$Label != 'Others',]
data$Label <- factor(data$Label,levels=c('LTR','LINE','SINE','DNA'))

data <- data[data$total.mut > 1000,]



data$g2a.rate <- data$ga / data$g.site.q
data$c2t.rate <- data$ct / data$c.site.q

data$g2a.c2t <- data$g2a.rate / data$c2t.rate
data$g2a.c2t.log <- log(data$g2a.c2t,2)
data$pval <- apply(data[,c(13,12,5,6)],1,fisherTest)
data$FER <- p.adjust(data$pval,method="bonferroni")
data <- data[order(data$g2a.c2t,decreasing=T),]

write.table(data,out.name,col.names=T,row.names=F,sep="\t",quote=F)


