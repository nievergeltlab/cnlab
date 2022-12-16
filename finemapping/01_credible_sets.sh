#Rscript: determine width of credible sets, N snps in set

#Read complete list of SNPs in risk loci (obtained from prior script)
d1 <- read.table('finemapLD.inform.TRANS.complete.txt',stringsAsFactors=F,header=T)

library(plyr)
library(data.table)


#Give the start/stop position of credible set
range2 <- function(x)
{
r1 <- range(x)
return(c(r1[1], r1[2]))
}

credsets1 <- ddply(d1, ~ Locus, colwise(min,"BP"))
credsets2 <- ddply(d1, ~ Locus, colwise(max,"BP"))


write.table(cbind(credsets1,credsets2),file="credible_set_windows.txt",row.names=F)

#Write the number of SNPs in each credible set
credsets2 <- ddply(d1, ~ Locus, colwise(length,"BP"))

write.table(credsets2,file="credible_set_nsnps.txt",row.names=F)

#Determine number of SNPs tested
 d2 <- read.table('finemapLD.inform.EUR.complete.allsnps.txt',stringsAsFactors=F,header=T)
 credsets3 <- ddply(d2, ~ Locus, colwise(length,"BP"))
 write.table(credsets3,file="all_set_nsnps.txt",row.names=F)


#Compare number of SNPs in credible set to number of SNPS in set with r2> 0.6 with leading variant, from FUMA
d3 <- fread('eur_snps.txt',data.table=F) #SNPS.txt from FUMA. Linked based on leading SNPs
d3$SNP <- d3$rsID

dm <- merge(d3,d2,by="SNP",suffixes=c("_fuma","_finemap"))

credsets4 <- ddply(dm, ~ Locus, colwise(length,"BP"))

write.table(credsets4,file="eur_set_nsnps_fuma.txt")