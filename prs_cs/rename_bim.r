library(data.table)
args <- commandArgs(trailingOnly = TRUE)
 infile <- args[1]
 outfile <- args[2]
 
 d1 <- fread(infile,data.table=F)
 try({
 d1[duplicated(d1$V2),]$V2 <- paste('RSrs',d1[duplicated(d1$V2),]$V4,sep='')
 d1[duplicated(d1$V2),]$V2 <- paste(d1[duplicated(d1$V2),]$V2,'_dup',sep='')
  d1[duplicated(d1$V2),]$V2 <- paste(d1[duplicated(d1$V2),]$V2,'_dup',sep='')
 d1[duplicated(d1$V2),]$V2 <- paste(d1[duplicated(d1$V2),]$V2,'_dup',sep='')}
,silent=TRUE)

 write.table(d1,file=outfile,quote=F,row.names=F,col.names=F)

