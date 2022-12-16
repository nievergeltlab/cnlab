 #PRS of MVP + EHR -> UKBB

#Format sumstats for training data
 zcat /mnt/ukbb/adam/tinnitus_gwas/MVP_tinnitus/eur_jul8_2021_allchr.any_tinnitus.maf01.ADD.resultsa.fuma.gz | awk '{ if (NR==1) {$1="SNP";$4="A1";$5="A2";BETA="BETA";$12="P"}; if(NR>1) BETA=log($9); print $1,toupper($4),toupper($5),BETA,$12}' > eur_jul8_2021_allchr.any_tinnitus.maf01.ADD.resultsa.fuma.gz.sst

#Get test data (UKBB) SNP list
 cat /mnt/ukbb/adam/tinnitus_gwas/bgeneur/ukbb_tinnitus_bgen_eur_unrel_*.bim > ukbb_snps.bim

#Run PRS-CS
 python /mnt/ukbb/adam/ptsd/cnv_aug_2021/PRScs-master/PRScs.py --ref_dir=/mnt/ukbb/adam/ptsd/prs_csx/ldblk_1kg_eur --bim_prefix=ukbb_snps \
 --sst_file=eur_jul8_2021_allchr.any_tinnitus.maf01.ADD.resultsa.fuma.gz.sst --n_gwas=308879   \
 --out_dir=ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25 --phi=1e-2  

#Combine all outputs (all chromosomes)
 cat ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25_pst_eff_a1_b0.5_phi1e-02_chr*.txt | awk '{print $2,$4,$6}' > ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25_pst_eff_a1_b0.5_phi1e-02_allchr.txt
 awk '{print $1}' ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25_pst_eff_a1_b0.5_phi1e-02_allchr.txt  >  ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25_pst_eff_a1_b0.5_phi1e-02_allchr.txt.snplist

#Reformat outputs for PLINK
for chr in {1..22}
do
cat ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25_pst_eff_a1_b0.5_phi1e-02_chr"$chr".txt | awk '{print $2,$4,$6}' > ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25_pst_eff_a1_b0.5_phi1e-02_chr"$chr".txt.plink
done


#Extract genotypes for the list of SNPs relevant to the PRS
 for chr in {1..22}
 do
  ../plink2 --bed /mnt/ukbb/adam/tinnitus_gwas/bgeneur/ukbb_tinnitus_bgen_eur_unrel_"$chr".bed \
 --bim /mnt/ukbb/adam/tinnitus_gwas/bgeneur/ukbb_tinnitus_bgen_eur_unrel_"$chr".bim \
 --fam /mnt/ukbb/adam/tinnitus_gwas/bgeneur/ukbb_tinnitus_bgen_eur_unrel_"$chr".fam \
 --extract   ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25_pst_eff_a1_b0.5_phi1e-02_allchr.txt.snplist  \
 --make-bed --out ukbb_mvpehr_pred_25/ukbb_mvpehr_"$chr"  
 done
 
#Sometimes there are duplicated markers. They are few, so unlikely to cause bias, so we'll just rename them. Ideally duplicates are removed from the .bim BEFORE running prs-cs!
 for chr in {1..22}
 do
 mv ukbb_mvpehr_pred_25/ukbb_mvpehr_"$chr".bim ukbb_mvpehr_pred_25/ukbb_mvpehr_"$chr".bim.bk
 done

 for chr in {1..22}
 do
 Rscript rename_bim.r ukbb_mvpehr_pred_25/ukbb_mvpehr_"$chr".bim.bk ukbb_mvpehr_pred_25/ukbb_mvpehr_"$chr".bim
 done

#Calculate PRS. PLINK1.9 used, but PLINK2 method shown and commented out (requires reformatting)
for chr in {1..22}
do

 #PLINK 1.9
 ../plink --bfile ukbb_mvpehr_pred_25/ukbb_mvpehr_"$chr"   --score ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25_pst_eff_a1_b0.5_phi1e-02_chr"$chr".txt.plink sum --out ukbb_mvpehr_pred_25/ukbb_$chr
 
  #PLINK2
 # ../plink2 --bfile ukbb_mvpehr_pred_25/ukbb_mvpehr_"$chr"   --score ukbb_mvpehr_pred_25/ukbb_mvpehr_pred_25_pst_eff_a1_b0.5_phi1e-02_chr"$chr".txt.plink --out ukbb_mvpehr_pred_25/ukbb_$chr
 #sed 's/#//g' ukbb_mvpehr_pred_25/ukbb_"$chr".sscore | awk -v CHR=$chr '{if(NR==1) TSCORE="TSCORE"CHR; if(NR>1) TSCORE=$3*$5; print $1,$2, TSCORE}' >  ukbb_mvpehr_pred_25/ukbb_"$chr".sscore.fixed
 
 
done

#You'll get a PRS (weighted sum across markers) for each chromosome. 
#The sum all chromosome PRS is the actual PRS. 
#Notice that this is based on sums and NOT averages. 
#Averages would need to be multiplied by N markers used prior to summation

#Anyway, now we'll sum the chromosome PRS to get the actual PRS, using some R tricks.  
#Here I join all files (subject as merge key), to produce a wide format file with a subject per row, chromosome specific PRS per column, then sum over rows.
#The alternative would be to concatenate all files into a long format, then use ddply to sum within subjects.
for chr in {1..22}
do
awk -v chr=$chr '{if(NR==1) $6="SCORESUM"chr; print}' ukbb_mvpehr_pred_25/ukbb_"$chr".profile >  ukbb_mvpehr_pred_25/ukbb_"$chr".profile.fixed
done

R
library(plyr)
library(data.table)

#read each file into a data frame with the same name
for (i in 1:22)
{
	assign(
		paste('ukbb_',i,sep=''), fread(paste('ukbb_mvpehr_pred_25/ukbb_',i,'.profile.fixed',sep=''), header=T, data.table=F)[,c(1,2,6)]
		) 
}

#parse the text list of data frame names as a list of data frames
data_list <- eval( 
			parse( 
				text=paste(
					"list(", paste(grep("ukbb",ls(),value=TRUE), collapse=','), ")" 
					)
				)
			)


#combine all data frames by id_visit (won't work for subjects missing this variable!!!!!)
datA <- join_all(data_list,by=c("FID","IID"), type="left", match="first")
#sum over all score columns

#datA$SCORE <- apply(datA[,-c(1,2)],1,sum) #plink2
datA$SCORE <- apply(datA[,-c(1:2)],1,sum) #plink1

#Save PRS to a file
write.table(datA,'ukbb_mvpehr_pred_25/ukbb_mvpehr_pred.profile',quote=F,row.names=F)



#Take the association between PRS and phenotype
R

#This Nagelkerke works with objects other than glm.
nagelkerke2 <- function(lintercept,lfull,n)
{
 #specify the log likelihood of the intercept model, the full model, and the n obs
 numer <- 1-(exp(2*(lintercept - lfull)/n))
 denom <- 1-exp(2*lintercept/n)
 nagelkerke=numer/denom
 return(as.numeric(nagelkerke))
}

library(data.table)
library(ordinal)

#Load PRS scores
prs <- fread('ukbb_mvpehr_pred_25/ukbb_mvpehr_pred.profile',data.table=F)

#Load principal components, phenotype data, and other covariates
pcs <- fread('/mnt/ukbb/adam/tinnitus_gwas/phenotype/UKB_tinnitus_eur_related_april2_2019.pheno',data.table=F)

#Merge PRS and phenotype related data
d2 <- merge(prs,pcs,by=c("FID","IID"),suffixes=c("","_na"))

#Rescale PRS to make it more interpretable
d2$SCALE <- scale(d2$SCORE)

#Example: ordinal logit for an ordinal phenotype
m1 <- clm(as.factor(f.4803.max_coding3) ~ f.22009.0.1+f.22009.0.2+ f.22009.0.3+ f.22009.0.4 +f.22009.0.5+ f.22009.0.6+ as.factor(f.22000.0.0) + as.factor(f.54.0.0) +SCALE,data=subset(d2))
m1a <- clm(as.factor(f.4803.max_coding3) ~ f.22009.0.1+f.22009.0.2+ f.22009.0.3+ f.22009.0.4 +f.22009.0.5+ f.22009.0.6+ as.factor(f.22000.0.0) + as.factor(f.54.0.0) ,data=subset(d2))
m1n <- clm(as.factor(f.4803.max_coding3) ~ 1,data=subset(d2))

#Here the r-square is the difference in nagelkerke r-square between the model with and without the PRS covariate
nagelkerke2(logLik(m1n),logLik(m1),172995) - nagelkerke2(logLik(m1n),logLik(m1a),172995)

# r-square = 0.007173332

#Example: binary dx
library(fmsb)
lm1 <- glm(f.4803.max_coding2 ~ f.22009.0.1+f.22009.0.2+ f.22009.0.3+ f.22009.0.4 +f.22009.0.5+ f.22009.0.6+ as.factor(f.22000.0.0) + as.factor(f.54.0.0) +SCALE,data=subset(d2),family='binomial')
lm1a <- glm(f.4803.max_coding2 ~f.22009.0.1+f.22009.0.2+ f.22009.0.3+ f.22009.0.4 +f.22009.0.5+ f.22009.0.6+ as.factor(f.22000.0.0) + as.factor(f.54.0.0) ,data=subset(d2),family='binomial')

#Can just use Nagelkerke function in fmsb library, if using GLM object.
NagelkerkeR2(lm1)$R2 - NagelkerkeR2(lm1a)$R2

# r-square = 0.007357539

#While we're at it, calculate liability scale transformation of Naglekerke R-2

#To user: Put in these values:
K=0.3676365 #Population prevalence
P=0.3676365  #Sample prevalence
h2=0.0096 #h2 on the observed scale
seh2=0.01 #standard error of h2 on the observed scale

#Do calculations
zv <- dnorm(qnorm(K))
h2_liab <- h2 * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2
var_h2_liab <- ( seh2 * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2) ^2

#Report liability scale h2snp and se 
h2_liab 
sqrt(var_h2_liab)



#Calculate quintie (5 group) odds ratios, with lowest quintile as reference
 d2$prsquantil <- cut(d2$SCORE,breaks=quantile(d2$SCORE,seq(0,1,0.2)))
 Jm1 <- as.data.frame(summary(glm(f.4803.max_coding2  ~ f.22009.0.1+f.22009.0.2+ f.22009.0.3+ f.22009.0.4 +f.22009.0.5+ f.22009.0.6+prsquantil,data=subset(d2),family='binomial'))$coefficients[8:11,1:2])

#Export regression coefficients to file for meta analysis
write.table(Jm1,file="prs_ors.txt")
