
#install polyfun
 git clone https://github.com/omerwe/polyfun
 cd polyfun
 sudo /usr/local/bin/anaconda3/condabin/conda env create -f polyfun.yml
 python test_polyfun.py #seems to run


#Start polyfun
conda activate polyfun
 
 
#Munge Sumstats

#Premunge - capitalize A1 and A2, name MAF and N columns, change X to 23.

 #FYI, The X chromosome is not present in default Polyfun reference data
 zcat ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz | awk '{if(NR==1) { $3="SNP"; $4="A1";$5="A2";$6="MAF";$7="N"; $8="Z";$9="P"}; if (NR>1) {$4=toupper($4); $5=toupper($5)};  print}' | grep -v X > ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.premunge

#Munge
 python /home/genetics/polyfun/munge_polyfun_sumstats.py --sumstats ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.premunge --out ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun

#Calculate priors - using approach 1 of pre-computed priors
 python /home/genetics/polyfun/extract_snpvar.py --sumstats ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun --out ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun_prior
 
 
#Remove the not found SNPs. Alternatively, just exclude them using a flag in the munge script..
R
library(data.table)
d1 <- fread('ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.premunge',data.table=F)
d2 <- fread('zcat ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun_prior.miss.gz', data.table=F)

d1a <- subset(d1,!(SNP %in% d2$SNP))

write.table(d1a,file='ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.premunge',quote=F,row.names=F)

#Now try to munge again
 python /home/genetics/polyfun/munge_polyfun_sumstats.py --sumstats ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.premunge --out ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun

#Again, attempt to calculate priors
 python /home/genetics/polyfun/extract_snpvar.py --sumstats ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun --out ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun_prior

#Make a directory for outputs
 mkdir output

#This file should be in .csv format, and contain the locus name, chromosome, start position, stop position, and name of the LD file that contains pairwise LD for all markers in the locus (see accompanying excel sheet)
 dos2unix tinnitus_finemaplist_withrefs.csv


IFS=$'\n'
for snpset in $(cat tinnitus_finemaplist_withrefs.csv )
do

snp=$(echo $snpset | awk 'BEGIN{FS=","}{print $1}')
chr=$(echo $snpset | awk 'BEGIN{FS=","}{print $2}')
start=$(echo $snpset | awk 'BEGIN{FS=","}{print $3}')
stop=$(echo $snpset | awk 'BEGIN{FS=","}{print $4}')
fileld=$(echo $snpset | awk 'BEGIN{FS=","}{print $5}')

 #Download LD information, warning: it's big
  wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB_LD/"$fileld".gz
  wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB_LD/"$fileld".npz
 
 #Single causal locus (no ld used)
   python /home/genetics/polyfun/finemapper.py \
   --sumstats ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun_prior \
   --chr $chr \
   --n 470336 \
   --start  $start 	  \
   --end $stop   \
   --method susie \
   --max-num-causal 1 \
   --out output/finemapLD_1.inform.EUR."$snp"."$chr"."$start"."$stop".gz
   
  #Two causal loci (LD used)
   python /home/genetics/polyfun/finemapper.py \
       --ld "$fileld" \
   --sumstats ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz.polyfun_prior \
   --chr $chr \
   --n 470336 \
   --start  $start 	  \
   --end $stop   \
   --method susie \
   --max-num-causal 2 \
   --out output/finemapLD_2.inform.EUR."$snp"."$chr"."$start"."$stop".gz
    
    
  #remove LD files if you're not going to use them anymore..   
   rm "$fileld".gz "$fileld".npz
   
      
done

#Given fine-mapping for all loci, put results into a concatenatable format
IFS=$'\n'
for snpset in $(cat tinnitus_finemaplist_withrefs.csv )
do

snp=$(echo $snpset | awk 'BEGIN{FS=","}{print $1}')
chr=$(echo $snpset | awk 'BEGIN{FS=","}{print $2}')
start=$(echo $snpset | awk 'BEGIN{FS=","}{print $3}')
stop=$(echo $snpset | awk 'BEGIN{FS=","}{print $4}')
fileld=$(echo $snpset | awk 'BEGIN{FS=","}{print $5}')

#Get just markers in credible sets
zcat output/finemapLD_2.inform.EUR."$snp"."$chr"."$start"."$stop".gz | awk -v leadsnp=$snp '{if (NR==1 || $15=="1") print leadsnp, $0}' > output_filtered/finemapLD.inform.EUR."$snp"."$chr"."$start"."$stop".credible

#Get results for all markers
zcat output/finemapLD_2.inform.EUR."$snp"."$chr"."$start"."$stop".gz | awk -v leadsnp=$snp '{print leadsnp, $0}' > output_filtered/finemapLD.inform.EUR."$snp"."$chr"."$start"."$stop".allsnps

done


#Concatenate credible SNPs only
cat output_filtered/*.credible | awk '{if(NR==1){$1="Locus"}; if(NR==1 || $2 != "CHR") print}' > finemapLD.inform.EUR.complete.txt

#Concatenate all SNPs
cat output_filtered/*.allsnps | awk '{if(NR==1){$1="Locus"}; if(NR==1 || $2 != "CHR") print}'   > finemapLD.inform.EUR.complete.allsnps.txt

