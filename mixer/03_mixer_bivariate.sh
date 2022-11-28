#!/bin/bash

 
 
 module load Boost/1.72.0-GCCcore-9.3.0-no_mpi
 module load  Tkinter/3.8.2-GCCcore-9.3.0
 #pip install matplotlib_venn --user
 #pip install numdifftools --user
 
 #When you start LISA:
 #may not need this step anymore
 #LD_LIBRARY_PATH=/home/maihofer/libraries/lib:$LD_LIBRARY_PATH
 


#Set path to where MiXeR files are stored!
cd  /home/maihofer/freeze3_gwas/mixer

##Convert sumstats

#Studies:
 # study1=eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz 
 # study2=eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz


studycomb="$study1"_"$study2"
python3 /home/maihofer/mixer/precimed/mixer.py fit2 \
      --trait1-file "$study1".mixer2.gz \
      --trait2-file "$study2".mixer2.gz \
      --trait1-params-file "$study1".rep${SLURM_ARRAY_TASK_ID}.json \
      --trait2-params-file "$study2".rep${SLURM_ARRAY_TASK_ID}.json \
      --out "$studycomb".rep${SLURM_ARRAY_TASK_ID} \
      --extract  /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps \
      --bim-file /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file  /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib  /home/maihofer/mixer/src/build/lib/libbgmg.so \
      

python3 /home/maihofer/mixer/precimed/mixer.py  test2 \
      --trait1-file "$study1".mixer2.gz \
      --trait2-file "$study2".mixer2.gz \
      --load-params-file "$studycomb".rep${SLURM_ARRAY_TASK_ID}.json \
      --out "$studycomb".test.rep${SLURM_ARRAY_TASK_ID} \
      --bim-file /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib  /home/maihofer/mixer/src/build/lib/libbgmg.so 
     
#Batch submission code:     
# sbatch --array=1-20  --time=16:05:00 --error errandout/bivariate_"$study1"_"$study2"_%a.e --output errandout/bivariate_"$study1"_"$study2"_%a.o --export=ALL,study1=$study1,study2=$study2 03_mixer_bivariate.sh 
