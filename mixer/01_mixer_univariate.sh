#!/bin/bash

 
 
 module load Boost/1.72.0-GCCcore-9.3.0-no_mpi
 module load  Tkinter/3.8.2-GCCcore-9.3.0
 #pip install matplotlib_venn --user
 #pip install numdifftools --user
 
 #When you start LISA:
 #may not need this step anymore
 #LD_LIBRARY_PATH=/home/maihofer/libraries/lib:$LD_LIBRARY_PATH
 
 
#Set path to where MiXeR files are stored!
 cd  /home/maihofer/ptsd_mr/results_filtered
 
#Limited SNP subset analysis
python3 /home/maihofer/mixer/precimed/mixer.py fit1 \
      --trait1-file "$study".mixer2.gz \
      --out "$study".rep${SLURM_ARRAY_TASK_ID} \
      --extract /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps \
      --bim-file /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib  /home/maihofer/mixer/src/build/lib/libbgmg.so \

#Full SNp analysis
python3 /home/maihofer/mixer/precimed/mixer.py test1 \
      --trait1-file "$study".mixer2.gz \
      --load-params-file "$study".rep${SLURM_ARRAY_TASK_ID}.json \
      --out "$study".mixer_results.test.rep${SLURM_ARRAY_TASK_ID} \
      --bim-file /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file /home/maihofer/mixer/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib  /home/maihofer/mixer/src/build/lib/libbgmg.so \

#Example code for running 1 study

#study=141_GCST90027161-EFO_0000274-Build38.f.tsv.gz.txt
#sbatch --array=1-20 --time=5:05:00 --error errandout/mixer_13_"$study"_%a.e --output errandout/mixer_13_"$study"_%a.o --export=ALL,study="$study" 01_mixer_univariate.sh 

#example code for running many 

# for study in ukbbcoding3relatednocov_mvpanytinnitusnocov_ssw1.tbl.fuma.gz f.4803.max_coding3_no_covar_related.bgen.stats.fuma.gz eur_jul8_2021_allchr.any_tinnitus.maf01.ADD.resultsa.gz #  eur_jul8_2021_allchr.icd_tin_use.maf01.ADD.resultsa.gz eur_jul8_2021_allchr.hvin_use.maf01.ADD.resultsa.gz 
# do
 # sbatch --array=1-20 --time=5:05:00 --error errandout/mixer_13_"$study"_%a.e --output errandout/mixer_13_"$study"_%a.o --export=ALL,study="$study" 01_mixer_univariate.sh 
# done
