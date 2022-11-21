#!/bin/bash

 
 module load 2019
 module load Boost.Python/1.67.0-intel-2019b-Python-3.6.6
 module load Tk/8.6.8-GCCcore-8.3.0 
 
 #When you start LISA:
 #may not need this step anymore
 #LD_LIBRARY_PATH=/home/maihofer/libraries/lib:$LD_LIBRARY_PATH
 
 
#Set path to where MiXeR files are stored!
 cd  /home/maihofer/freeze3_gwas/mixer
 

 for study in pgc-bip2021-all.vcf.tsv.gz PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz pgc-mdd2022-no23andMe-eur-v3.49.24.11.pgc.gz  #  eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz # eur_ptsd_pcs_v4_aug3_2021.fuma.gz eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz TotalPCL_MVP_eur.gz
 do

#Generate plots using just the subset snps 
python3 /home/maihofer/mixer/precimed/mixer_figures.py combine --json "$study".rep@.json --out "$study".mixer_results.fit
python3 /home/maihofer/mixer/precimed/mixer_figures.py one --json "$study".mixer_results.fit.json --out "$study".mixer_results.fit.plots --statistic mean std
 
#Generate plots based on whole model 
python3 /home/maihofer/mixer/precimed/mixer_figures.py combine --json "$study".mixer_results.test.rep@.json --out "$study".mixer_results.fitB
python3 /home/maihofer/mixer/precimed/mixer_figures.py one --json "$study".mixer_results.fitB.json --out "$study".mixer_results2.test.plots --statistic mean std

done


#Batch code
 # sbatch --time=1:05:00 --error errandout/mixeruniplots_13_%a.e --output errandout/mixeruniplots_13_%a.o --export=ALL 02_mixer_univariate_plots.sh 


