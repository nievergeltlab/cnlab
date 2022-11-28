
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

 study1=eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz 
 study2=eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz
 

# study2=pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz
  
#study1=


#

study2=pgc-mdd2022-no23andMe-eur-v3.49.24.11.pgc.gz

eur_ptsd_pcs_v4_aug3_2021.fuma.gz 
eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz
 TotalPCL_MVP_eur.gz
 
 #Singular matrix for EHR, case/control
 
for study1 in eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz   eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz
do
echo $study1
studycomb="$study1"_"$study2"

#Bivarate plot based on fit
  python /home/maihofer/mixer/precimed/mixer_figures.py combine --json "$studycomb".rep@.json --out      "$studycomb".rep.fit         

#Bivariate plot based on total    
  python /home/maihofer/mixer/precimed/mixer_figures.py combine --json "$studycomb".test.rep@.json --out "$studycomb".test.rep.fit

#Make plot
  python /home/maihofer/mixer/precimed/mixer_figures.py two --json-fit "$studycomb".rep.fit.json --json-test "$studycomb".test.rep.fit.json  --out "$studycomb".test.rep.combine  --statistic mean std
     
done

#Run this in the login shell, no need to submit job..
