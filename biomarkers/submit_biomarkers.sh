#!/bin/bash 
#
module load R/3.6.1
module load plink2
module load zstd
today=`date +%m%d%y`
dirLog="$mr/code/multisnpnet/biomarkers/log"  # Path to save the log files
[ -d "$dirLog" ] || mkdir -p "$dirLog"

rank=(5 10 20 35)

for idx in "${!rank[@]}"; do
  sbatch -p stat,mrivas \
       --time=6:00:00 --mem=200000 \
       --output=$dirLog/mr_${today}_${rank[$idx]}.out\
       --error=$dirLog/mr_${today}_${rank[$idx]}.err \
       --job-name="multisnpnet" \
       --wrap="Rscript ./compute_metric_biomarkers.R ${rank[$idx]}"
done
