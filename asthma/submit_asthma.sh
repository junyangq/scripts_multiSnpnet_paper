#!/bin/bash 
#
module load R/3.6.1
module load plink2
module load zstd
today=`date +%m%d%y`
dirLog="$mr/code/multisnpnet/log/asthma"  # Path to save the log files
[ -d "$dirLog" ] || mkdir -p "$dirLog"

rank=(2 3 4 5 6 7 8)

for idx in "${!rank[@]}"; do
  sbatch -p stat,mrivas \
       --time=4:00:00 --mem=200000 \
       --output=$dirLog/mr_${today}_${rank[$idx]}.out\
       --error=$dirLog/mr_${today}_${rank[$idx]}.err \
       --job-name="multisnpnet" \
       --wrap="Rscript ./compute_metric_asthma.R ${rank[$idx]}"
done
