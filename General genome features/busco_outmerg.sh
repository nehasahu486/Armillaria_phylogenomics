#!/bin/bash
#SBATCH --job-name=busco_output_merger
#SBATCH -N 
#SBATCH -w 
#SBATCH --ntasks-per-node
#SBATCH --no-requeue

MERGED="/directory/where/.txt_files/will_be_saved/"
for OUT_TEXT in /directory/with/busco_output_folders/*/;
do
 for OUPUT in $OUT_TEXT*.txt;
 do
  cp $OUPUT $MERGED/$(echo $(basename $OUT_TEXT)"_"$(basename $OUPUT))
  done
done