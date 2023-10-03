#!/bin/bash

#SBATCH --job-name=run_pipe
#SBATCH --constraint=A100
#SBATCH --time=0-24:00:00 
#SBATCH --ntasks-per-node=10
#SBATCH --nodes=1
#SBATCH --output=/share/nas2/ela/ASKAP/scripts/logs/slurm_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=emma.alexander@manchester.ac.uk

nvidia-smi 
echo ">>>start"

## TODO: add in multiple SBs at once? Probably simpler to run SB by SB 
##SBlist=("10007" "10040")

pwd;
conda activate jupyterenv3.8
which python
export PYTHONPATH="${PYTHONPATH}:/share/nas2/ela/ASKAP/scripts/RM-Tools/"
export PATH="${PATH}:/share/nas2/ela/ASKAP/scripts/RM-Tools/"
echo ">>>running"

SB="9962"
## run the download and cutout script (note things are currently hardcoded in them)
## TODO take SB/ urllist as input 
python /share/nas2/ela/ASKAP/scripts/askap_download_fromlist.py
python /share/nas2/ela/ASKAP/scripts/make_cutouts.py

DATADIR="/share/nas2/ela/ASKAP/data/"
RMTOOLDSIR="/share/nas2/ela/ASKAP/scripts/RM-Tools/"

## declare an array variable
declare -a ARRAY=()

## loop over the csv file of extended sources and add to the bash array
while IFS="," read -r srcname
do
  ARRAY+=($srcname)
done < <(cut -d "," -f8 /share/nas2/ela/ASKAP/data/${SB}/${SB}_extended.csv | tail -n +2)

## loop over the bash array to run rmsynth on everything 
for i in "${ARRAY[@]}"
do
   echo Working on "$i"
   cd ${DATADIR}/${SB}/"$i"
   ##python ${RMTOOLDSIR}RMtools_3D/make_freq_file.py "$i"_icube.fits "$i"_freqs.dat
   ##python ${RMTOOLDSIR}RMtools_3D/do_fitIcube.py "$i"_icube.fits -freqFile "$i"_freqs.dat -o "$i"
   ##python ${RMTOOLDSIR}RMtools_3D/do_RMsynth_3D.py -v -t -r -o "$i" -i "$i"_model.i.fits "$i"_qcube.fits "$i"_ucube.fits "$i"_freqs.dat
   python ${RMTOOLDSIR}RMtools_3D/do_RMclean_3D.py -v -o "$i" -c 0.0002 "$i"FDF_tot_dirty.fits "$i"RMSF_tot.fits
   echo "Done??"
done
echo "All Done??"

## TODO: post rm-synth plotting and processing script 