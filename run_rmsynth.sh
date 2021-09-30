#!/bin/bash 
## declare an array variable
##declare -a arr=("1709-2226" "1712-2435" "1713-2502")
##declare -a arr=("2043-5059" "A3716" "2046-5107" "2058-4838")
## now loop through the above array
##declare -a arr=("2018-5540")
declare -a arr=("2151-5520")
for i in "${arr[@]}"
do
   echo Working on "$i"

   cd "$i"

   ##mkdir old_fdf_files
   ##mv FDF_*.fits old_fdf_files
   ##mv RMSF_*fits old_fdf_files

   ##mkdir clean_0.002
   ##mv FDF_CC* clean_0.002/
   ##mv FDF_clean* clean_0.002/

   python ../../RM-Tools-master/RMtools_3D/do_RMsynth_3D.py -i "$i"_i.fits -d 5 -l 2000 -v "$i"_q.fits "$i"_u.fits ../freqs_Hz_SB10040.dat
   python ../../RM-Tools-master/RMtools_3D/do_RMclean_3D.py -v -c 0.0001 FDF_tot_dirty.fits RMSF_tot.fits

   cd ..
done

python /Users/emma/OneDrive/PhD/possum_tools/faraday_spectrum_fit.py 
python /Users/emma/OneDrive/PhD/possum_tools/possum_plots_boogaloo.py
