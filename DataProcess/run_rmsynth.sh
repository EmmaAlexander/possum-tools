#!/bin/bash 
## declare an array variable
declare -a arr=("2019-5254")

for i in "${arr[@]}"
do
   echo Working on "$i"

   cd "$i"

   python /Users/emma/GitHub/RM-Tools/RMtools_3D/make_freq_file.py "$i"_POSSUM_i.fits "$i"_freqs.dat
   python /Users/emma/GitHub/RM-Tools/RMtools_3D/do_fitIcube.py "$i"_i.fits "$i"_freqs.dat
   python /Users/emma/GitHub/RM-Tools/RMtools_3D/do_RMsynth_3D.py -v -i Imodel.fits -n Inoise.dat -w 'variance' "$i"_POSSUM_q.fits "$i"_POSSUM_u.fits "$i"_freqs.dat
   python /Users/emma/GitHub/RM-Tools/RMtools_3D/do_RMclean_3D.py -v -c 0.0002 FDF_tot_dirty.fits RMSF_tot.fits

   cd ..
done

