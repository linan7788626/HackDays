#!/bin/bash
#for i in $( ls ./finalcats ); do

#python ./galsim_no_pixel_run.py -i ./$1/$2 -o ./$1/$2 -f g  
#python ./galsim_no_pixel_run.py -i ./$1/$2 -o ./$1/$2 -f r 
#python ./galsim_no_pixel_run.py -i ./$1/$2 -o ./$1/$2 -f i 
#python ./galsim_no_pixel_run.py -i ./$1/$2 -o ./$1/$2 -f z

python ./galsim_no_pixel_run.py -i ./lens_gals.cat -o gals_lens -f g
python ./galsim_no_pixel_run.py -i ./lens_gals.cat -o gals_lens -f r
python ./galsim_no_pixel_run.py -i ./lens_gals.cat -o gals_lens -f i
python ./galsim_no_pixel_run.py -i ./lens_gals.cat -o gals_lens -f z

