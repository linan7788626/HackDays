#ds9 -geometry 1024x768 $fitsfile -scale limits 0.0 4.0 -zoom to fit -view colorbar no -saveimage png $fitsfile:r.png -exit

#mencoder mf://*.png -mf w=800:h=600:fps=5:type=png -ovc raw -oac copy -o output.avi

#import numpy as np
import subprocess as sp

#files_dir = "./output_fits/"
#videos_dir = "./output_fits/output_movie/"
#cmd1 = "ls "+files_dir+" |grep output_double.fits"
#filename_list = sp.check_output(cmd1,shell=True)
#filenames = filename_list.split("\n")[:-1]

#for i in filenames:
    #cmd2 = "ds9 -geometry 1024x768 "+files_dir+i+" -scale limits 0.0 4.0 -zoom to fit -view colorbar no -saveimage png " + videos_dir+ i[:-5]+".png -exit"
    #sp.call(cmd2,shell=True)

cmd3 = "cd ./output_fits/output_movie/ && mencoder mf://*.png -mf w=800:h=600:fps=5:type=png -ovc raw -oac copy -o output.avi"

sp.call(cmd3,shell=True)
