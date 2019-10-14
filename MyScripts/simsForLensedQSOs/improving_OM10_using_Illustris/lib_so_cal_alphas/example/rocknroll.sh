make clean 
make 


#NCC=2048
#BSZ=61.44 #arcsec # 2.5 Mpc/h
#ZZL=0.659 # g = shear/(1-kappa)
#ZZS=2.481

NCC=4096
BSZ=368.64 
ZZL=0.55 
ZZS=2.0

./cal_alphas ../../lenses/${1}.bin ${NCC} ${BSZ} ${ZZL} ${ZZS} \
      ../../output/${1}_posx1.bin \
      ../../output/${1}_posx2.bin \
      ../../output/${1}_alpha1.bin \
      ../../output/${1}_alpha2.bin \
      ../../output/${1}_shear1.bin \
      ../../output/${1}_shear2.bin \
      ../../output/${1}_kappa.bin \
      ../../output/${1}_mu.bin \
	  ${2} lensed_${2}

