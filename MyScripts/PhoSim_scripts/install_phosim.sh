wget https://bitbucket.org/phosim/phosim_release/get/d912a460f878.zip
unzip d912a460f878.zip

#cfitsio
##=============================
##Install FFTW3
##=============================
#wget http://www.fftw.org/fftw-3.3.4.tar.gz
#tar xzvf fftw-3.3.4.tar.gz
#cd fftw-3.3.4
#./configure --prefix=$GAL_DIR/deps CFLAGS='-fPIC'
#make
#make install
#cd ../
##=============================

./configure
make

./phosim examples/star -c examples/nobackground
