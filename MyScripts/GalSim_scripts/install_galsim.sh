#mkdir galsim
#cd galsim
#mkdir deps
#wget http://sourceforge.net/projects/scons/files/scons/2.3.4/scons-2.3.4.tar.gz
#wget http://sourceforge.net/projects/boost/files/boost/1.58.0/boost_1_58_0.tar.gz
#wget http://www.fftw.org/fftw-3.3.4.tar.gz
#
#tar xzvf scons-2.3.4.tar.gz
#tar xzvf boost_1_58_0.tar.gz
#tar xzvf fftw-3.3.4.tar.gz
#
#svn checkout http://tmv-cpp.googlecode.com/svn/tags/v0.72/ tmv-cpp
git clone https://github.com/GalSim-developers/GalSim.git

#GAL_DIR=`pwd -P`
##GAL_DIR=/global/project/projectdirs/hacc/nanli/galsim
#echo $GAL_DIR

##=============================
##Install Scons
##=============================
#
#cd scons-2.3.4
#python setup.py install --prefix=$GAL_DIR/deps
#cd ../
#
##=============================
##Install FFTW3
##=============================
#
#cd fftw-3.3.4
#./configure --prefix=$GAL_DIR/deps CFLAGS='-fPIC'
#make
#make install
#cd ../
#
##=============================
##Install boost_1_58_0
##=============================
#
#cd boost_1_58_0
#./bootstrap.sh --with-python=/Users/uranus/anaconda/bin/python
#./b2 link=shared
#./b2 --prefix=$GAL_DIR/deps link=shared install
#cd ../
#
##=============================
##Install TMV
##=============================
#
#cd tmv-cpp
#../deps/bin/scons install PREFIX=$GAL_DIR/deps
#cd ../
#
##=============================
##Install GalSim
##=============================
#
#cd GalSim
#$GAL_DIR/deps/bin/scons TMV_DIR=$GAL_DIR/deps FFTW_DIR=$GAL_DIR/deps BOOST_DIR=$GAL_DIR/deps
#
#$GAL_DIR/deps/bin/scons install PREFIX=$GAL_DIR/deps PYPREFIX=$HOME/.local/lib/python2.7/site-packages
#
#$GAL_DIR/deps/bin/scons tests
#$GAL_DIR/deps/bin/scons examples
#cd ../
#
#install_name_tool -id /anaconda/lib/libpython2.7.dylib /anaconda/lib/libpython2.7.dylib 
#brew install boost-python boost

cd GalSim
scons install PYTHON=/Users/uranus/anaconda/bin/python PYPREFIX=/Users/uranus/.local/lib/python2.7/site-packages
scons tests
scons examples

#=============================
#Uninstall GalSim
#=============================
#
#scons uninstall
#scons -c
