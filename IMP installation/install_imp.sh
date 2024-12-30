#!/bin/bash

install_path=${1:?"Please specify the path of the directory where you would like to install IMP."}
ncores=${2:?"Please specify the number of cores you would like to use for the build."}

cd $install_path ;
mkdir imp-clean ;

pip3 install --user matplotlib numpy scipy scikit-learn ;

sudo dnf install boost-devel gperftools-devel CGAL-devel graphviz gsl-devel cmake doxygen hdf5-devel swig fftw-devel opencv-devel  gcc-c++ ninja-build  python-devel ;
sudo dnf install openmpi-devel ;
sudo dnf install environment-modules ;
sudo dnf install cereal-devel
source /usr/share/Modules/init/bash

echo 'module load mpi/openmpi-x86_64' >> ~/.bash_profile ;
source ~/.bash_profile ;

cd $install_path/imp-clean/ ;

git clone -b main https://github.com/salilab/imp.git ;
cd imp && ./setup_git.py ;
cd ../
mkdir build ;
cd build ;

cmake ../imp/ -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc -DCMAKE_BUILD_TYPE=Release -DIMP_MAX_CHECKS=NONE -DIMP_MAX_LOG=SILENT -G Ninja -DCGAL_DIR=/usr/share/CGAL/cmake ;
ninja -j $ncores
