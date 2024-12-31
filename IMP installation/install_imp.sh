#!/bin/bash

install_path=${1:?"Please specify the path of the directory where you would like to install IMP."}
ncores=${2:?"Please specify the number of cores you would like to use for the build."}

cd $install_path ;
mkdir imp-clean ;

pip3 install --user matplotlib numpy scipy scikit-learn ;

sudo dnf install boost-devel gperftools-devel CGAL-devel graphviz gsl-devel cmake doxygen hdf5-devel swig fftw-devel opencv-devel  gcc-c++ ninja-build  python-devel ;
sudo dnf install openmpi-devel ;
sudo dnf install environment-modules ;
sudo dnf install cereal-devel ;
source /usr/share/Modules/init/bash ;

echo 'module load mpi/openmpi-x86_64' >> ~/.bash_profile ;
source ~/.bash_profile ;

cd $install_path/imp-clean/ ;
git clone -b main https://github.com/salilab/imp.git ;
# git clone -b main git@github.com:salilab/imp.git ;
cd imp && ./setup_git.py ;
cd ../ ;
mkdir build ;
cd build ;

# set cgal directory based on fedora version
# ref: https://integrativemodeling.org/2.22.0/download/IMP.spec

if grep -q "Fedora" /etc/os-release; then

    fedora_version=$(cat /etc/os-release | grep VERSION_ID | awk -F '=' '{print $2}')

    if [ $fedora_version -ge 29 ] && [ $fedora_version -lt 32 ]; then
        cgal_dir="%{_libdir}/cmake/CGAL"

    elif [ $fedora_version -eq 32 ]; then
        cgal_dir="/usr/share/CGAL/cmake"

    elif [ $fedora_version -ge 33 ]; then
        cgal_dir="/usr/share/cmake/CGAL"

    fi

else
    echo "This is not a Fedora system"
    exit 1

fi

cmake ../imp/ -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc -DCMAKE_BUILD_TYPE=Release -DIMP_MAX_CHECKS=NONE -DIMP_MAX_LOG=SILENT -G Ninja -DCGAL_DIR=$cgal_dir ;
ninja -j $ncores
