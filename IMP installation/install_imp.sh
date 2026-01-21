#!/bin/bash

install_path=${1:?"Please specify the path of the directory where you would like to install IMP."}
ncores=${2:?"Please specify the number of cores you would like to use for the build."}
install_mode=${3:?"Please specify the installation source: 'github' or 'tarball'."}
if [ "$install_mode" != "github" ] && [ "$install_mode" != "tarball" ]; then
    echo "Invalid installation mode. Please specify 'github' or 'tarball'."
    exit 1
fi
imp_version=${4:?"Please provide the IMP version to install (e.g., 2.23.0). or branch name if installing from github."}

cwd_=$(pwd) ;

cd $install_path ;
if [[ -d "./imp-clean" ]]; then
    echo "The imp-clean directory already exists.";
else
    mkdir imp-clean ;
fi

pip3 install --user matplotlib numpy scipy scikit-learn ;

sudo dnf install boost-devel gperftools-devel CGAL-devel graphviz gsl-devel cmake doxygen hdf5-devel swig fftw-devel opencv-devel  gcc-c++ ninja-build  python-devel ;
sudo dnf install openmpi-devel ;
sudo dnf install environment-modules ;
sudo dnf install cereal-devel ;
source /usr/share/Modules/init/bash ;

echo 'module load mpi/openmpi-x86_64' >> ~/.bash_profile ;
source ~/.bash_profile ;

cd $install_path/imp-clean/ ;

if [ $install_mode == "tarball" ]; then
    echo "Downloading IMP $imp_version tarball";

    tarball_name="imp-$imp_version"

    if [[ -f "$tarball_name.tar.gz" ]]; then
        echo "The IMP tarball already exists, will not download again.";
    else
        echo "Downloading IMP tarball version $imp_version";
        wget https://integrativemodeling.org/$imp_version/download/$tarball_name.tar.gz ;
        wget https://integrativemodeling.org/$imp_version/download/SHA256SUM ;
        if [ $? -ne 0 ]; then
            echo "Error downloading the IMP tarball. Please check the version number and your internet connection."
            exit 1
        fi
        python $cwd_/check_sha256sum.py --file_path $tarball_name.tar.gz --sha256_file SHA256SUM ;
        if [ $? -eq 0 ]; then
            echo "SHA256 checksum verification passed."
        else
            echo "SHA256 checksum verification failed. Exiting."
            exit 1
        fi
        wait 2 ;
    fi

    if [[ -d "./$tarball_name" ]]; then
        echo "The IMP source directory already exists, removing it.";
        rm -rf $tarball_name ;
    fi

    tar -xvzf $tarball_name.tar.gz ;
    cd $tarball_name ;

    imp_src_path=$install_path/imp-clean/$tarball_name ;

elif [ $install_mode == "github" ]; then

    if [[ -d "./imp" ]] ; then
        echo "The imp repository already exists, will not clone again.";

    else

        echo "Cloning IMP from github repository.";

        {
            git clone --branch $imp_version git@github.com:salilab/imp.git ;
        } || {
            git clone --branch $imp_version https://github.com/salilab/imp.git ;
        }
        if [ $? -ne 0 ]; then
            echo "Error cloning the IMP repository. Please check your internet connection or the repository URL."
            exit 1
        fi
        wait 2 ;

    fi

    cd imp ;
    ./setup_git.py ;

    imp_src_path=$install_path/imp-clean/imp ;

else 
    echo "Invalid installation mode. Exiting.";
    exit 1 ;

fi

read -p "Do you want to install your own module and recompile IMP (y/n)?" choice

case "$choice" in
y|Y ) echo "yes"; custom_module_install="yes";;
n|N ) echo "no"; custom_module_install="no";;
* ) echo "invalid"; custom_module_install="no";;
esac

if [ $custom_module_install == "yes" ]; then

    echo "Installing custom module and recompiling IMP";

    cd $imp_src_path/modules ;

    read -p "Please provide the full path to the custom module: " custom_module_path

    if [[ ! -d $custom_module_path ]]; then
        echo "The custom module path does not exist. Exiting..."
        exit 1
    fi

    echo "Copying custom module to imp-clean"

    if [[ -d $imp_src_path/modules/$(basename $custom_module_path) ]]; then
        echo "The custom module already exists in imp"
        read -p "Do you want to overwrite the existing module (y/n)?" choice
        case "$choice" in
          y|Y ) echo "yes";;
          n|N ) echo "no"; exit 1;;
          * ) echo "invalid"; exit 1;;
        esac
        rm -rf $imp_src_path/modules/$(basename $custom_module_path)
    fi

    cp -r $custom_module_path $imp_src_path/modules
fi

cd $install_path/imp-clean ;

if [[ -d "./build" ]]; then
    echo "The build directory already exists, removing it for fresh compilation.";
    rm -rf build ;
fi
mkdir build ;
cd build ;

# set cgal directory based on fedora version
# ref: https://integrativemodeling.org/2.22.0/download/IMP.spec

if grep -q "Fedora" /etc/os-release; then

    fedora_version=$(cat /etc/os-release | grep VERSION_ID | awk -F '=' '{print $2}')
    cgal_dir="/usr/share/cmake/CGAL"

else
    echo "This is not a Fedora system"
    exit 1

fi

cmake $imp_src_path -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc -DCMAKE_BUILD_TYPE=Release -DIMP_MAX_CHECKS=NONE -DIMP_MAX_LOG=SILENT -G=Ninja -DCGAL_DIR=$cgal_dir ;
ninja -j $ncores 1> recompile.log 2> recompile.err;
