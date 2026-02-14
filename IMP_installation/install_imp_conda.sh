#!/usr/bin/bash

# user input

usage() {
	echo >&2 \
	"usage: install_imp_conda.sh -e conda_env_name -i install_path -s install_mode -I imp_version [-c num_cores] [-P python_version] [-m custom_module_path] [-y] [-h]

Required arguments:
-e conda_env_name: the name of the conda environment to be created or used
-i install_path: the path where IMP will be installed
-s install_mode: 'github' or 'tarball'
-I imp_version: IMP version to install (e.g., 2.24.0) or branch name if installing from github

Optional arguments:
-c num_cores: number of cores to use for compilation (default: half the number of cores available on the system)
-P python_version: Python version to use for installation (e.g., 3.12)
-m custom_module_path: the path to a custom IMP module that you want to install along with IMP.
-y: flag to skip confirmation prompts and use default values for environment name and installation choices
-h: display this message"
}

#defaults
conda_env_name=
install_path=
install_mode=
imp_version=
ncores=$(($(nproc) / 2))
python_version=3.12
custom_module_path=""
confirm_flag="false"

if grep -q "Fedora" /etc/os-release; then

    fedora_version=$(cat /etc/os-release | grep VERSION_ID | awk -F '=' '{print $2}') ;
    cgal_dir="/usr/share/cmake/CGAL" ;

else
    echo "This is not a Fedora system" ;
    echo "This installation script is designed for Fedora systems. Exiting..." ;
    exit 1 ;

fi

while getopts "e:i:s:I:c:P:yhm:" opt;
do
    case "$opt" in
        e) conda_env_name="$OPTARG" ;;
        i) install_path="$OPTARG" ;;
        s) install_mode="$OPTARG" ;;
        I) imp_version="$OPTARG" ;;
        c) ncores="$OPTARG" ;;
        P) python_version="$OPTARG" ;;
        m) custom_module_path="$OPTARG" ;;
        y) confirm_flag="true" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

if [ -z "$conda_env_name" ] || [ -z "$install_path" ] || [ -z "$install_mode" ] || [ -z "$imp_version" ]; then
    echo "Missing required arguments. Please provide values for conda_env_name, install_path, install_mode, and imp_version."
    usage
    exit 1
fi

if [ "$install_mode" != "github" ] && [ "$install_mode" != "tarball" ]; then
    echo "Invalid installation mode. Please specify 'github' or 'tarball'."
    exit 1
fi

cwd_=$(pwd) ;
script_dir_="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" ;
install_path=$(realpath $install_path) ;

# conda environment setup

# Check for a flag to use the default value for env_choice
if [[ "$confirm_flag" == "true" ]]; then

    choice=$conda_env_name

else

    read -p """
    This script will install IMP from the source code.

    A new conda environment $conda_env_name will be created for the installation.
    All the python dependencies required for IMP will be installed in this environment.

    Please confirm the environment name if you want to continue?: """ choice

fi

if [ $choice == $conda_env_name ]; then

    if conda info --envs | grep -q "^$choice\b" ; then

        # Check for a flag to use the default value for env_choice
        if [[ "$confirm_flag" == "true" ]]; then

            env_choice="y"

        else
            echo "environment already exists."
            read -p """

            WARNING: The environment already exists.

            Do you want to continue installation in the existing environment (y/n)?: """ env_choice

            read -p ""
        fi

        case "$env_choice" in

        y|Y ) echo "yes"; to_continue="yes";;
        n|N ) echo "no"; to_continue="no";;
        * ) echo "invalid"; to_continue="no";;

        esac

        if [ $to_continue == "no" ]; then
            echo "Exiting..." ;
            exit 1 ;

        elif [ $to_continue == "yes" ]; then
            echo "Installing IMP dependencies in existing conda environment $conda_env_name" ;
            source ~/.bashrc ;
            conda activate $conda_env_name ;
            conda install python=$python_version matplotlib numpy scipy scikit-learn -y ;
        fi

    else
        echo "Installing IMP dependencies in a new conda environment $conda_env_name" ;

        conda create -n $conda_env_name python=$python_version matplotlib numpy scipy scikit-learn -y ;
        source ~/.bashrc ;

        if conda info --envs | grep -q "^$choice\b" ; then
            echo "environment exists." ;
        else
            echo "There was an error creating the conda environment. Exiting..." ;
            exit 1 ;
        fi
    fi

else
    echo "The environment name does not match. Exiting..." ;
    exit 1 ;

fi

echo $CONDA_DEFAULT_ENV ;

# install required system packages

cd $install_path ;
if [[ -d "./imp-clean" ]]; then
    echo "The imp-clean directory already exists.";
else
    mkdir imp-clean ;
fi

sudo dnf install boost-devel gperftools-devel CGAL-devel graphviz gsl-devel cmake doxygen hdf5-devel swig fftw-devel opencv-devel gcc-c++ ninja-build python-devel ;
sudo dnf install openmpi-devel ;
sudo dnf install environment-modules ;
sudo dnf install cereal-devel ;
source /usr/share/Modules/init/bash ;

echo 'module load mpi/openmpi-x86_64' >> ~/.bash_profile ;
source ~/.bash_profile ;

cd $install_path/imp-clean/ ;

# IMP source code download

if [ $install_mode == "tarball" ]; then
    echo "Downloading IMP $imp_version tarball";

    tarball_name="imp-$imp_version" ;

    if [[ -f "$tarball_name.tar.gz" ]]; then
        echo "The IMP tarball already exists, will not download again.";
    else
        echo "Downloading IMP tarball version $imp_version";
        wget https://integrativemodeling.org/$imp_version/download/$tarball_name.tar.gz ;
        wget https://integrativemodeling.org/$imp_version/download/SHA256SUM ;
        if [ $? -ne 0 ]; then
            echo "Error downloading the IMP tarball. Please check the version number and your internet connection." ;
            exit 1 ;
        fi
        python $script_dir_/check_sha256sum.py --file_path $tarball_name.tar.gz --sha256_file SHA256SUM ;
        if [ $? -eq 0 ]; then
            echo "SHA256 checksum verification passed." ;
        else
            echo "SHA256 checksum verification failed. Exiting." ;
            exit 1 ;
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
        whether_to_install="no" ;
        echo "The imp repository already exists. Checking version." ;
        github_imp_version=$(cat $install_path/imp-clean/imp/VERSION) ;

        if [ "$github_imp_version" == "$imp_version" ]; then
            echo "The cloned GitHub IMP version ($github_imp_version) matches the requested version ($imp_version). Will not clone again." ;
            whether_to_install="no" ;

        else
            echo "The cloned GitHub IMP version ($github_imp_version) does not match the requested version ($imp_version). Re-cloning the repository." ;
            whether_to_install="yes" ;
            sudo rm -rf imp ;

        fi

    else
        whether_to_install="yes" ;
    fi

    if [ "$whether_to_install" == "yes" ]; then

        echo "Cloning IMP from github repository.";

        {
            git clone --recurse-submodules --branch $imp_version git@github.com:salilab/imp.git ;
        } || {
            git clone --recurse-submodules --branch $imp_version https://github.com/salilab/imp.git ;
        }
        if [ $? -ne 0 ]; then
            echo "Error cloning the IMP repository. Please check your internet connection or the repository URL." ;
            exit 1 ;
        fi
        wait 2 ;

    fi

    cd imp ;
    git submodule update --init --recursive ;
    ./setup_git.py ;

    imp_src_path=$install_path/imp-clean/imp ;

else 
    echo "Invalid installation mode. Exiting.";
    exit 1 ;

fi

# custom module installation

if [[ -n "$custom_module_path" ]]; then

    echo "Custom module path provided as an argument. Will attempt to install the custom module from the provided path." ;

    if [[ ! -d $custom_module_path ]]; then
        echo "The custom module path does not exist. Exiting..." ;
        exit 1 ;
    fi

    custom_module_install="yes" ;

else

    read -p "Do you want to install your own module and recompile IMP (y/n)?" choice

    case "$choice" in
    y|Y ) echo "yes"; custom_module_install="yes";;
    n|N ) echo "no"; custom_module_install="no";;
    * ) echo "invalid"; custom_module_install="no";;
    esac

    if [ $custom_module_install == "yes" ]; then

        read -p "Please provide the full path to the custom module: " custom_module_path
    
        if [[ ! -d $custom_module_path ]]; then
            echo "The custom module path does not exist. Exiting..." ;
            exit 1 ;
        fi
    fi

fi

if [ $custom_module_install == "yes" ]; then

    echo "Installing custom module and recompiling IMP" ;

    cd $imp_src_path/modules ;

    echo "Copying custom module to imp-clean" ;

    if [[ -d $imp_src_path/modules/$(basename $custom_module_path) ]]; then
        echo "The custom module already exists in imp" ;
        echo "Removing the existing custom module from imp for fresh installation" ;
        rm -rf $imp_src_path/modules/$(basename $custom_module_path) ;
    fi

    cp -r $custom_module_path $imp_src_path/modules ;
fi

echo "Custom module installation complete." ;

# Recompiling IMP

cd $install_path/imp-clean ;

if [[ -d "./build" ]]; then
    echo "The build directory already exists, removing it for fresh compilation.";
    rm -rf build ;
fi
mkdir build ;
cd build ;

if [ "$$CONDA_DEFAULT_ENV" != "$conda_env_name" ]; then
    conda activate $conda_env_name ;
fi

cmake $imp_src_path -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc -DCMAKE_BUILD_TYPE=Release -DIMP_MAX_CHECKS=NONE -DIMP_MAX_LOG=SILENT -G=Ninja -DCGAL_DIR=$cgal_dir ;

echo "Compiling IMP" ;
ninja -j $ncores 1> recompile.log 2> recompile.err ;

sudo chown -R $USER:$USER $install_path/imp-clean ;