#!/usr/bin/bash

# user input

usage() {
	echo >&2 \
	"usage: install_imp_conda.sh -e conda_env_name -i install_path -s install_mode -I imp_version [-c num_cores] [-P python_version] [-m custom_module_path] [--prism] [--pmi_analysis] [--pyrmsd] [-y] [-h]

    Required arguments:
    -e conda_env_name: the name of the conda environment to be created or used
    -i install_path: the path where IMP will be installed
    -s install_mode: 'github' or 'tarball'
    -I imp_version: IMP version to install (e.g., 2.24.0) or branch name if installing from github

    Optional arguments:
    --prism: flag to indicate whether to clone the PRISM repository along with IMP (default: no)
    --pmi_analysis: flag to indicate whether to clone the PMI analysis (default: no)
    --pyrmsd: flag to indicate whether to clone the PyRMSD (default: no)
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
clone_prism="no"
clone_pmi_analysis="no"
clone_pyrmsd="no"

if grep -q "Fedora" /etc/os-release; then

    fedora_version=$(cat /etc/os-release | grep VERSION_ID | awk -F '=' '{print $2}') ;
    cgal_dir="/usr/share/cmake/CGAL" ;

else
    echo "This is not a Fedora system" ;
    echo "This installation script is designed for Fedora systems. Exiting..." ;
    exit 1 ;

fi

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -e) conda_env_name="$2"; shift ;;
        -i) install_path="$2"; shift ;;
        -s) install_mode="$2"; shift ;;
        -I) imp_version="$2"; shift ;;
        -c) ncores="$2"; shift ;;
        -P) python_version="$2"; shift ;;
        -m) custom_module_path="$2"; shift ;;
        --prism) clone_prism="yes" ;;
        --pmi_analysis) clone_pmi_analysis="yes" ;;
        --pyrmsd) clone_pyrmsd="yes" ;;
        -y) confirm_flag="true" ;;
        -h) usage; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; usage; exit 1 ;;
    esac
    shift
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

START_TIME=$(date +%s) ;

cwd_=$(pwd) ;
script_dir_="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" ;
install_path=$(realpath $install_path) ;
commit_hash_file_path="$install_path/commit_hash.txt" ;

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
            conda install -n $conda_env_name python=$python_version -y ;
            conda run -n $conda_env_name pip install matplotlib numpy scipy scikit-learn ;
        fi

    else
        echo "Installing IMP dependencies in a new conda environment $conda_env_name" ;

        conda create -n $conda_env_name python=$python_version -y ;
        source ~/.bashrc ;
        conda run -n $conda_env_name pip install matplotlib numpy scipy scikit-learn ;

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
        sleep 2 ;
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
        sleep 2 ;

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

# clone PrISM if the flag is set
if [[ "$clone_prism" == "yes" ]]; then

    echo "Cloning PRISM repository along with IMP." ;
    cd $install_path ;

    if [[ -d "./prism" ]] ; then
        echo "The PRISM repository already exists." ;
    else
        {
            git clone git@github.com:isblab/prism.git ;
        } || {
            git clone https://github.com/isblab/prism.git ;
        }
        if [ $? -ne 0 ]; then
            echo "Error cloning the PRISM repository. Please check your internet connection or the repository URL." ;
            exit 1 ;
        fi
        sleep 2 ;
        cd $install_path/prism ;
        last_commit_hash=$(git rev-parse HEAD) ;
        echo "prism commit hash: $last_commit_hash - $(date)" >> $commit_hash_file_path ;
    fi

fi

# clone PMI analysis if the flag is set
if [[ "$clone_pmi_analysis" == "yes" ]]; then

    echo "Cloning PMI analysis repository along with IMP." ;
    cd $install_path ;

    if [[ -d "./PMI_analysis" ]] ; then
        echo "The PMI analysis repository already exists." ;
    else
        {
            git clone git@github.com:salilab/PMI_analysis.git ;
        } || {
            git clone https://github.com/salilab/PMI_analysis.git ;
        }
        if [ $? -ne 0 ]; then
            echo "Error cloning the PMI analysis repository. Please check your internet connection or the repository URL." ;
            exit 1 ;
        fi

        sleep 2 ;

        cd $install_path/PMI_analysis ;
        last_commit_hash=$(git rev-parse HEAD) ;
        echo "PMI_analysis commit hash: $last_commit_hash - $(date)" >> $commit_hash_file_path ;
    fi  
fi

# clone PyRMSD if the flag is set
if [[ "$clone_pyrmsd" == "yes" ]]; then

    echo "Cloning PyRMSD repository along with IMP." ;
    cd $install_path ;

    if [[ -d "./pyRMSD" ]] ; then
        echo "The PyRMSD repository already exists." ;
    else
        {
            git clone git@github.com:salilab/pyRMSD.git ;
        } || {
            git clone https://github.com/salilab/pyRMSD.git ;
        }
        if [ $? -ne 0 ]; then
            echo "Error cloning the PyRMSD repository. Please check your internet connection or the repository URL." ;
            exit 1 ;
        fi

        sleep 2 ;

        cd $install_path/pyRMSD ;
        last_commit_hash=$(git rev-parse HEAD) ;
        echo "PyRMSD commit hash: $last_commit_hash - $(date)" >> $commit_hash_file_path ;
    fi
fi

if [[ -d "./pyRMSD" ]] ; then

    if [[ "$confirm_flag" == "true" ]]; then

        build_pyrmsd="yes"

    else

        echo "Do you want to build pyRMSD (y/n)?" ;
        read choice ;

        case "$choice" in
            y|Y ) echo "yes"; build_pyrmsd="yes";;
            n|N ) echo "no"; build_pyrmsd="no";;
            * ) echo "invalid"; build_pyrmsd="no";;
        esac

    fi

    if [ $build_pyrmsd == "yes" ]; then

        echo "Making some changes to the default.conf file for pyRMSD build." ;

        OLD_LINES='"CUDA_BASE": "/usr/local/cuda-4.2",
                "CUDA_ARCHITECHTURE": "sm_11",
                "PYTHON_LIBRARY_FOLDER": "AUTO",
                "PYTHON_LIBRARY" : "python2.7",'

        NEW_LINES='"CUDA_BASE": "/usr/local/cuda-12.8",
                "CUDA_ARCHITECHTURE": "sm_50",
                "PYTHON_LIBRARY_FOLDER": "AUTO_ALT",
                "PYTHON_LIBRARY" : "python3.12",'

        mapfile -t old_lines < <(echo "$OLD_LINES" | sed 's/^[ \t]*//') ;
        mapfile -t new_lines < <(echo "$NEW_LINES" | sed 's/^[ \t]*//') ;

        for i in "${!old_lines[@]}"; do
            
            old_line=${old_lines[i]}
            new_line=${new_lines[i]}

            echo " 
                - $old_line 
                + $new_line" ;

            sed -i "s@$old_line@$new_line@g" $install_path/pyRMSD/build_conf/default.conf ;

        done

        cd $install_path/pyRMSD ;
        echo "Building pyRMSD." ;
        conda run -n $conda_env_name python build.py --build --clean ;

    fi

fi

END_TIME=$(date +%s) ;
ELAPSED_TIME=$((END_TIME - START_TIME)) ;

echo "Installation complete. IMP and selected repositories have been installed at $install_path. Total installation time: $ELAPSED_TIME seconds." ;