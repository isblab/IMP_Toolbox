# Description: This script will install a custom module and recompile IMP.
echo """

This script will install a custom module and recompile IMP.
The current build directory (if present) will be replaced with a new build directory.

"""

install_path=${1:?"Please provide the installation path where you installed imp (/full/path/to/imp_install_dir)"}
ncores=${2:?"Please provide the number of cores to use for the build (12)"}

# check for presence of imp-clean directory
if [ ! -d $install_path/imp-clean ]; then

    echo "imp-clean directory not found
    Please run install_imp.sh first to clone the imp repository"
    exit 1

fi

read -p "Do you want to install your own module and recompile IMP (y/n)?" choice

case "$choice" in
y|Y ) echo "yes"; custom_module_install="yes";;
n|N ) echo "no"; custom_module_install="no";;
* ) echo "invalid"; custom_module_install="no";;
esac

if [ $custom_module_install == "yes" ]; then

    echo "Installing custom module and recompiling IMP";

    cd $install_path/imp-clean/imp/modules ;

    read -p "Please provide the full path to the custom module: " custom_module_path

    echo "Copying custom module to imp-clean"

    if [ ! -d $custom_module_path ]; then
        echo "The custom module path does not exist. Exiting..."
        exit 1
    fi

    if [ -d $install_path/imp-clean/imp/modules/$(basename $custom_module_path) ]; then
        echo "The custom module already exists in imp
        Please remove the existing module and try again"
        exit 1
    fi

    cp -r $custom_module_path $install_path/imp-clean/imp/modules
fi

read -p "Do you want to install sampcon module (y/n)?" choice

case "$choice" in
  y|Y ) echo "yes"; sampcon_install="yes";;
  n|N ) echo "no"; sampcon_install="no";;
  * ) echo "invalid"; sampcon_install="no";;
esac

if [ $sampcon_install == "yes" ]; then

    echo "Installing sampcon module"
    cd $install_path ;

    git clone https://github.com/salilab/imp-sampcon.git
    echo "sampcon module fetched at $install_path"
    echo "Waiting for 10 seconds ..."

    sleep 10

    echo "Copying sampcon module to imp-clean"

    if [ -d $install_path/imp-clean/imp/modules/sampcon ]; then
        echo "The sampcon module already exists in imp"
        read -p "Do you want to overwrite the existing module (y/n)?" choice

        case "$choice" in
          y|Y ) echo "yes";;
          n|N ) echo "no"; exit 1;;
          * ) echo "invalid"; exit 1;;
        esac

        rm -rf $install_path/imp-clean/imp/modules/sampcon
    fi

    cp -r $install_path/imp-sampcon $install_path/imp-clean/imp/modules
    mv $install_path/imp-clean/imp/modules/imp-sampcon $install_path/imp-clean/imp/modules/sampcon

fi

cd $install_path/imp-clean ;

echo "Removing the existing build directory"
if [ -d $install_path/imp-clean/build ]; then
    rm -rf $install_path/imp-clean/build
fi

# set cgal directory based on fedora version

if grep -q "Fedora" /etc/os-release; then

    fedora_version=$(cat /etc/os-release | grep VERSION_ID | awk -F '=' '{print $2}')
    cgal_dir="/usr/share/cmake/CGAL"

else
    echo "This is not a Fedora system"
    exit 1

fi

mkdir build ;
cd build ;
cmake $install_path/imp-clean/imp -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc -DCMAKE_BUILD_TYPE=Release -DIMP_MAX_CHECKS=NONE -DIMP_MAX_LOG=SILENT -G Ninja -DCGAL_DIR=$cgal_dir ;

echo "Recompiling IMP" ;
ninja -j $ncores 1> recompile.log 2> recompile.err;
