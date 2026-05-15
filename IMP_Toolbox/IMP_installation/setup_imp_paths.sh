#!/usr/bin/bash

# user input

usage() {
	echo >&2 \
    "Usage: $0 -i install_path [--prism] [--pmi_analysis] [--pyrmsd]

    Required arguments:
    -i install_path: The path to the directory where IMP and optional tools are installed. This should be the parent directory containing 'imp-clean', 'prism', 'PMI_analysis', and/or 'pyRMSD' directories.

    Optional arguments:
    --prism: Add PrISM to the environment variables.
    --pmi_analysis: Add PMI analysis to the environment variables.
    --pyrmsd: Add PyRMSD to the environment variables.
    -h: Show this help message and exit."
}

add_prism="no"
add_pmi_analysis="no"
add_pyrmsd="no"

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -i) install_path="$2"; shift ;;
        --prism) add_prism="yes" ;;
        --pmi_analysis) add_pmi_analysis="yes" ;;
        --pyrmsd) add_pyrmsd="yes" ;;
        -h) usage; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; usage; exit 1 ;;
    esac
    shift
done

if [ -z "$install_path" ]; then
    echo "Missing required arguments. Please provide a value for install_path."
    usage
    exit 1
fi

install_path=$(realpath $install_path) ;

if [[ -d "$install_path/imp-clean" ]]; then
    echo "IMP installation found at $install_path/imp-clean. Adding to PATH and environment variables..."
else
    echo "IMP installation not found at $install_path/imp-clean. Please check the path and try again."
    exit 1
fi

lines_to_add=(
    "export PATH=\$PATH:$install_path/imp-clean/build/bin"
    "export IMP_BIN_DIR=$install_path/imp-clean/build/bin"
    "export IMP_TMP_DIR=$install_path/imp-clean/build/tmp"
    "export IMP_DATA=$install_path/imp-clean/build/data"
    "export IMP_EXAMPLE_DATA=$install_path/imp-clean/build/doc/examples"
    "export SAMPCON_PATH=$install_path/imp-clean/build/lib/IMP/sampcon"
    "export PYTHONPATH=\$PYTHONPATH:$install_path/imp-clean/build/lib"
)

# add PrISM if the flag is set
if [[ "$add_prism" == "yes" ]]; then
    if [[ -d "$install_path/prism" ]] ; then
        echo "PrISM installation found at $install_path/prism. Adding to environment variables..."
    else
        echo "PrISM installation not found at $install_path/prism. Please check the path and try again."
        exit 1
    fi
    lines_to_add+=(
        "export PRISM_PATH=$install_path/prism"
    )

fi

# add PMI analysis if the flag is set
if [[ "$add_pmi_analysis" == "yes" ]]; then

    if [[ -d "$install_path/PMI_analysis" ]] ; then
        echo "PMI analysis installation found at $install_path/PMI_analysis. Adding to environment variables..."
    else
        echo "PMI analysis installation not found at $install_path/PMI_analysis. Please check the path and try again."
        exit 1
    fi

    lines_to_add+=(
        "export PMI_ANALYSIS_PATH=$install_path/PMI_analysis"
        "export PYTHONPATH=\$PYTHONPATH:$install_path/PMI_analysis/pyext/src"
    )

fi


# add PyRMSD if the flag is set
if [[ "$add_pyrmsd" == "yes" ]]; then

    if [[ -d "$install_path/pyRMSD" ]] ; then
        echo "PyRMSD installation found at $install_path/pyRMSD. Adding to environment variables..."
    else
        echo "PyRMSD installation not found at $install_path/pyRMSD. Please check the path and try again."
        exit 1
    fi

    lines_to_add+=(
        "export PYTHONPATH=\$PYTHONPATH:$install_path/pyRMSD"
    )

fi

for line in "${lines_to_add[@]}"; do
    grep -qxF "$line" ~/.bashrc || echo "$line" >> ~/.bashrc
done