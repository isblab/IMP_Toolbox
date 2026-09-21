# IMP_Toolbox

> [!NOTE]
> If you can not find any of the scripts you have previously used from IMP_Toolbox in the `main` branch
> , please check in the [`archive`](https://github.com/isblab/IMP_Toolbox/tree/archive) branch.

- A toolbox to aid integrative modeling with [IMP](https://github.com/salilab/imp).

- See [Modeling with IMP](https://docs.google.com/document/d/1gaG83RsEBQNemuWwhra0TP0c0jg_WwdeZI1yD_S3RPQ/edit?usp=sharing) for helpful tips while starting.

## Pre-requisites

- Conda or Miniconda (https://www.anaconda.com/docs/getting-started/miniconda/install/linux-install)
- Fedora OS

## Installation

- To use IMP_Toolbox, clone this repository
  ```bash
  git clone https://github.com/isblab/IMP_Toolbox.git
  ```

- Add the path to `~/.bash_profile`. You can then use it as a module.
  ```bash
  # replace the path and add it at the end of ~/.bash_profile
  export PYTHONPATH=/path/to/IMP_Toolbox:$PYTHONPATH
  ```

- If you want to install IMP, follow the instructions given in [IMP installation](#imp-installation) and activate the conda environment.
  ```bash
  conda activate imp_omg
  ```

- Alternatively, if you only want to use [IMP_Toolbox](./) without installing IMP, create a new conda environment. You need `Python >=3.12`.
  ```bash
  conda create -n imp_toolbox python=3.12
  conda activate imp_toolbox
  ```

- Install the required packages from [`requirements.txt`](./requirements.txt).
  ```bash
  pip install -r requirements.txt
  ```

### Optional

- If you plan to use [`sequence`](./IMP_Toolbox/sequence/) module, you need to install `EMBOSS` package as follows:
  ```bash
  sudo dnf install EMBOSS
  ```

- If you plan to use [`structure.burial`](./IMP_Toolbox/structure/burial.py), you need
  to download [DSSP](https://pdb-redo.eu/dssp/download).
  ```bash
  sudo dnf install dssp
  ```

- If you plan to use [`structure.burial.get_burial_info`](./IMP_Toolbox/structure/burial.py)
  with `include_residue_depth` set to True, you also need [MSMS](https://ccsb.scripps.edu/msms/downloads/) installed.

## IMP Installation

- Go to [`IMP_installation`](./IMP_Toolbox/IMP_installation/) directory

- Run the following command in terminal. (change the installation path and conda environment)
  ```bash
  bash install_imp_conda.sh \
    -e imp_omg \
    -i ~/IMP_OMG \
    -s tarball \
    -I 2.24.0 \
    -c 16 \
    -P 3.12
  ```

- Use `-h` flag to see more details.

## Additional Information

**License:** [GPLv3](./LICENSE)

**Testable:** Yes