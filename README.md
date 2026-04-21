# IMP_Toolbox

> [!NOTE]
> If you can not find any of the scripts you have previously used from IMP_Toolbox in the `main` branch
> , please check in the [`archive`](https://github.com/isblab/IMP_Toolbox/tree/archive) branch.

- See [Modeling with IMP](https://docs.google.com/document/d/1gaG83RsEBQNemuWwhra0TP0c0jg_WwdeZI1yD_S3RPQ/edit?usp=sharing) for helpful tips while starting.

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

- Install the required packages from `requirements.txt`. You need `Python >=3.12`.
  ```bash
  conda install --file requirements.txt
  ```

- If you plan to use `sequence` module, you need to install `EMBOSS` package as follows:
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

- Use the following command in terminal. (change the installation path and conda environment)
  ```bash
  bash install_imp_conda.sh \\
    -e imp_omg \\
    -i ~/IMP_OMG \\
    -s tarball \\
    -I 2.24.0 \\
    -c 16 \\
    -P 3.12
  ```

- Use `-h` flag to see more details.

## Additional Information

**License:** [GPLv3](./LICENSE)

**Testable:** Yes