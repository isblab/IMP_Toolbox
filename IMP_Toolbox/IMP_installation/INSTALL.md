# Installing IMP and related packages

- Install anaconda or miniconda. You can follow [these](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) instructions.

- The following steps are meant for building IMP from source code, which is necessary if you
  want to use a custom module (e.g. module for MPDBR restraint) inside IMP. If you do not have
  such a module and want to use standard IMP, the easiest way to install it is via one of the
  packaged versions from official website of IMP (https://integrativemodeling.org).
    - <img src="https://integrativemodeling.org/images/linux.svg" width="20"/> [Packages for Linux distributions](https://integrativemodeling.org/download-windows.html)
    - <img src="https://integrativemodeling.org/images/apple.svg" width="20"/> [Package for Mac](https://integrativemodeling.org/download-windows.html)
    - <img src="https://integrativemodeling.org/images/windows.svg" width="20"/> [Binary for Windows](https://integrativemodeling.org/download-windows.html)
    - <img src="https://integrativemodeling.org/images/anaconda-symbol.svg" width="20"> [Anaconda package](https://integrativemodeling.org/download-anaconda.html)

- Alternatively, follow the steps mentioned here except for the custom module installation part.

> [!NOTE]
> The following steps are only meant for Fedora Linux distribution. If you want
> to build IMP from source code on other OS, follow the [installation guide](https://integrativemodeling.org/2.24.0/doc/manual/installation.html).

- Use the following command in terminal.

    ```bash
    install_imp_conda.sh \
    -e conda_env_name \
    -i install_path \
    -s install_mode \
    -I imp_version \
    -c num_cores \
    -P python_version \
    --prism \
    --pmi_analysis \
    --pyrmsd
    ```

- A conda environment of the provided name will be created if it does not already exists.
  All the python packages required for modeling with IMP will be installed in this environemnt.

- `install_path` is the path to the directory where you would like to install IMP.

- `install_mode` specifies two options:
    - `github`: Installing IMP from the [github repository](https://github.com/salilab/imp).
    - `tarball`: Install from tarball specified in [releases](https://github.com/salilab/imp/releases).

- `imp_version` specifies the version of IMP to install. Check the
  [website](https://integrativemodeling.org/download.html) or
  [releases](https://github.com/salilab/imp/releases) to see
  the most recent stable version of IMP.

- `num_cores` specifies the number of cores to be used to build IMP.

- `python_version` specifies the python version to install in the specified conda environemnt.

- `--prism`, `--pmi_analysis`, and `--pyrmsd` flags are used to clone the respective repositories for
  [PrISM](https://github.com/isblab/prism), [PMI_analysis](https://github.com/salilab/PMI_analysis),
  and [pyRMSD](https://github.com/salilab/pyRMSD).

- If you want to install a custom module use the following flag and provide the path
  to the module.
    ```
    -m /path/to/custom_module
    ```

- There are two additional optional flags.
    - `-y`: To skip confirmation prompts and use default values for environment name and installation choices.
    - `-h`: To display help message.

> [!TIP]
> If using "github" mode, either the tag (tag ~ IMP release version) or branch name
> can also be used for `imp_version`.

> [!NOTE]
> There is an alternate script which does not require conda. But, it does not include all the options
> as above. Use it at your own discretion.
>
>  ```bash
>  bash install_imp.sh </path/to/imp-install/> <num-cores> <mode> <imp_version>>
>  # e.g,
>  bash install_imp.sh /path/$USER/IMP_OMG 16 tarball 2.24.0
>  # `mode` can be "github" or "tarball".
>  # If custom module is being used, you will prompted to provide it's path.
>  # Currently only supports installing one custom module.
>  ```
