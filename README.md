# IMP_Toolbox

- See [Modeling with IMP](https://docs.google.com/document/d/1gaG83RsEBQNemuWwhra0TP0c0jg_WwdeZI1yD_S3RPQ/edit?usp=sharing) for helpful tips while starting.
- To use IMP_Toolbox, clone this repository and add the path to `~/.bash_profile`. You can then use its subdirectories as a module (e.g. `from af_pipeline import AFInput`).
```bash
# replace the path and add it at the end of ~/.bash_profile
export PYTHONPATH=/path/to/IMP_Toolbox:$PYTHONPATH
```
- Install the required packages from `requirements.txt`. You need `Python >=3.12`.
```bash
conda install --file requirements.txt
```

**May the light of Durin's day shine upon those who choose this path.**

## Installation

- Use the following command in terminal.\
 (`bash install_imp.sh </path/to/imp-install/> <num-cores> <mode> <imp_version>`)
```bash
bash install_imp.sh /path/$USER/IMP_OMG 16 tarball 2.23.0

# or if using conda
bash install_imp_conda.sh imp_omg /path/$USER/IMP_OMG 16 tarball 2.23.0
```

- `mode` can be "github" or "tarball".

- If using "github" mode, either the tag (tag ~ IMP release version) or branch name can be used for `imp_version`.

- If custom module is being used, you will prompted to provide it's path. Currently only supports installing one custom module.

## Pre-processing

Note: the classes for pre-processing and AFpipeline are in `preprocessing` and `af_pipeline`. You need to import and use them for your own systems. See `examples` for example scripts. 


## Modeling (not what you may think at first)

## Analysis

## Validation

## License
[CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0
International License.  
