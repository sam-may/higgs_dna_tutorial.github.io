# HiggsDNA Tutorial
Welcome to the tutorial for the HiggsDNA (Higgs to Diphoton NanoAOD) framework! The [HiggsDNA](https://gitlab.cern.ch/HiggsDNA-project/HiggsDNA) repository provides tools for running a Higgs to diphoton analysis.

# 1. Introduction
TODO

# 2. Setup
The installation procedure consists in the following steps:

**1. Clone this repository**
```  
git clone --recursive https://gitlab.cern.ch/HiggsDNA-project/HiggsDNA
cd HiggsDNA
```
**2. Install dependencies**

The necessary dependencies (listed in ```environment.yml```) can be installed manually, but the suggested way is to create a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/mana
ge-environments.html) by running:
``` 
conda env create -f environment.yml
```
the conda env can become pretty large (multiple GB), so you may want to specify an installation location in the above step with
```
conda env create -f environment.yml -p <path to install conda env>
```
then activate the environment with
```
conda activate higgs-dna
```
If you specified an installation path, your activation command will instead be:
```
conda activate <path to install conda env>/higgs-dna
```

** lxplus-specific notes **
If you are running on `lxplus` you may run into permission errors, which can be fixed with:
```
chmod 755 -R /afs/cern.ch/user/<your_username_first_initial>/<your_username>/.conda
```
You may also want to increase your disk quota at [this link](https://resources.web.cern.ch/resources/Manage/EOS/Default.aspx), otherwise you may run out of space while installing your `conda` environment.

Please note that the field ```python>=3.6``` will create an environment with the most recent stable version of Python. Change it to suite your needs (but still matching the requirement of Python>=3.6).

One additional package, `correctionlib`, must be installed via `pip`, rather than `conda`. Run
```
setup.sh
```
to install this script.

**3. Install ```higgs_dna```**

**Users** can install the package by simply running:
``` 
python setup.py install
```
(when a stable version will be available, it will be uploaded to the PyPI and it will be possible to install with just ```pip install higgs_dna``` without the need to clone the repository).


For **developers**, the suggested way to install is:
``` 
pip install -e .
```
this prevents the need to run the installation step every time a change is performed.

If you notice issues with the ```conda pack``` command for creating the tarball, try updating and cleaning your environment with (after running ```conda activate higgs-dna```):
```
conda env update --file environment.yml --prune
``` 

# 3. Using HiggsDNA for physics analysis 
The likely starting point for the majority of users is the `scripts/run_analysis.py` script. This script can be used to run a selection (i.e. a sequence of `higgs_dna.taggers.tagger.Tagger` objects) and apply scale factors and corrections (i.e. a set of `higgs_dna.systematics.systematic.Systematic` objects) over a list of samples (specified through a `json` file), and create a set of ntuple-like outputs (in the `parquet` format) with a specified set of fields (or ``branches''), merging the outputs and calculating scale1fb and scaling the normalization of weights for MC samples, if requested.

We can illustrate most of the functionality through an example: suppose I want to develop a ttH analysis.  

