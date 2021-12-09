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
## 3.1 Overview
The likely starting point for the majority of users is the `scripts/run_analysis.py` script. This script can be used to run a selection (i.e. a sequence of `higgs_dna.taggers.tagger.Tagger` objects) and apply scale factors and corrections (i.e. a set of `higgs_dna.systematics.systematic.Systematic` objects) over a list of samples (specified through a `json` file), and create a set of ntuple-like outputs (in the `parquet` format) with a specified set of fields (or ``branches''), merging the outputs and calculating scale1fb and scaling the normalization of weights for MC samples, if requested.

The `scripts/run_analysis.py` script can be used for a variety of purposes, including:
- **Developing an analysis** : creating ntuple-like files to form the starting point for developing an analysis (making data/MC plots, yield tables, training MVAs, etc).
- **Systematics studies** : creating ntuple-like files to form the starting point for performing systematics studies (e.g. deriving a correction/scale factor for photons in a Z->ee selection)
- **Running a full analysis** : running a selection, possibly with multiple tags, and propagating relevant systematicsand their associated uncertainties.

### 3.1.1 Creating a config file
The script can be configured to perform any of these tasks through a single `json` config file, where there are 5 main fields the user will want to pay attention to:
1. `"tag_sequence"` : a list which specifies a set of `higgs_dna.taggers.tagger.Tagger` objects and dictates the event selection.
2. `"systematics"` : a dictionary with two keys, `"weights"` and `"independent_collections"`, under which a dictionary is specified for each `higgs_dna.systematics.systematic.WeightSystematic` and `higgs_dna.systematics.systematic.SystematicWithIndependentCollection`, respectively.
3. `"samples"` : a dictionary which points to a catalog of samples, and specifies which samples and years should be processed from this catalog.
4. `"variables_of_interest"` : a list of variables which should be saved in the output `parquet` files
5. `"branches"` : a list of branches which should be read from the input nanoAODs.

There are two additional fields, which likely do not need to be changed by the user, but should be present:
6. `"name"` : a string which can identify this analysis
7. `"function"` : specifies which function should be run for each job. At the time of writing, this is always the `higgs_dna.analysis.run_analysis` function, but is left to be configurable for potential future needs.

A basic config file to run the diphoton preselection without systematics would look like this:
```json
{
    "name" : "tth_preselection",
    "function" : {
        "module_name" : "higgs_dna.analysis",
        "function_name" : "run_analysis"
    },
    "tag_sequence" : [
            {
                "module_name" : "higgs_dna.taggers.diphoton_tagger",
                "tagger" : "DiphotonTagger"
            }
    ],
    "systematics" : {
        "weights" : {},
        "independent_collections" : {},
    },
    "samples" : {
        "catalog" : "metadata/samples/tth_tutorial.json",
	"sample_list" : ["Data", "ttH_M125", "ggH_M125"],
	"years" : ["2016", "2017", "2018"]
    },
    "variables_of_interest" : [
        ["Diphoton", "mass"], ["Diphoton", "pt"], ["Diphoton", "eta"], ["Diphoton", "phi"], ["Diphoton", "dR"],
        ["LeadPhoton", "pt"], ["LeadPhoton", "eta"], ["LeadPhoton", "phi"], ["LeadPhoton", "mass"], ["LeadPhoton", "mvaID"], ["LeadPhoton", "genPartFlav"], ["LeadPhoton", "pixelSeed"],
        ["SubleadPhoton", "pt"], ["SubleadPhoton", "eta"], ["SubleadPhoton", "phi"], ["SubleadPhoton", "mass"], ["SubleadPhoton", "mvaID"], ["SubleadPhoton", "genPartFlav"], ["SubleadPhoton", "pixelSeed"],
        "GenHggHiggs_pt", "GenHggHiggs_eta", "GenHggHiggs_phi", "GenHggHiggs_mass", "GenHggHiggs_dR",
        "weight_central",
        "event"
    ],
    "branches" : [
            "Photon_pt", "Photon_eta", "Photon_phi", "Photon_mass",
            "Photon_pixelSeed", "Photon_mvaID", "Photon_electronVeto",
            "Photon_sieie",
            "Photon_r9", "Photon_hoe", "Photon_pfRelIso03_chg", "Photon_pfRelIso03_all",
            "Photon_isScEtaEB", "Photon_isScEtaEE",
            "Photon_trkSumPtHollowConeDR03", "Photon_photonIso", "Photon_chargedHadronIso",
            "Photon_genPartFlav",
            "genWeight", "run", "event", "fixedGridRhoAll"
    ]
}
```

## 3.2 Example : ttH analysis
We can illustrate most of the functionality through an example: suppose I want to develop a ttH analysis.
As a first step, we will run the standard diphoton preselection and a loose ttH preselection to produce a `parquet` file which we can use to make data/MC plots and yield tables and train MVAs.

### 3.2.1 Constructing a tag sequence
As an initial step, we can construct a `json` file which will perform the diphoton preselection.
To start, 

### 3.2.2 Adding systematics

### 3.2.3 Defining a list of samples

### 3.2.4 Assessing the outputs
