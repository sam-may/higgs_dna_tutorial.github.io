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
    "name" : "diphoton_preselection",
    "function" : {
        "module_name" : "higgs_dna.analysis",
        "function_name" : "run_analysis"
    },
    "tag_sequence" : [
            {
                "module_name" : "higgs_dna.taggers.diphoton_tagger",
                "tagger" : "DiphotonTagger",
		"kwargs" : {
	            "name" : "diphoton_tagger",
		    "options" : {
			"gen_info" : { "calculate" : true }
		    }
		}
            }
    ],
    "systematics" : {
        "weights" : {},
        "independent_collections" : {}
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
 	    "GenPart_eta", "GenPart_genPartIdxMother", "GenPart_mass", "GenPart_pdgId", "GenPart_phi", "GenPart_pt", "GenPart_status","GenPart_statusFlags",
            "genWeight", "run", "event", "fixedGridRhoAll"
    ]
}
```

## 3.2 Example : ttH analysis
We can illustrate most of the functionality through an example: suppose I want to develop a ttH analysis.
As a first step, we will run the standard diphoton preselection and a loose ttH preselection to produce a `parquet` file which we can use to make data/MC plots and yield tables and train MVAs.

### 3.2.1 Constructing a tag sequence
As an initial step, we can construct a `json` file which will perform the diphoton preselection. We can start with the the file above, which is also located in `higgs_dna/metadata/tutorial/diphoton.json`.

In this example, the `DiphotonTagger` would be run over all events (with the default options listed in `higgs_dna/taggers/diphoton_tagger.py`), saving the variables listed under `"variables_of_interest"` into parquet files. All samples listed under `"samples"` will be processed for all years listed under `"years"`. To start, we can run over just the ttH sample and run over just a couple files:
```
python run_analysis.py
--log-level "DEBUG"
--config "metadata/tutorial/diphoton.json" 
--merge_outputs # merge all parquets into 1 file
--short # run just one job per sample/year
--sample_list "ttH_M125" # run over only these samples (csv)
--output_dir "tutorial_dipho" # directory to save outputs to
```

After verifying that there are no errors, we can remove the `--short` option and run over the full set of Run 2 ttH M125 MC.

Since the samples are accessed via `xrootd`, HiggsDNA will check for a valid proxy first and raise an error if it does not find one:
```
WARNING  [misc_utils : check_proxy] We were not able to find grid proxy or proxy was found to be      misc_utils.py:127
         expired. Output of 'voms-proxy-info': ['', "Couldn't find a valid proxy.", '']                                
ERROR    [CondorManager : prepare_inputs] We were not able to find grid proxy or proxy was found  sample_manager.py:114
         to be expired. Since you are accessing files through xrd, a valid proxy is necessary.                         
         NoneType: None 
```

After ensuring we have a valid proxy, HiggsDNA will first print out some summary info about the run options we selected with our config file. Next, it will run the `SampleManager` class to loop through each of the samples and years and grab the relevant files:
```
DEBUG    [SampleManager : get_samples] For sample 'ttH_M125', year '2018', found xs of 0.507100 pb sample_manager.py:73
         and bf of 0.002270                                                                                            
DEBUG    [SampleManager : get_samples] For sample 'ttH_M125', year '2018', we interpreted the     sample_manager.py:109
         specified files 'root://redirector.t2.ucsd.edu//store/user/hmei/nanoaod_runII/HHggtautau                      
         /ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_RunIIAutumn18MiniAOD-102X_upgrade201                      
         8_realistic_v15-v1_MINIAODSIM_v0.6_20201021/' as a directory to be accessed via xrd                           
DEBUG    [SampleManager : get_files_from_xrd] We will find files for dir 'root://redirector.t2.uc sample_manager.py:270
         sd.edu//store/user/hmei/nanoaod_runII/HHggtautau/ttHJetToGG_M125_13TeV_amcatnloFXFX_mads                      
         pin_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_v0.6_20201                      
         021/' with the command:                                                                                       
          xrdfs root://redirector.t2.ucsd.edu/ ls /store/user/hmei/nanoaod_runII/HHggtautau/ttHJe                      
         tToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_real                      
         istic_v15-v1_MINIAODSIM_v0.6_20201021/                                                                        
DEBUG    [SampleManager : get_files_from_xrd] Found 17 files in dir 'root://redirector.t2.ucsd.ed sample_manager.py:277
         u//store/user/hmei/nanoaod_runII/HHggtautau/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_p                      
         ythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_v0.6_20201021/' 
```

The `Manager` classes in `higgs_dna/job_management/managers` will then split the files into jobs which may either be submitted locally (default) or to `HTCondor`:
```
INFO     [Task: create_jobs] Task 'ttH_M125_2016' : splitting 10 input files into 4 jobs                     task.py:63
100%|███████████████████████████████████████████████████████████████████████████████████| 4/4 [00:00<00:00, 137.10it/s]
INFO     [Task: create_jobs] Task 'ttH_M125_2017' : splitting 10 input files into 4 jobs                     task.py:63
100%|███████████████████████████████████████████████████████████████████████████████████| 4/4 [00:00<00:00, 139.74it/s]
INFO     [Task: create_jobs] Task 'ttH_M125_2018' : splitting 17 input files into 6 jobs                     task.py:63
100%|███████████████████████████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 218.09it/s]
INFO     [AnalysisManager : run] Running 3 tasks with 37 total input files split over 3 total jobs.     analysis.py:283
```

NOTE: the current job management tools are planned to be superseded by `coffea` style tools.

The `Manager` will then periodically print out the progress of the jobs:
```
[JobsManager : summarize] Summarizing task progress.                                                          

INFO     [JobsManager : summarize] Task 'ttH_M125_2016' :                                                managers.py:81
                 1/1 (100.00 percent) jobs running                                                                     
INFO     [JobsManager : summarize] Task 'ttH_M125_2017' :                                                managers.py:81
                 1/1 (100.00 percent) jobs running                                                                                                                                                                                       
INFO     [JobsManager : summarize] Task 'ttH_M125_2018' :                                                managers.py:81
                 1/1 (100.00 percent) jobs completed
``` 
until all jobs are completed. Once all jobs are completed, it will then merge each of the individual `parquet` files. In the cases of merging files of different years or samples, new fields `year` and `process_id` will be added to the `parquet` files with relevant values. Merging will also calculate scale1fb values for MC samples (provided a cross section is given) based on the number of successfully processed events. This means that in the case of a single corrupt MC file that does not run successfully, your normalization will still be calculated correctly. The `weight_central` branch will be scaled according to the scale1fb and luminosity for the given year (though the unscaled values are still available in the unmerged files, if desired).

HiggsDNA lets us know that all jobs have completed and tells us where the merged `parquet` file is located:
```
INFO     [JobsManager : merge_outputs] For syst_tag 'nominal' merging 3 files into merged file          managers.py:126
         '/home/users/smay/HiggsDNA/scripts/tutorial_dipho/merged_nominal.parquet'.                                    
DEBUG             merging 539177 events                                                                 managers.py:130
INFO     [AnalysisManager : run] Finished running analysis 'diphoton_preselection'. Elapsed time:       analysis.py:302
         0:03:51.303396 (hours:minutes:seconds). 
```


 

### 3.2.2 Adding systematics

### 3.2.3 Defining a list of samples

### 3.2.4 Assessing the outputs

# awkward Arrays and Columnar Analysis
