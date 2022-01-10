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

It also prints out some useful information about the physics content and technical aspects of running:
```
[AnalysisManager : run] Task 'ttH_M125_2017', PERFORMANCE summary                                             
DEBUG             [PERFORMANCE : ttH_M125_2017] Processed 414189 total events in 0:01:12.598349         analysis.py:306
         (hours:minutes:seconds) of total runtime (5705.21 Hz).                                                        
DEBUG             [PERFORMANCE : ttH_M125_2017] Fraction of runtime spent on load: 91.28 percent        analysis.py:308
DEBUG             [PERFORMANCE : ttH_M125_2017] Fraction of runtime spent on syst: 0.00 percent         analysis.py:308
DEBUG             [PERFORMANCE : ttH_M125_2017] Fraction of runtime spent on taggers: 7.48 percent      analysis.py:308
DEBUG    [AnalysisManager : run] Task 'ttH_M125_2017', PHYSICS summary                                  analysis.py:310
DEBUG             [PHYSICS : ttH_M125_2017] events set 'nominal' has eff. of 49.59 percent (all         analysis.py:312
         taggers)                                                                                                      
DEBUG             [PHYSICS : ttH_M125_2017] With a cross section times BF of 0.001151117 pb and a sum   analysis.py:314
         of weights of 409277.240008000, scale1fb for this sample is 0.000002813
```

We can then explore this output file in several ways. One useful tool is located under `higgs_dna/bonus/assess.py` and will be described later in the tutorial, but for quick exploration, the `python` interpreter serves nicely:
```
python
Python 3.9.7 | packaged by conda-forge | (default, Sep 29 2021, 19:20:46) 
[GCC 9.4.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import awkward
>>> events = awkward.from_parquet("tutorial_dipho/merged_nominal.parquet")
```

We can check the total number of raw events and the fields which are saved for each event:
```
>>> len(events)
539177
>>> events.fields
['Diphoton_mass', 'Diphoton_pt', 'Diphoton_eta', 'Diphoton_phi', 'Diphoton_dR', 'LeadPhoton_pt', 'LeadPhoton_eta', 'LeadPhoton_phi', 'LeadPhoton_mass', 'LeadPhoton_mvaID', 'LeadPhoton_genPartFlav', 'LeadPhoton_pixelSeed', 'SubleadPhoton_pt', 'SubleadPhoton_eta', 'SubleadPhoton_phi', 'SubleadPhoton_mass', 'SubleadPhoton_mvaID', 'SubleadPhoton_genPartFlav', 'SubleadPhoton_pixelSeed', 'GenHggHiggs_pt', 'GenHggHiggs_eta', 'GenHggHiggs_phi', 'GenHggHiggs_mass', 'GenHggHiggs_dR', 'event', 'weight_central_initial', 'weight_central', 'process_id', 'year']
```

We can check the total number of weighted events for all Run 2 and for an individual year:
```
>>> awkward.sum(events.weight_central)
77.62231871412764
>>> awkward.sum(events.weight_central[events.year == 2018])
33.89958685493177
```

We can check the diphoton pT and compare it to the gen-level pT:
```
>>> events.Diphoton_pt
<Array [91.8, 106, 110, 271, ... 13, 67.9, 230] type='539177 * float32'>
>>> events.GenHggHiggs_pt
<Array [91.2, 108, 109, 270, ... 13.4, 69, 224] type='539177 * float64'>
```

Next, we can add a loose ttH preselection that one might use as a starting point for developing a ttH analysis: data/MC studies, training MVAs, etc.

For a loose ttH preselection we will define two orthogonal channels, hadronic and leptonic, targeting fully hadronic and semi/di-leptonic decays of the ttbar system, respectively.

The hadronic channel will require:
- exactly 0 leptons
- at least 4 jets

while the leptonic channel will require:
- at least 1 lepton
- at least 2 jets 

This is implemented in a `higgs_dna.taggers.tagger.Tagger` class under `higgs_dna/taggers/tutorial.py`. We will walk through this example.

The first (optional) step is creating a constructor for the tagger. This is not necessary, but in the example below we define the constructor such that the default options can be updated through the config `json`. In this way, to update a single cut value, for example the electron pT cut, this single value can be included in the config `json` and that value will be updated while leaving all other values as their default values:
```
def __init__(self, name, options = {}, is_data = None, year = None):
        super(TTHPreselTagger, self).__init__(name, options, is_data, year)

        if not options:
            self.options = DEFAULT_OPTIONS
        else:
            self.options = misc_utils.update_dict(
                    original = DEFAULT_OPTIONS,
                    new = options
            )
```

Next, we can define the function which actually performs the selection, `calculate_selection`. We will first select our physics objects: electrons, muons, and jets. The helper functions `select_electrons`, `select_muons`, and `select_jets` located under `higgs_dna/selections/lepton_selections.py` and `higgs_dna/selections/jet_selections.py` will be useful for this.

Each of these functions has default options which can be seen in their respective files. For example, for electrons, we see:
```json
DEFAULT_ELECTRONS = {
        "pt" : 10.0,
        "eta" : 2.4,
        "dxy" : 0.045,
        "dz" : 0.2,
        "id" : "WP90", # this means that the 90% eff working point ID will be used
        "dr_photons" : 0.2, # requires that all selected electrons be at least dR > 0.2 from selected photons
        "veto_transition" : True
}
```
Perhaps we would like to deviate from the default options, requiring higher pT for our electrons and muons. We can do this by setting some default options for our tagger (this could also be done inside the config `json`):
```json
DEFAULT_OPTIONS = {
    "electrons" : {
        "pt" : 20.0,
        "dr_photons" : 0.2
    },
    "muons" : {
        "pt" : 20.0,
        "dr_photons" : 0.2
    }
}   
```
In doing this, the pT cuts for the electrons and muons will be 20.0, while all other cuts will remain the same as the default options listed under `leptons_selections.py`. We also define the option `"dr_photons"` so that we can ensure the selected leptons are not overlapping with the selected photons.

We first calculate the cut which gives us the selected electrons from the input Electron collection:
```python
electron_cut = lepton_selections.select_electrons(
	electrons = syst_events.Electron,
	options = self.options["electrons"],
	clean = {
	    "photons" : {
		"objects" : syst_events.Diphoton.Photon,
		"min_dr" : self.options["electrons"]["dr_photons"]
	    }
	},
	name = "SelectedElectron",
	tagger = self
)
```
and then add these electrons to the events array with the name "SelectedElectron":
```python
electrons = awkward_utils.add_field(
	events = syst_events,
	name = "SelectedElectron",
	data = syst_events.Electron[electron_cut]
)
it is important to add the selected electrons to the events array if other Taggers or Systematics will need to access these (for example, if we wanted to calculate apply an electron scale factor, this would likely be applied on SelectedElectron, rather than the nominal Electron collection in the nanoAODs). We can then do the same for the muons and the jets.

We might also want to sort the jets by b-tag score (rather than pT by default) in order to do things like count the number of jets passing a given working point, or getting the highest b-tag score in an event. We can obtain a resorted version of the jets with:
```
bjets = jets[awkward.argsort(jets.btagDeepFlavB, axis = 1, ascending = False)]
```

we can then finish our preselection by making requirements on the number of leptons and jets:
```python
# Preselection
n_electrons = awkward.num(electrons)
n_muons = awkward.num(muons)
n_leptons = n_electrons + n_muons

n_jets = awkward.num(jets)

# Hadronic presel
hadronic = (n_leptons == 0) & (n_jets >= 4)

# Leptonic presel
leptonic = (n_leptons >= 1) & (n_jets >= 2)

presel_cut = hadronic | leptonic

return presel_cut, syst_events
```



### 3.2.2 Adding systematics

### 3.2.3 Defining a list of samples

### 3.2.4 Assessing the outputs

# awkward Arrays and Columnar Analysis
