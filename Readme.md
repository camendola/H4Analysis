H4Analysis
==========
[<img src="https://gitlab.cern.ch/spigazzi/H4Analysis/badges/master/pipeline.svg">](https://gitlab.cern.ch/spigazzi/H4Analysis/-/pipelines/latest)

# [Documentation](https://h4analysis.web.cern.ch/)
  - This repository aims to provide a fast reconstruction of data
    acquired with the H4DAQ.
  - The main executable is =bin/H4Reco=. The data processing is 
    steered through the input cfg file (examples in cfg/).

# Requirements
## Code:
  - c++17
  - ROOT
## Docs:
  - python=3.9
  - sphinx
  - breathe
  - sphinx\_rtd\_theme
  - sphinxcontrib-srclinks
## Conda environment:
  The above requirements can be easily satisfied by creating a dedicated conda environment:
  
  conda env create -f environment.yml

## LXPLUS

  source script/setup.sh(csh) (only on lxplus, or any machine with afs and EOS access). 

# Install and run          
   ```
   git clone --recursive https://github.com/simonepigazzini/H4Analysis.git
   cd H4Analysis
   make
   ```
   to run the reconstruction locally
   ```
   bin/H4Reco cfg/Oct2015_timing.cfg 4443
   ```
#  to lounch reco on grid
   ```
   python scripts/submitBatch.py -r [run_number] -c cfg/ECAL_H4_Oct2021/ECAL_H4_Phase2_base_100GeV_HP.cfg -v [production_tag] -s [output_directory] --notar --spills-per-job 1 -q tomorrow
   ```
   find examples in `scripts/submit2021_HPfit_Escan.sh`
     
