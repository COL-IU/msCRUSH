# SPLASH 


## Introduction
SPLASH (standing for SPectrum clustering using LocAlity Sensitive Hashing) was developed by [Lei Wang](mailto:wang558@indiana.edu), [Sujun Li](https://scholar.google.com/citations?user=y4keCocAAAAJ&hl=en) and [Haixu Tang*](https://www.sice.indiana.edu/all-people/profile.html?profile_id=308), for the purpose of clustering large-scale tandem mass spectra (MS/MS spectra) and then generating high quality *consensus spectra* for clusters of similar MS/MS spectra. Multithreading is enabled in this package.
SPLASH can take as input multiple MGF files (regular expression is supported) of spectra of multiple charge states (including spectra without charge). If the input MS/MS spectra come with multiple charge states, then a clustering result file will be created for each charge state. 
## Prerequisites
g++ with version 5.1.0+ is required.

## Installation
1. `cd` to the main directory of SPLASH.
2. Type `./install.sh`
3. Two executable files will be placed under *bin* directory: *splash_on_general_charge* for clustering similar spectra and *generate_consensus_spectrum_for_splash* for generating consensus spectra.

## Cluster similar spectra
This process will cluster similar MS/MS spectra and then output to files the clusters containing the titles of spectra inside. If the input MS/MS spectra contain some spectra without charge, then these no-charge spectra will be clutered to spectra with charge (in the order of charge 2+, 3+, 4+, 5+, 1+, 6+, 7+, 8+, 9+, 10+ by default in SPLASH) if their similarity is above threshold in each clustering iteration. 
1. `cd bin`
2. Usage: `./splash_on_general_charge threads_to_use hash_func_num iteration min_similarity result_prefix mgf_file(s).`
3. Typical example: `./splash_on_general_charge 40 15 100 0.65 ../mgf/clusters ../mgf/D0*.mgf`. You will find 5 clustering result files with prefix *clusters* under directory `../mgf`.

## Generate consensus spectra
Note that writing consensus spectra (cs) to disk in MGF format can be time consuming, so if consensus spectra is needed, run scripts below.
1. `cd bin`
2. Usage: `./generate_consensus_spectrum_for_splash cs_title_prefix cs_path_prefix splash_clusters_name(s) mgf_file(s).`
3. Typical example: `./generate_consensus_spectrum_for_splash D01CS ../mgf/D01 ../mgf/clusters-c*.txt ../mgf/D01*part*.mgf`. You will find 5 MGF files, each of which match to a clustering file of a specific charge state.

## Citation
```latex
  @article {Wang308627,
    author = {Wang, Lei and Li, Sujun and Tang, Haixu},
    title = {SPLASH: fast tandem mass spectra clustering using locality sensitive hashing},
    year = {2018},
    doi = {10.1101/308627},
    publisher = {Cold Spring Harbor Laboratory},
    URL = {https://www.biorxiv.org/content/early/2018/04/25/308627},
    eprint = {https://www.biorxiv.org/content/early/2018/04/25/308627.full.pdf},
    journal = {bioRxiv}
  }

```

## Questions
Please contact Lei Wang (wang558@indiana.edu) for assistance.
## Acknowledgement
This work was supported by the NIH grant 1R01AI108888 and the Indiana University Precision Health Initiative.
