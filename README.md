# msCRUSH 


## Introduction
msCRUSH (standing for **m**ass **s**pectrum **C**luste**R**ing **U**sing locality **S**ensitive **H**ashing) was developed by [Lei Wang](mailto:wang558@indiana.edu), [Sujun Li](https://scholar.google.com/citations?user=y4keCocAAAAJ&hl=en) and [Haixu Tang*](https://www.sice.indiana.edu/all-people/profile.html?profile_id=308), for the purpose of clustering large-scale tandem mass (MS/MS) spectra and then generating high quality *consensus spectra* for clusters of similar MS/MS spectra. Multithreading is enabled in this package.
msCRUSH can take as input multiple MGF files (regular expression is supported) of spectra of multiple charge states (including spectra without charge). If the input MS/MS spectra come with multiple charge states, then a clustering result file will be created for each charge state. 
## Prerequisites
g++ with version 5.1.0+ is required.

## Installation
1. `cd` to the main directory of msCRUSH.
2. Type `./install.sh`
3. Two executable files will be placed under *bin* directory: *mscrush_on_general_charge* for clustering similar spectra and *generate_consensus_spectrum_for_mscrush* for generating consensus spectra.

## Cluster similar spectra
- This process will cluster similar MS/MS spectra and then output to files the clusters containing the titles of spectra inside. If the input MS/MS spectra contain some spectra without charge, then these no-charge spectra will be clustered to spectra with charge (in the order of charge 2+, 3+, 4+, 5+, 1+, 6+, 7+, 8+, 9+, 10+ by default in msCRUSH) if their similarity is above threshold in each clustering iteration. 
  1. `cd bin`
  2. Usage: `./mscrush_on_general_charge <threads_to_use> <hash_func_num> <iteration> <min_similarity> <min_mz> <max_mz> <result_prefix> <mgf_file(s)>.`
  3. Typical example: `./mscrush_on_general_charge 40 15 100 0.65 200 2000 ../mgf/clusters ../mgf/D01*part*.mgf`. You will find 5 clustering result files with prefix *clusters* under directory `../mgf`.

- It is very important to select proper paramters for msCRUSH to cluster similar spectra with high sensitivity and specificity. Specifically:
  1. `hash_func_num` controls the collision probability of two spectra in a single hash table, the more hash functions, the smaller collision probabiilty, the more clusters, the faster msCRUSH program. 10 is a good starting point for datasets with size less than 1 million, while 15 is a good starting point for datasets with size around 10 million.
  2. `iteration` controls the number of hash tables to use with the aim to increase the probabiilty of two similar spectra to collide in at least one hash table, the more iterations, the higher collision probability, the less clusters, the slower mcCRUSH program. 100 iterations is a good starting point to play with.
  3. `min_similarity` is the minimum similarity to consider two spectra to be similar and merge into a consensus spectrum, the higher min similarity, the more clusters, the slower msCRUSH. 0.5 through 0.65 is recommended.
  4. `min_mz` is the minimum mz to consider a peak. 
  5. `max_mz` is the maximum mz to consider a peak.

## Generate consensus spectra
Note that writing consensus spectra (cs) to disk in MGF format can be time consuming, so if consensus spectra is needed, run scripts below.
1. `cd bin`
2. Usage: `./generate_consensus_spectrum_for_mscrush <cs_title_prefix> <cs_path_prefix> <mscrush_clusters_name(s)> <mgf_file(s)>.`
3. Typical example: `./generate_consensus_spectrum_for_mscrush D01CS ../mgf/D01 ../mgf/clusters-c*.txt ../mgf/D01*part*.mgf`. You will find 5 MGF files, each of which match to a clustering file of a specific charge state.

## Paper Pointer
[Link](https://www.biorxiv.org/content/early/2018/05/11/308627)
<!---
## Citation
```latex
    @article {Wang308627,
      author = {Wang, Lei and Li, Sujun and Tang, Haixu},
      title = {msCRUSH: fast tandem mass spectra clustering using locality sensitive hashing},
      year = {2018},
      doi = {10.1101/308627},
      publisher = {Cold Spring Harbor Laboratory},
      URL = {https://www.biorxiv.org/content/early/2018/04/25/308627},
      eprint = {https://www.biorxiv.org/content/early/2018/04/25/308627.full.pdf},
      journal = {bioRxiv}
    }
```
-->
## Questions
Please contact Lei Wang (wang558@indiana.edu) for assistance.
## Acknowledgement
This work was supported by the NIH grant 1R01AI108888 and the Indiana University Precision Health Initiative.
