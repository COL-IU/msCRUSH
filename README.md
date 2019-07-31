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
- This process will cluster similar MS/MS spectra and then output to files the clusters containing the titles of spectra inside. If the input MS/MS spectra contain some spectra without charge, then these no-charge spectra will be clustered to spectra with charge (in the order of charge 2+, 3+, 4+, 5+, 1+, 6+, 7+, 8+, 9+, 10+ by default in msCRUSH) if their similarity is above preset clustering threshold in each respective clustering iteration. 
  1. `cd bin`
  2. Usage: `./mscrush_on_general_charge -f mgf_files(s) [-t threads_to_use] [-n hash_func_num] [-i iteration] [-s min_similarity] [-l min_mz] [-r max_mz] [-c clustering_prefix] .`
  3. Typical example: `./mscrush_on_general_charge -f ../mgf/D01*part*.mgf -t 40 -n 15 -i 100 -s 0.65 -l 200 -r 2000 -c ../clusters/clusters `. You will find clustering result files under the path specified through `-c` flag, in our case, it is 5 clustering result files with name prefix *clusters* under `../clusters`.
  4. Description
  Type `./mscrush_on_general_charge -h` to see full list of command options.
        * -f,    --files (required)
            > MGF files to cluster.
        
        * -i,    --iteration
            > Clusbering iteration.
            This parameter is optional. The default value is '100'.

        * -n,    --hash
            > Hash functions per hash table.
     This parameter is optional. The default value is '15'.

        * -t,    --thread
            > Threads to use.
     This parameter is optional. The default value is '20'.

        * -l,    --min_mz
            > Minimum cosine similalrity for clustering.
     This parameter is optional. The default value is '200'.

        * -r,    --max_mz
            > Minimum cosine similalrity for clustering.
     This parameter is optional. The default value is '2000'.

        * -s,    --similarity
            > Minimum cosine similalrity for clustering.
     This parameter is optional. The default value is '0.65'.

        * -c,    --clustering_prefix
            > Clustering result file prefix.
     This parameter is optional. The default value is 'clustering'.

        * -d,    --delimiter
            > Delimiter to separate MS2 titles in clusters.
     This parameter is optional. The default value is '|'.


- It is very important to select proper paramters for msCRUSH to cluster similar spectra with high sensitivity and specificity. Specifically:
  1. `hash_func_num` controls the collision probability of two spectra in a single hash table, the more hash functions, the smaller collision probabiilty, the more clusters, the faster msCRUSH program. 10 is a good starting point for datasets with size less than 1 million, while 15 is a good starting point for datasets with size around 10 million.
  2. `iteration` controls the number of hash tables to use with the aim to increase the probabiilty of two similar spectra to collide in at least one hash table, the more iterations, the higher collision probability, the less clusters, the slower mcCRUSH program. 100 iterations is a good starting point to play with.
  3. `min_similarity` is the minimum similarity to consider two spectra to be similar and merge into a consensus spectrum, the higher min similarity, the more clusters, the slower msCRUSH. 0.5 through 0.65 is recommended.
  4. `min_mz` is the minimum mz to consider a peak. 
  5. `max_mz` is the maximum mz to consider a peak.


- The output file locate in the `--clustering_prefix` folder, which holding all clusters with respective charge and each file should contains one line per cluster with the format like:   
  ID    Titles    
  0    042913-SludgeD1_13-18, Cmpd 25, +MSn(815.43), 0.07 min    
  1    042913-SludgeD1_13-18, Cmpd 21300, +MSn(814.31), 47.09 min    
  2    042913-SludgeD1_19-24, Cmpd 38786, +MSn(1351.9), 91.14 min    
  3    042913-SludgeD1_13-18, Cmpd 37999, +MSn(1357.62), 90.13 min   
  4    042913-SludgeD1_25-30, Cmpd 8384, +MSn(866.89), 19.79 min   
  5    042913-SludgeD1_25-30, Cmpd 10087, +MSn(732.39), 22.75 min    
          ...     
  60    042913-SludgeD1_25-30, Cmpd 23895, +MSn(810.04), 51.89 min    
  61    042913-SludgeD1_1-6, Cmpd 15070, +MSn(595.83), 35.72 min    
  62	  042913-SludgeD1_19-24, Cmpd 6029, +MSn(931.71), 15.07 min   
  63	  042913-SludgeD1_13-18, Cmpd 10135, +MSn(591.77), 23.32 min    
  64	  042913-SludgeD1_25-30, Cmpd 21901, +MSn(639.34), 46.86 min|042913-SludgeD1_25-30, Cmpd 21979, +MSn(639.33), 47.07 min   
  65    042913-SludgeD1_25-30, Cmpd 4893, +MSn(882.46), 12.13 min     
          ...     
  in which each line start with the cluster-id and a `--delimiter` separated list of the title of the spectra that belong to that cluster.    
  
  
## Generate consensus spectra
Note that writing consensus spectra (cs) to disk in MGF format can be time consuming, so if consensus spectra is needed, run scripts below.
1. `cd bin`
2. Usage: `./generate_consensus_spectrum_for_mscrush -c mscrush_cluster(s) -f mgf_files(s) [-t consensus_title] [-p consensus_path_prefix] [-d decimal_place].`
3. Typical example: `./generate_consensus_spectrum_for_mscrush -c ../clusters/clusters-c*.txt -f ../mgf/D01*part*.mgf -d 7 -t CONSENSUS -p ../consensus/consensus`. You will find consensus spectra files, each of which matches to a clustering file of a specific charge state, in our case, it is 5 files with name prefix *consensus* under dir `../consensus`
4. Description
  Type `./generate_consensus_spectrum_for_mscrush -h` to see full list of command options.
    * -c,    --clusters      (required)
        > Clustering files by msCRUSH.

    * -f,    --files (required)
        > MGF files.
        
    * -d,    --decimal
        > Decimal places for numbers.
   This parameter is optional. The default value is '3'.

    * -s,    --separator
        > Delimiter to separate MS2 titles in clusters
   This parameter is optional. The default value is '|'.

    * -t,    --consensus_title
        > Consensus spectrum title prefix.
   This parameter is optional. The default value is 'CONSENSUS'.

    * -p,    --consensus_path_prefix
        > Consensus result file prefix to write into
   This parameter is optional. The default value is 'consensus'.

## Citation
[msCRUSH: Fast Tandem Mass Spectral Clustering Using Locality Sensitive Hashing](https://pubs.acs.org/doi/10.1021/acs.jproteome.8b00448)

## Questions
Please contact Lei Wang (wang558@indiana.edu) for assistance.
## Acknowledgement
This work was supported by the NIH grant 1R01AI108888 and the Indiana University Precision Health Initiative.

