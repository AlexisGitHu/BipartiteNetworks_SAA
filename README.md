# BipartiteNetworks_SAA
This repository is provided to keep all R scripts related to Spectral Analysis Normalization on Bipartite Networks in Ecology

## Reproducibility
### Versions
All of the following libraries with their version would be needed to execute the code in "spectral_distance_code_trials". If no version is specified, then its the latest version at the moment:
dplyr 

forcats 

ggplot2==3.4.2

ggrepel

igraph

kcorebip (included in files)

patchwork==1.1.2

purrr

tidyr

vegan

### Repository
As this work adds content to Javier García Algarra original  repositories Kcorebip and Spectnull, this repository provides a submodule for both.

However, all the changes to those repositories are included in the spectull_changed directory (Kcorebip has not been modified). It is worth to mention that the integration with the Kcorebip module is transparent.

### Sources files
All sources files specifications within R scripts provide the necessary functions collected from other libraries. Some functions have been modified, so use the ones provided and avoid sourcing Spectnull original repository (it may work, but this submodule is not intended to change in the future). 

### Comments
If you clone the repository and don't install the kcorebip package, it won't work.
If you clone the repository and the kcorebip directory does not have as its name 'kcorebip' it may not work.


