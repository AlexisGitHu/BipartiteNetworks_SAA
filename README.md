# BipartiteNetworks_SAA
This repository is provided to keep all R scripts related to Spectral Analysis Normalization on Bipartite Networks in Ecology

## Packages
### Versions
All of the following libraries with their version would be needed to execute the code in "spectral_distance_code_trials". If no version is specified, then its the latest version at the moment:
dplyr 
```
forcats 

ggplot2==3.4.2

ggrepel

igraph

kcorebip (included in files)

patchwork==1.1.2

purrr

tidyr

vegan
```
Newest versions of the packages still work with the implementations but avoid them unless specified by kcorebip package installation.

### Repository
As this work adds content to Javier García Algarra original  repositories Kcorebip and Spectnull, this repository provides a submodule for both.

However, all the changes to those repositories are included in the spectull_changed directory (Kcorebip has not been modified).

### Sources files
To execute the key calculations of the project, the only package that needs to be sourced or executed from terminal or using an IDE is ```spectral_distance_code_trials.R```. This will generate all the folders needed to save results and make calculations.

### Kcorebip installation
Even if source files are included as a submodule, the easiest way to install the package is using the following command in an R termianl or IDE:
```R
remove.packages("kcorebip")
library("devtools")
install_github("jgalgarra/kcorebip")
```

### Comments
If you clone the repository and don't install the kcorebip package, it won't work.
If you clone the repository and the kcorebip directory does not have as its name 'kcorebip' it may not work.


