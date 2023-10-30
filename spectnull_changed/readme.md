# Instructions to reproduce results

# Brief description of idea

The key idea behind the results is to prove with a measurable scale-free spectral distance that bipartite networks, both weighted and non-weighted, can be comparable with one another.

Therefore, to address the necessity of this comparison when working with bipartite networks generated through null models, there have been developed some functions that help us achieve this.

It is vital to mention that weighted networks and non-weighted need to be consider as different cases. Otherwise, comparisons will not yield accurate results when comparing different null model generated bipartite networks. Also, spectral distances between Laplacian matrices and Adjacency matrices should be consider different cases for the same reason.

All functions to address these purposes can be found in the "spectral_distance_code_trials.R" which uses cran packages and two non-cran packages (kcorebip and spectull) following the results of J. G. Algarra [cite].

The generation of null models is addressed in J. G. Algarra "spectnull" repository. Therefore, this will be needed to execute the script mentioned above.

From this point on, the methods and their implementation will be addressed.

# Methods

1. Boxplots generation with different normalization criteria. These criteria include:
    1. Min-Max normalization
        1. **By null model type – single network**. In these calculations, we only use the Min-Max criterion to normalize spectral distances across each single network. This holds that the maximum and minimum spectral distance is provided by the network analyzed.
            1. The method to do the box-plots graphs with the distances normalized by this criterion is called ```min_max_normalization()```.
        2. **By maximum and minimum value across all null model types – single networks.** In these calculations, we also use the Min-Max criterion to normalize but we normalize the spectral distances of each network by the maximum and minimum value found across all networks and null model types.
            1. The method to do the box-plots graphs with the distances normalized by this criterion is called ```min_max_normalization_model_independant()```.
        3. **By maximum and minimum value across all null model types – all networks**. In these calculations, we also use the Min-Max criterion to normalize but we normalize considering all spectral distances across all null model types among all networks.
            1. The method to do the box-plots graphs with the distances normalized by this criterion is called ```min_max_normalization_all_distances()```.
    2. Bounds Normalization: in this case, we know, by definition, that the spectral distance between two matrices is the squared root of the differences squared of the sorted ascendant eigenvalues of those two matrices. As we are only considering bipartite networks, specifically Adjacency and Laplacian matrices, we can calculate analytically the spectral distance bounds for each type and normalize all distances by the upper and lower bounds.
        1. For each network analyzed, we must calculate these bounds in the same manner, this is done by calling the ```calculate_bounds_adjacency()``` for the bounds for Adjacency matrices and 'save\_parameters\_laplacian' for the Laplacian matrices.
        2. These bounds will vary depending on the networks given, therefore we must update those values till we reach the maximum and minimum among all the networks.
        3. Once those are calculated, we can then calculate the boxplots with the methods ```bounds_normalization()``` for Adjacency matrix and ```bounds_normalization_laplacian()``` for the Laplacian matrix.
    3. IQR Normalization: it is important to mention that in this case the spectral distance abandons it geometric meaning and take a variability explanation. This normalization transforms the spectral distance into a variability rate of null models. Methods II and IV do not have much explainable meaning, but methods I and II do have it but in different ways.
        1. **By null model type – single network**. This is the robust normalization only considering the IQR by null model type spectral distances for each null model type. Therefore, the distances are independently normalized by network and by null model type.
            1. This is achieved by calling the ```iqr_normalization()```.
        2. **IQR for each null model type across all networks**. In this case, we must normalize with robust normalization by calculating the IQR of the spectral distances of each null model across all networks and then normalize the values. This is explained in the next section:
            1. For each network, there are going to exist distances associated with each null model type. If we consider the whole group of distances with respect to one single model type across all networks and calculate its IQR, then we can normalize each distance belonging to a network and a null model by that IQR.
            2. This is achieved by calling ```iqr_normalization_all_distances_bymodeltype()```.
        3. **IQR for all null model types in each single network.** In this case, it is not considered the separation by null model type, and we calculate the IQR for all spectral distances generated between all null models for a single network. We then normalize those values by robust normalization of all spectral distances of that single network.
            1. This is done by calling the ```iqr_normalization_all_distances_bynetwork()``` method.
        4. **All null model types across all networks.** In this case, it is not considered the separation by null model type as we take the distances among all networks. Therefore, we calculate the robust normalization of spectral distances by considering all the spectral distances of all networks.
            1. This is achieved by calling the ```iqr_normalization_all_distances()```.
    4. Covariance Normalization: with this normalization happens the same as with IQR normalization, as eigenvalues from the covariance matrix represents variability among null model types.
        1. In this context, it will have no sense to calculate this among all networks, as we want to calculate the variability with the covariance matrix of the network considering the random variables as the spectral distances vectors outputted by each null model.
        2. The method to calculate this is called ```covariance_eigenvalues_normalization()```.
2. Scatterplots of mean values of Adjacency and Laplacian spectral distances. This is done for non-weighted networks and for the weighted ones.
    1. MinMax II type normalization.
        1. As a result of the boxplots obtained by all the normalization criteria, the only one who preserves the notion of geometric distance, is scale-free and keeps the non-normalized relative separation of the null models with the real network is the Min-Max II.
        2. As a result, to compare which null model type generates the most closed to reality network, compared with the real one, we generate a scatterplot. This scatterplot confronts the mean values of the spectral distance of each network and each null model type of the Adjacency matrix and the Laplacian matrix.
    2. Covariance type normalization.
        1. One could argue that those mean values are not representative as we do not know the variability of the data neither we can assure it has a normal distribution as the histograms are asymmetric. To account for that, we use the covariance normalization.
        2. In this way, we can compute the scatterplots of the rates of variability by taking the mean of the rate of variability of the spectral distances of the Adjacency Matrix and the Laplacian Matrix.
        3. Following this way, we can see in an experimental way a more realistic approach on how the variability of the means is clustered without any assumption.
    3. Both scatterplots are calculated using the same method ```plot_clustering_spd_means()``` by providing the name of the normalization wanted. Inside the method, the separation between the two types of networks is considered.
3. Histograms of mean values of Adjacency and Laplacian spectral distances.
    1. MinMax II type normalization.
        1. To provide the hypothesis with more experimental data that supports it, we also provide the histograms of the means of the adjacency spectral distance and Laplacian spectral distance normalized by this criterion.
        2. These histograms agree with the clustered data as we can see trends on them for each null model type.
    2. Covariance type normalization.
        1. In the same way as before, covariance histograms also support the trend of the rate of variability of the means.
    3. All histograms are calculated using the same method ```plot_means_histogram()``` by providing the name of the normalization wanted. Inside the method, the separation between the two types of networks is considered.
4. Null models networks generation: all this work is based on the functions provided by the package 'spectnull' from J.G. Algarra repository.

# Global variables

There are mainly four key global variables that will hold all the spectral distances. These four variables are explained below:

1. Data lists type: both lists hold the spectral distances saved by null model type. These lists are refreshed per network, so they will only contain the spectral distances associated with each null model type of one network at a time.
    1. ```data_lists```: results of the spectral distances associated with the Adjacency matrices.
    2. ```data_lists_laplacian```: results of the spectral distances associated with the Laplacian matrices.
2. Data lists of all graphs type: both lists hold the spectral distances save by null model type and by network name. These lists are not refreshed as they are intended to keep all the spectral distances calculated and grouped by network and null model.
    1. ```data_lists_all_graphs```: results of the spectral distances associated with the Adjacency matrices.
    2. ```data_lists_all_graphs_laplacian```: results of the spectral distances associated with the Laplacian matrices.

# Data

All data has been collected from "The web of life: ecological database" which can be found at this URL: [web of life](https://arxiv.org/abs/1403.2575). Specifically, all pollination networks have been retrieved and the ones that are considered radically different from the rest, have been removed.

The script takes all this networks as '.csv' files that can be found in the following path: "/spectnull\_changed/data".

An example of the names provided are:

- M\_PL\_001.csv
- M\_PL\_060\_01.csv
- M\_SD\_001.csv
- NOVELLA\_2019\_MAT.csv
- PINEDA\_2020\_MAT.csv
- RA\_HP\_001.csv

It is important to mention that these csv are not yet filtered by weighted or non-weighted network types. This is done in the script by checking if any of the edges between nodes in the network have a value greater than 0.

The results of the plots generated by the methods described above can be found in these directories correctly classified:

- base_path <- "/spectnull_changed/box_plots_data_per_network"
- cluster_path <- "/spectnull_changed/clustering_data_graph"
- histogram_path <- "/spectnull_changed/histogram_data_graph"

# Null models

The null models used to generate the bipartite networks and their matrices to calculate their corresponding spectral distances are divided into two groups, depending on the network weightness type.

1. Weighted networks: the following null models are used to calculate from the original network the network to be compared.
    1. RND: Erdos-Renyi null model for bipartite weighted networks.
    2. SYTR: stochastic generative model of world trade network based on the work published in [https://pubmed.ncbi.nlm.nih.gov/31811254/](https://pubmed.ncbi.nlm.nih.gov/31811254/).
    3. Shuffle.
    4. VAZ.
    5. SWAP
2. Non-weighted networks: the following null models are used to calculate from the original network the network to be compared.
    1. RND: Erdos-Renyi null model for bipartite networks non-weighted networks.
    2. MGEN.
    3. Shuffle.
    4. VAZ.
    5. DU\_05.
