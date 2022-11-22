# README
This directory consists of a few R scripts, a Sampling-schemes directory and a Stopping-point-code directory. 

## *master.R*
The *master.R* script should be executed first during each session, before any of the other code is executed in R. The master script imports the packages you will need for each of the scripts. You may also execute instances a sampling scheme using the *master.R* script, as it is able to source the *data_preparation.R* and the sampling scheme R scripts in the Sampling-schemes directory.
Not all of the packages in the *master.R* script are relevant to the scripts in this repository.

## Sampling-schemes (directory)
This directory contains eight R scripts. One for each of the sampling schemes used in the research. These schemes are sourced from the *master.R*, the *walking_distance_simulation_{sp_scheme}.R* scripts, and from the *finding_optimal_k_for_kmeans.R* and *k_to_stopping_point_mapping_for_kNN.R* script.

## Stopping-point-code (directory)
This directory contains two R scripts. One for each of the stopping point algorithms used in the research. These algorithms are sourced from the *data_preparation.R* script.

## *tst_speedcomp_srs.R*
This script uses the SRS sampling scheme to compare the performance of the eight TST algorithms considered in this research. Be sure to isntallt he Concorde and Linkern algorithms onto your machine (using the links in the research paper), and then source the paths of these instalations in you R script.
The data returned by this script can be used to calculate the TST comparison graphs, as well as the performance tables in the research paper.

## *data_preparation.R*
The *data_preparation.R* script prepares the village data, and processes it into a form to be used by the stopping point algorithms, and by the sampling schemes.

## *k_to_stopping_point_mapping_for_kNN.R*