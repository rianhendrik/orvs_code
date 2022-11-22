#Master script
packages <- c("spatstat", "geometry", "rgl", "dplyr", "plot3D", "ggplot2",
             "optrees", "sf", "TSP", "rgdal", "osmdata", "leaflet", "maptools",
              "tictoc", "Gmedian", "stringr", "geometry", "openxlsx", 'matrixStats', 
             'maptools', 'fpc')
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

exclude_nodoghouses <- T

tsp_solver <- "cheapest"
miplot <- 0 #Global parameter for plotting. 0 => no plotting, 1 => plotting.


#Set location of your Concorde installation for Concorde TST algorithm.

#Set working directory to the location of this master script.

curr_village <- "Machochwe" #Machochwe, Morotonga, etc.
sp_algorithm <- "kmeans" #knn, kmeans.
theta <- 50 #theta is the parameter needed by each respective stopping point algorithm
walking_radius <- 200 #as used in core paper
edge_correction = T
bandwidth <- 50 #bandwidth for KWSS scheme

#Preparing the data
source("data_preparation.R") #Make sure this script is in the root directory, with master.R.