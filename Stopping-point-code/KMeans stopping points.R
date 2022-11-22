# This script is sourced from the data_preparation.R script.
# Master script should have been run first, to import packages.
d <- data.frame(cbind(x, y)) # the vectors x and y come from the data preparation script.
km <- kmeans(na.omit(d), k)
sp <- km$centers # The sp variable is then used further in the data_preparation.R script.


