# Calculating walking distance simulations for kNN based cases.
save_boot_dists_in_sample_script <- F #We do not want to save the walking distances for each simulation indivually, but rather in one database
print_progress <- F #This is the progress in the sampling scripts themselves
write_progress <- T #This is the progress in the sampling scripts themselves
output_png_bootdists <- F #We do not want the sampling script to output nsim*length(thetas) plots!
curr_village <- "Morotonga"
sp_algorithm <- "knn" #knn only
sampling_scheme <- "SpatStratS"
tsp_solver <- "farthest_insertion"
coverage_goals <- c('30%', '40$', '50%', '60%', '70%')
thetas <- c(10, 7, 5, 3, 1)
bandwidth <- 50
print_progress <- F #This is the progress in the sampling scripts themselves
parms = data.frame(cbind(coverage_goals, thetas))

for (i in 1:nrow(parms)){
  
  row <- parms[i,]
  
  coverage_goal <- row$coverage_goals
  theta <- as.integer(row$thetas)
  
  process_code <- paste0(sampling_scheme, '10000_coverage_', coverage_goal)
  
  source("R code/data_preparation.R")
  nsim <- 10000
  nsamp <- round(0.7*nrow(ghu)) #We want to sample 70% of the accessible houses. Maybe we can play around with this value?
  source(paste0("R code/Sampling schemes/", sampling_scheme,".R"))  #Here is where the sampling scheme is determined
  coverages_dogs <- coverage[,1]
  coverages_houses <- coverage[,2]
  all_boot_dists <- boot_dists
  coverages_houses
  
  write.xlsx(data.frame(distances_walked = as.numeric(unlist(boot_dists))), file = "your/path/here")
  write.xlsx(data.frame(coverage_houses = as.numeric(unlist(coverages_houses))), file = "your/path/here")
  write.xlsx(data.frame(coverage_dogs = as.numeric(unlist(coverages_dogs))), file = "your/path/here")
}

