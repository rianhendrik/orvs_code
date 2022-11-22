save_boot_dists_in_sample_script <- F #We do not want to save the walking distances for each simulation indivually, but rather in one database
print_progress <- F #This is the progress in the sampling scripts themselves
write_progress <- T #This is the progress in the sampling scripts themselves
output_png_bootdists <- F #We do not want the sampling script to output nsim*length(thetas) plots!

#for data preparation run
exclude_nodoghouses <- T
curr_village <- "Machochwe"
sp_algorithm <- "kmeans" #Only kmeans for this script
walking_radius <- 200 #as used in core paper
#for sampling scheme
sampling_scheme <- "USS"
tsp_solver <- "farthest_insertion"
coverage_goal <- '40%'
coverage_goal_num <- 0.699
# coverage_goal_num <- 0.6
process_code <- paste0(sampling_scheme, '10000_coverage_', coverage_goal, '_test')

bandwidth = 50 #for a second round of KWSS tests

# Loop to find the random seed which gives us our desired coverage!
thetas <- seq(20,30) 
coverages_houses <- 0
nrandseeds <- 1000
temp_res <- matrix(NA, nrow = nrandseeds, ncol = 3)
fin_res = c()
colnames(temp_res) = c('theta', 'seed', 'coverage_houses')
coverage_target_met = F

for (theta in thetas){
  if (length(fin_res) == 0){
    for (ii in 1:nrandseeds){
      set.seed(ii) #setting seed for kmeans to give X% coverage
      source("R code/data_preparation.R")
      nsim <- 1
      nsamp <- round(0.7*nrow(ghu)) #We want to sample 70% of the accessible houses. Maybe we can play around with this value?
      source(paste0("R code/Sampling schemes/", sampling_scheme,".R"))  #Here is where the sampling scheme is determined
      coverages_houses <- coverage[,2]
      temp_res[ii,] <- c(theta, ii, coverages_houses)
      if (ii%%990 == 0){print(max(temp_res[,3], na.rm = T))}
      if (coverages_houses >= coverage_goal_num){
        print(paste0("Warning! The goal is ", coverage_goal_num, ", but you obtained a coverage of ",
                     coverages_houses))
      }
    }
    fin_res = rbind(fin_res,temp_res[which(temp_res[,3] >= (coverage_goal_num)),])
    print(theta)
   }else break 
}
cat(paste0('Scheme: ', sampling_scheme, '\n',
       'Coverage goal: ', coverage_goal_num, '\n', 
       'Theta: ', fin_res[1,1], '\n',
       'Seed: ', fin_res[1,2], '\n',
       'Coverage achieved: ', round(fin_res[1,3],3), '\n'))
fin_res


#Find the probability of being close to 30
#Select one of those random seeds as the best seed to get n% coverage with theta stopping points.
#Finding the optimal theta for each sampling scheme
