#Speed and accuracy comparison for TST algorithms
tst_methods <- c("nn", "repetitive_nn","nearest_insertion", "farthest_insertion",
                 "cheapest_insertion", "arbitrary_insertion", "concorde",
                 "2-opt")

concordepath <- "your/path"
linkernpath <- 'your/path'

nresults <- 100000 #we want 100 000 results
tst_speed  <- matrix(NA, nrow = nresults, ncol = length(tst_methods))
tst_distance <-  matrix(NA, nrow = nresults, ncol = length(tst_methods))
nhouses <- matrix(NA, nrow = nresults, ncol = 1)
colnames(tst_speed) <- tst_methods
colnames(tst_distance) <- tst_methods
colnames(nhouses) <- 'number of houses'

#We are using the SRS sampling scheme to perform the sampling here. Any scheme can be used, however.
srs_coverage <- function(ghu){
  samp_idx <- sample(seq_len(nrow(ghu)), 0.7*nrow(ghu))
  samp <- ghu[samp_idx,]
  remainder <- ghu[-samp_idx,]
  return(list("samp" = samp, "remainder" = remainder))
}

counter <- 0
while(is.na(nhouses[nresults]) == T){
  srs_func <- srs_coverage(ghu)
  srs <- srs_func$samp
  distances <- matrix(NA, max(srs$grp_num), 1)
  for (i in sort(unique(srs$grp_num))){
    srs_gi <- srs %>%
      filter(srs$grp_num == i)
    srs_gi$idx <- seq_len(nrow(srs_gi))
    nodes <- c(srs_gi[, 1], nrow(srs_gi)+1) #Making space for stopping point at end of list
    locs1 <- cbind(nodes, rbind(sp[i, 1:2], srs_gi[, 2:3])) #adding the stopping point

    dist_mat <- dist(
      locs1 %>% select(x, y),
      method = 'euclidean'
    )*metres

    if (nrow(locs1)>3){
      counter <- counter + 1 #only count here, so that the result matrix is only appended if there is a result to record (nhouses>=4)
      nhouses[counter] <- nrow(locs1)
      tsp_prob <- TSP(dist_mat)
      tsp_solver_counter <- 0
      for (tsp_solver in tst_methods){

        #To set correct path for concorde and linkern
        if (tsp_solver == 'concorde'){
          concorde_path(concordepath)
        }else{concorde_path(linkernpath)}

        tsp_solver_counter <- tsp_solver_counter + 1
        tst_start <- Sys.time()
        tour <- invisible(suppressWarnings(solve_TSP(
          tsp_prob,
          method = tsp_solver,
          start = 1
        )))
        #Saving the results
        tst_distance[counter, tsp_solver_counter] <- tour_length(tour)
        tst_speed[counter, tsp_solver_counter] <- Sys.time() - tst_start
      }
    }
  }

  progress <- ((nresults - sum(is.na(nhouses)))/nresults)*100
  progress_file <- file(paste0("output_progress/tst_speed_test_intial.txt"))
  writeLines(paste0("Process is ",round(progress, 2),"% complete."), progress_file)
  close(progress_file)
}


# Write results to na excell file
# The below data is then used to compare the performance of each algorithm. The TST comparison plots 
# in the masters document was generated using the below data.
write.csv(tst_distance, 'tst_distance.csv')
write.csv(tst_speed, 'tst_speed.csv')
write.csv(nhouses, 'nhouses.csv')