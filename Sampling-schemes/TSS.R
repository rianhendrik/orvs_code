#TradStratS - Traditional stratified sampling
method <- "TradStratS"
method_full <- "Traditional stratified sampling"

nsamp <- round(0.7*nrow(ghu))

tradstrat_coverage <- function(ghu){
  prop_goal <- 0.72
  cond = F
  prop_counter = 0
  while (cond == F){
    prop <- 1
    ghu_tss <- ghu #tss = tradtional stratified sampling
    sps <- 1:max(ghu_tss$grp_num)
    counter = 0
    while (prop > prop_goal){
      remove <- sample(sps, 1)
      ghu_tss <- ghu_tss %>%
        filter(ghu_tss$grp_num != remove) 
      prop <- nrow(ghu_tss)/nrow(ghu)
      counter <- counter + 1
      # print(prop)
    }
    if (round(prop,3) > 0.7){cond = T}else{prop_counter = prop_counter + 1}
    if (prop_counter > 250){prop_goal = prop_goal + 0.05}
  }
  # plot(ghu_tss$x, ghu_tss$y)
  # print(prop)
  return(ghu_tss)
}


ghu_tss <- tradstrat_coverage(ghu)

diff <- nrow(ghu_tss) - nsamp
if (diff > 0){
  random_sample <- sample(seq_len(nrow(ghu_tss)), diff)
  ghu_tss <- ghu_tss[-random_sample,]
}

#nsim is defined in master script
boot_dists <- matrix(NA, nsim, 1) #Les distances are the walking distance for each iteration
if (mst == T){boot_dists <- matrix(NA, nsim, 2)
  colnames(boot_dists) <- c("TSP distances", "MST distances")}
coverage <- matrix(NA, nsim, 2)
counter <- 0
for (iter in 1:nsim){
  counter <- counter + 1
  tstrat <- tradstrat_coverage(ghu)
  
  diff <- nrow(tstrat) - nsamp
  if (diff > 0){
    random_sample <- sample(seq_len(nrow(ghu_tss)), diff)
    tstrat <- tstrat[-random_sample,]
  }
  
  distances<-matrix(NA, max(tstrat$grp_num), 1)
  colnames(distances)[1]<-"TSP distances"
  if (mst == T){distances<-matrix(NA, max(srs$grp_num), 2)
    colnames(distances)<-c("TSP distances", "MST distances")}
  for (i in sort(unique(tstrat$grp_num))){
    tstrat_gi <- tstrat %>%
      filter(tstrat$grp_num == i) 
    tstrat_gi$idx <- seq_len(nrow(tstrat_gi))
    nodes <- c(tstrat_gi[, 1], nrow(tstrat_gi)+1) #Making space for stopping point at end of list
    locs1 <- cbind(nodes, rbind(sp[i, 1:2], tstrat_gi[, 2:3])) #adding the stopping point

    dist_mat<-dist(
      locs1 %>% dplyr::select(x, y),
      method = 'euclidean'
    )*metres

    tsp_prob <- TSP(dist_mat)

    tour <- solve_TSP(
      tsp_prob,
      method = tsp_solver,
      start = 1
    )
    distances[i,1]  <- tour_length(tour)

    if (mst == T){
      nodes_mst <- locs1[,1]
      arcs <- gtools::combinations(max(nodes_mst), 2, nodes_mst)
      d_mst <- round(matrix(dist_mat, length(dist_mat),1),1)
      arcs <- cbind(arcs, d_mst)
      mst_i <- getMinimumSpanningTree(nodes_mst, arcs, algorithm  = "Prim",
                                      show.data  = F, show.graph  = F)
      distances[i,2] <- mst_i$weight
    }
  }
  boot_dists[iter,] <- colSums(distances, na.rm<-T)/1000 #some entries will be NA, if a group did not appear in the sample.
  coverage[iter,] <- c((sum(tstrat$ndogs)/total_dogs), (nrow(tstrat)/n))
  progress <- (counter/nsim)*100
  if (print_progress == T){
    print(paste0("Process is ",progress,"% complete."))
  }
  #To see progress when console gets clogged with Concorde's messages.
  if (write_progress == T){
    progress_file <- file(paste0("output_script", process_code, ".txt"))
    writeLines(paste0("Process is ",progress,"% complete."), progress_file)
    close(progress_file)
  }
}


if (save_boot_dists_in_sample_script == T){
  saveRDS(boot_dists, file = paste0("location/where/you/want/to/save/walking/distance/output"))
}