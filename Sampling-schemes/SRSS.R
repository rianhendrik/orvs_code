#SRSS - systematic regular spatial sampling
method <- "SRSS-constant"
method_full <- "Systematic Regular Spatial Sampling Constant"

srss_coverage <- function(ghu){
  
  cond1 <- F
  counter <- 0
  grid_size <- sqrt(nsamp)
  points_to_nsamp_prop <- c()
  while (cond1 == F){
    counter <- counter + 1
    regpp <- rsyst(ch, grid_size)
    points_to_nsamp_prop[counter] <- regpp$n/nsamp #We want this to be as close to 1 as possible.
    if (points_to_nsamp_prop[counter] > 1){cond1 <- T}else{grid_size <- grid_size + 1}
  }
  
  diff <- regpp$n - nsamp
  if (diff > 0){
    random_sample <- sample(seq_len(regpp$n), diff)
    regpp <- regpp[-random_sample]
  }
  
  srss <- as.data.frame(regpp)
  ghu_srss <- ghu
  
  if (miplot == 1){
    plot(village, main = paste0(curr_village))
    points(regpp, col = "red", pch = 19)
    points(ghu_srss[,2:3], col = "blue", pch = 19)
  }
  
  sample <- matrix(NA, nrow = nrow(srss), ncol = 2) #First col is for house number, and second row is for house loc in ghu
  colnames(sample) <- c("house num", "location in ghu")
  inter <- rbind(data.frame("idx" = 0, x = 0, y = 0), ghu_srss[, 1:3])
  for (i in seq_len(nrow(srss))){
    inter[1,] <- c(0, srss[i,])
    rownames(inter) <- inter[, 1] #This is because some rownames were greater than 712, which causes a problem once we extract the sampled houses from ghu
    dists <- as.matrix(dist(
      inter %>% 
        dplyr::select(x, y),
      method = 'euclidean'))[,1]*metres
    if (miplot == 1){
      plot(village, main = paste0(curr_village, '_', i))
      # points(ghu_srss[,2:3], col = "blue", pch = 19)
      points(inter[2:nrow(inter),2:3], col = "blue", pch = 19)
      points(inter[1,2:3], col = "red", pch = 19)
      points(inter[which(dists == min(dists[2:length(dists)])), 2:3], col = "green", pch = 19)
    }
    
    sampled_house <- as.integer(names(which(dists == min(dists[2:length(dists)]))[1])) #extract thi first distance, in case there are two houses tied for the minimum distance
    sample[i,] <- c(sampled_house, which(ghu[, 1] == sampled_house)) #the sampled point (point closest to the ith uniform point.) - the house number as it appears in ghu_srss
    inter <- inter[-which(inter[,1] == sampled_house),]
  }
  return(sample)
}

#nsim is defined in master script
boot_dists <- matrix(NA, nsim, 1) #Les distances are the walking distance for each iteration
if (mst == T){boot_dists <- matrix(NA, nsim, 2)
colnames(boot_dists) <- c("TSP distances", "MST distances")}
coverage <- matrix(NA, nsim, 2)
counter <- 0
for (iter in 1:nsim){
  counter <- counter + 1
  srss_sample <- srss_coverage(ghu) # Why is a Null being added in last row?
  srss <- ghu[srss_sample[, 2],]
  distances<-matrix(NA, max(srss$grp_num), 1)
  colnames(distances)[1]<-"TSP distances"
  if (mst == T){distances<-matrix(NA, max(srss$grp_num), 2)
  colnames(distances)<-c("TSP distances", "MST distances")}
  for (i in sort(unique(srss$grp_num))){
    srss_gi <- srss %>% #gi denoted group i
      filter(srss$grp_num == i) 
    srss_gi$idx <- seq_len(nrow(srss_gi))
    nodes <- c(srss_gi[, 1], nrow(srss_gi)+1) #Making space for stopping point at end of list
    locs1 <- cbind(nodes, rbind(sp[i, 1:2], srss_gi[, 2:3])) #adding the stopping point
    
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
  coverage[iter,] <- c((sum(srss$ndogs)/total_dogs), (nrow(srss)/n))
  colnames(coverage) = c('%dogs', '%houses')
  progress <- (counter/nsim)*100
  
  if (print_progress == T){
    print(paste0("Process is ",progress,"% complete."))
  }
  #To see progress when console gets clogged with Concorde's messages.
  if (write_progress == T){
    progress_file <- file(paste0("location/where/you/want/to/save/walking/distance/output"))
    writeLines(paste0("Process is ",progress,"% complete."), progress_file)
    close(progress_file)
  }
}