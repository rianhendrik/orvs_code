#SnaSS - Systematically non-aligned spatial sampling

method <- "SnaSS"
method_full <- "Systematically non-aligned spatial sampling"
snapp <- spsample(ch_sp, nsamp, "nonaligned")@coords #systematic non-aligned point pattern

#Plotting is only to see if things run as they should
if (miplot == 1){
  plot(ch_sp, main = paste0(curr_village), cex.main = 3)
  points(snapp, col = "red", pch = 19)

  plot(ch_sp, main = paste0(curr_village), cex.main = 3)
  points(village, main = paste0(curr_village), cex.main = 3)
  points(snapp, col = "red", pch = 19)
}

snass <- as.data.frame(snapp)
ghu_snass <- ghu

snass_coverage <- function(ghu){
  
  nsamp_snapp <- nsamp
  crits <- c(0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.90, 0.89, 0.88)
  move_along1 <- F
  move_along2 <- F
  j = 0
  big_counter <- 0
  redo = T
  while (redo == T){
    big_counter <- 0
    j = j + 1
    crit <- crits[j]
    if (is.na(crit) == T){print("We need a lower crit")}
    if (is.na(crit) == T)break
    cond2 <- F
    while (cond2 == F){
      cond3 = F
      points_to_nsamp_prop = c()
      counter = 0
      grid_size <- ceiling(sqrt(nsamp))  
      while (cond3 == F){
        counter <- counter + 1
        snapp <- spsample(ch_sp, nsamp_snapp, "nonaligned")@coords
        points_to_nsamp_prop[counter] <- nrow(snapp)/nsamp #We want this to be as close to 1 as possible.
        if (points_to_nsamp_prop[counter] > 1){cond3 <- T}
        # nsamp_snapp = nsamp_snapp + 1
      }
      if (length(points_to_nsamp_prop) == 1){
        # print("Move along to special case 1.")
        move_along1 = T
      }
      if (length(points_to_nsamp_prop) == 1)break
      if (points_to_nsamp_prop[length(points_to_nsamp_prop)-1] > crit){cond2 = T}
      big_counter <- big_counter + 1
      if (big_counter > 250)break
    }
    if (big_counter < 250)break
  }
  
  diff <- nrow(snapp) - nsamp
  if (diff > 0){
    random_sample <- sample(seq_len(nrow(snapp)), diff)
    snapp <- snapp[-random_sample,]
  }
  
  snass <- as.data.frame(snapp)
  ghu_snass <- ghu
  
  if (miplot == 1){
    plot(village, main = paste0(curr_village))
    points(snapp, col = "red", pch = 19)
    points(ghu_snass[,2:3], col = "blue", pch = 19)
  }
  
  sample <- matrix(NA, nrow = nrow(snass), ncol = 2) #First col is for house number, and second row is for house loc in ghu
  colnames(sample) <- c("house num", "location in ghu")
  inter <- rbind(data.frame("idx" = 0, x = 0, y = 0), ghu_snass[, 1:3])
  for (i in seq_len(nrow(snass))){
    inter[1,] <- c(0, snass[i,])
    rownames(inter) <- inter[, 1] #This is because some rownames were greater than 712, which causes a problem once we extract the sampled houses from ghu
    dists <- as.matrix(dist(
      inter %>% 
        dplyr::select(x, y),
      method = 'euclidean'))[,1]*metres
    if (miplot == 1){
      png(paste0("GIFs/SnaSS/images/", i, ".png"))
      plot(village, main = paste0(method_full, "\n Iteration: ", i))
      points(inter[2:nrow(inter),2:3], col = "blue", pch = 19)
      points(inter[1,2:3], col = "red", pch = 19)
      points(inter[which(dists == min(dists[2:length(dists)])), 2:3], col = "green", pch = 19)
      legend("topright", bty = "n", pch = 19, cex = 1.1,
             col=c("red", "green"), c("Sampled point", "Nearest house"))
      dev.off()
    }
    
    sampled_house <- as.integer(names(which(dists == min(dists[2:length(dists)]))[1])) #extract thi first distance, in case there are two houses tied for the minimum distance
    sample[i,] <- c(sampled_house, which(ghu[, 1] == sampled_house)) #the sampled point (point closest to the ith uniform point.) - the house number as it appears in ghu_snass
    inter <- inter[-which(dists == min(dists[2:length(dists)]))[[1]],] #removing the sampled point (to avoid it being selected again.)
  }
  return(sample)
}


tictoc::tic()
#nsim is defined in master script
boot_dists <- matrix(NA, nsim, 1) #Les distances are the walking distance for each iteration
if (mst == T){boot_dists <- matrix(NA, nsim, 2)
  colnames(boot_dists) <- c("TSP distances", "MST distances")}
coverage <- matrix(NA, nsim, 2)
counter <- 0
for (iter in 1:nsim){
  counter <- counter + 1
  snass_sample <- snass_coverage(ghu)
  snass <- ghu[snass_sample[, 2],]
  distances<-matrix(NA, max(snass$grp_num), 1)
  colnames(distances)[1]<-"TSP distances"
  if (mst == T){distances<-matrix(NA, max(srs$grp_num), 2)
    colnames(distances)<-c("TSP distances", "MST distances")}
  for (i in sort(unique(snass$grp_num))){
    snass_gi <- snass %>% #gi denoted group i
      filter(snass$grp_num == i) 
    snass_gi$idx <- seq_len(nrow(snass_gi))
    nodes <- c(snass_gi[, 1], nrow(snass_gi)+1) #Making space for stopping point at end of list
    locs1 <- cbind(nodes, rbind(sp[i, 1:2], snass_gi[, 2:3])) #adding the stopping point

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
  coverage[iter,] <- c((sum(snass$ndogs)/total_dogs), (nrow(snass)/n))
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