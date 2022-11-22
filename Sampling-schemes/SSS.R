#SpatStratS - Spatial stratified sampling
method <- "SpatStratS"
method_full <- "Spatial Stratified Sampling"

spatstrat_coverage <- function(ghu){
  ghu_sstrat <- ghu
  sample <- matrix(NA, nrow = nrow(ghu), ncol = 2) #First col is for house number, and second row is for house loc in ghu
  colnames(sample) <- c("house num", "location in ghu")
  redo = T
  prop_sampled = 0.7
  while (redo == T){
    inner_counter <- 0
    for (i in 1:max(ghu$grp_num)){
      curr <- ghu_sstrat %>%
        filter(ghu_sstrat$grp_num == i) #The first house is also the stopping point.
        ncurr <- nrow(curr)
        duplicates <- any(duplicated(curr[,2:3]))
      if (ncurr > 0 & duplicates == F){ #This any(duplciated) searches to see if any house locations are duplicated, resulting in a radius of 0
        if(ncurr > 1){
          radius <- max(as.matrix(dist(
            curr %>% 
              dplyr::select(x, y),
            method = 'euclidean'))[,1])
          
          cw <- disc(radius, c(curr[1, 2], curr[1, 3]))
          
          rand_points <- as.data.frame(runifpoint(ceiling(round(prop_sampled*ncurr)), cw))
          
          inter <- rbind(data.frame("idx" = 0, x = 0, y = 0), curr[, 1:3])
          for (j in seq_len(nrow(rand_points))){
            ###
            inner_counter <- inner_counter + 1
            ###
            inter[1,] <- c(0, rand_points[j,])
            rownames(inter) <- inter[, 1] #This is because some rownames were greater than 712, which causes a problem once we extract the sampled houses from ghu
            dists <- as.matrix(dist(
              inter %>% 
                dplyr::select(x, y),
              method = 'euclidean'))[,1]*metres
            
            if (miplot == 1){
              png(paste0("GIFs/spatstrat/images/", j, ".png"))
              plot(curr[,2:3], main = paste0(method_full, "\n Iteration: ", j))
              points(inter[2:nrow(inter),2:3], col = "blue", pch = 19)
              points(inter[1,2:3], col = "red", pch = 19)
              points(inter[which(dists == min(dists[2:length(dists)])), 2:3], col = "green", pch = 19)
              dev.off()
            }
            
            #points(inter[which(dists == min(dists[2:length(dists)])), 2:3], col = "green", pch = 19)
            sampled_house <- as.integer(names(which(dists == min(dists[2:length(dists)]))[1])) #extract thi first distance, in case there are two houses tied for the minimum distance
            sample[inner_counter,] <- c(sampled_house, which(ghu[, 1] == sampled_house)) #the sampled point (point closest to the ith uniform point.) - the house number as it appears in ghu_sstrat
            inter <- inter[-c(which(dists == min(dists[2:length(dists)]))[[1]]),] #removing the sampled point (to avoid it being selected again.)
          }
        }else{
          inner_counter <- inner_counter + 1
          sample[inner_counter,] <- c(curr[[1]], which(ghu[, 1] == curr[1]))}
      }
    }
    sample <- na.omit(sample) #remove all the unused rows from sample
    diff = nsamp - nrow(sample)
    if (diff > 0){redo = T
                  prop_sampled <- prop_sampled + 0.01
                  sample <- matrix(NA, nrow = nrow(ghu) + abs(diff), ncol = 2)}
    if (diff < 0){sample = sample[-c(round(runif(abs(diff), 1, nrow(sample)))),]
                  redo = F}
    if (diff == 0){redo = F}
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
  sstrat_sample <- spatstrat_coverage(ghu)
  sstrat <- ghu[sstrat_sample[, 2],]
  find_nulls <- which(is.na(sstrat[, 1]) == T)
  if (length(which(is.na(sstrat[, 1]) == T)) > 0){
    sstrat <- sstrat[-which(is.na(sstrat[, 1]) == T),]
  }
  distances<-matrix(NA, max(sstrat$grp_num), 1)
  colnames(distances)[1]<-"TSP distances"
  if (mst == T){distances<-matrix(NA, max(srs$grp_num), 2)
    colnames(distances)<-c("TSP distances", "MST distances")}
  for (i in sort(unique(sstrat$grp_num))){
    sstrat_gi <- sstrat %>% #gi denoted group i
      filter(sstrat$grp_num == i) 
    sstrat_gi$idx <- seq_len(nrow(sstrat_gi))
    nodes <- c(sstrat_gi[, 1], nrow(sstrat_gi)+1) #Making space for stopping point at end of list
    locs1 <- cbind(nodes, rbind(sp[i, 1:2], sstrat_gi[, 2:3])) #adding the stopping point

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
  coverage[iter,] <- c((sum(sstrat$ndogs)/total_dogs), (nrow(sstrat)/n))
  progress <- (counter/nsim)*100
  if (print_progress == T){
    print(paste0("Process is ",progress,"% complete."))
  }
  #To see progress when console gets clogged with Concorde's messages.
  if (write_progress == T){
    progress_file <- file(paste0("output_script_", process_code, ".txt"))
    writeLines(paste0("Process is ",progress,"% complete."), progress_file)
    close(progress_file)
  }
}

if (save_boot_dists_in_sample_script == T){
  saveRDS(boot_dists, file = paste0("location/where/you/want/to/save/walking/distance/output"))
}