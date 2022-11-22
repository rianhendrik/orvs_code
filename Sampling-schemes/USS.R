#Uniform spatial sampling (USS)
method <- "USS"
method_full <- "Uniform Spatial Sampling"

uss_coverage <- function(ghu){
  
  unifpp <- runifpoint(nsamp, ch)
  
  if (miplot == 1){
    plot(village, main = paste0(curr_village))
    points(unifpp, col = "red", pch = 19)
    points(ghu[,2:3], col = "blue", pch = 19)
  }
  
  uss <- as.data.frame(unifpp)
  ghu_uss <- ghu
  remainder <- ghu
  sample <- matrix(NA, nrow = nrow(uss), ncol = 2) #First col is for house number, and second row is for house loc in ghu
  colnames(sample) <- c("house num", "location in ghu")
  inter <- rbind(data.frame("idx" = 0, x = 0, y = 0), ghu_uss[, 1:3])
  for (i in seq_len(nrow(uss))){
    inter[1,] <- c(0, uss[i,])
    rownames(inter) <- inter[, 1] #This is because some rownames were greater than 712, which causes a problem once we extract the sampled houses from ghu
    dists <- as.matrix(dist(
      inter %>% 
        dplyr::select(x, y),
      method = 'euclidean'))[,1]*metres
    if (miplot == 1){
      png(paste0("GIFs/uss/images/", i, ".png"))
      plot(village, main = paste0(method_full, "\n Iteration: ", i))
      points(inter[2:nrow(inter),2:3], col = "blue", pch = 19)
      points(inter[1,2:3], col = "red", pch = 19)
      points(inter[which(dists == min(dists[2:length(dists)])), 2:3], col = "green", pch = 19)
      legend("topright", bty = "n", pch = 19, cex = 1.1,
             col=c("red", "green"), c("Sampled point", "Nearest house"))
      dev.off()
    }
    sampled_house <- as.integer(names(which(dists == min(dists[2:length(dists)]))[1])) #extract thi first distance, in case there are two houses tied for the minimum distance
    sample[i,] <- c(sampled_house, which(ghu[, 1] == sampled_house)) #the sampled point (point closest to the ith uniform point.) - the house number as it appears in ghu_uss
    inter <- inter[-which(dists == min(dists[2:length(dists)]))[[1]],] #removing the sampled point (to avoid it being selected again.)
    remainder <- remainder[-which(remainder[, 1] == sampled_house),]
    deliver <- list("sample" = sample, "remainder" = remainder)
  }
  return(deliver)
}

#nsim is defined in master script
boot_dists <- matrix(NA, nsim, 1) #Les distances are the walking distance for each iteration
if (mst == T){boot_dists <- matrix(NA, nsim, 2)
  colnames(boot_dists) <- c("TSP distances", "MST distances")}
coverage <- matrix(NA, nsim, 2)
counter <- 0
for (iter in 1:nsim){
  counter <- counter + 1
  func_out <- uss_coverage(ghu)
  uss_sample <- func_out$sample
  uss_remainder <- func_out$remainder #Use uss_remainder for the SVP scheme
  uss <- ghu[uss_sample[, 2],]
  distances<-matrix(NA, max(uss$grp_num), 1)
  colnames(distances)[1]<-"TSP distances"
  if (mst == T){distances<-matrix(NA, max(uss$grp_num), 2)
    colnames(distances)<-c("TSP distances", "MST distances")}
  for (i in sort(unique(uss$grp_num))){
    uss_gi <- uss %>% #gi denoted group i
      filter(uss$grp_num == i) 
    uss_gi$idx <- seq_len(nrow(uss_gi))
    nodes <- c(uss_gi[, 1], nrow(uss_gi)+1) #Making space for stopping point at end of list
    locs1 <- cbind(nodes, rbind(sp[i, 1:2], uss_gi[, 2:3])) #adding the stopping point

    dist_mat<-dist(
      locs1 %>% dplyr::select(x, y),
      method = 'euclidean'
    )*metres

    tsp_prob <- TSP(dist_mat)

    tour <- solve_TSP(
      tsp_prob,
      method = "cheapest",
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
  coverage[iter,] <- c((sum(uss$ndogs)/total_dogs), (nrow(uss)/n))
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