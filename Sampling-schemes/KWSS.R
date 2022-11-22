method <- "KWSS"
method_full <- "Kernel weighted spatial sampling"

kwss_coverage <- function(ghu){
  ghu_kwss <- ghu
  x <- ghu_kwss$x #get x and y coordinates from ghu (accessible houses)
  y <- ghu_kwss$y

  W <- owin(xrange = c(min(x), max(x)), yrange = c(min(y), max(y))) #Define a new window around just the ghu points
  village_pp <- ppp(x, y, window = W) #create a new point pattern in this window
  ch <- convexhull(village_pp) #Make a convex hull window for this new point pattern
  
  # Adding edge corrections to the window.
  bandwidth_adj <- (bandwidth*4)/metres
  xbuff <- as.numeric(ch[["bdry"]][[1]][["x"]])
  ybuff <- as.numeric(ch[["bdry"]][[1]][["y"]])
  buff <- data.frame(cbind(xbuff, ybuff))
  buff <- rbind(data.frame(cbind(xbuff, ybuff)), buff[1,]) # Adding the first vertex again at the end, to close the polygon.
  p = st_polygon(list(as.matrix(buff)))
  pbuf = st_buffer(p, bandwidth_adj)
  ch_buff <- as.owin(pbuf)
  ch_sp <- as(ch_buff, "SpatialPolygons") #converting ch to Spatial Polygon for sp package for some sampling schemes.
  ghu_pp <- ppp(x, y, ch_buff)
  
  if (edge_correction == F){
    ch_sp <- as(ch, "SpatialPolygons") #converting ch to Spatial Polygon for sp package for some sampling schemes.
    ghu_pp <- ppp(x, y, ch)
  }
  
  kde <- density.ppp(ghu_pp, sigma = bandwidth/metres)

  if (miplot == 1){plot(kde)}
  samp_pixels <- which(is.na(kde$v) == F) #The pixels that are not null (those in the convex hulle from which we may sample)
  probs <- kde$v[samp_pixels] - min(kde$v[samp_pixels]) #Technically adding min(kde) here, since the minimum is a negative value.
  probs <- probs/sum(probs) #All probabilities add up to 1

  sample_space <- cbind(samp_pixels, probs) #defining the sampling space, and the probability of sampling each pixel
  kde_sample <- sample(sample_space[,1], round(0.7*nrow(ghu)), prob = sample_space[,2])
  loc <- data.frame(column = cbind((kde_sample %% 128) + 1, row = ceiling(kde_sample/128))) #adding 1 to column and ceiling row, in order to avoid having an index of 0.
  kwss_points <- cbind(kde$xcol[loc[,2]], kde$yrow[loc[,1]])

  if (miplot == 1){
    plot(ch_sp)
    points(ghu_pp)
    points(kwss_points, pch = 19, col = "red")
  }

  sample <- matrix(NA, nrow = nrow(kwss_points), ncol = 2) #First col is for house number, and second row is for house loc in ghu
  colnames(sample) <- c("house num", "location in ghu")
  inter <- rbind(data.frame("idx" = 0, x = 0, y = 0), ghu_kwss[, 1:3])
  for (i in seq_len(nrow(kwss_points))){
    inter[1,] <- c(0, kwss_points[i,])
    rownames(inter) <- inter[, 1] #This is because some rownames were greater than 712, which causes a problem once we extract the sampled houses from ghu
    dists <- as.matrix(dist(
      inter %>%
        select(x, y),
      method = 'euclidean'))[,1]*metres
    if (miplot == 1){
      png(paste0("GIFs/kwss/images/", i, ".png"))
      plot(village, main = paste0(method_full, "\n Iteration: ", i))
      points(inter[2:nrow(inter),2:3], col = "blue", pch = 19)
      points(inter[1,2:3], col = "red", pch = 19)
      points(inter[which(dists == min(dists[2:length(dists)])), 2:3], col = "green", pch = 19)
      legend("topright", bty = "n", pch = 19, cex = 1.1,
             col=c("red", "green"), c("Sampled point", "Nearest house"))
      dev.off()
    }

    sampled_house <- as.integer(names(which(dists == min(dists[2:length(dists)]))[1])) #extract the first distance, in case there are two houses tied for the minimum distance
    sample[i,] <- c(sampled_house, which(ghu[, 1] == sampled_house)) #the sampled point (point closest to the ith uniform point.) - the house number as it appears in ghu_shss
    inter <- inter[-which(dists == min(dists[2:length(dists)]))[[1]],] #removing the sampled point (to avoid it being selected again.)
  }
return(sample)
}

boot_dists <- matrix(NA, nsim, 1) #Les distances are the walking distance for each iteration
if (mst == T){boot_dists <- matrix(NA, nsim, 2)
  colnames(boot_dists) <- c("TSP distances", "MST distances")}
coverage <- matrix(NA, nsim, 2)
counter <- 0

for (iter in 1:nsim){
  counter <- counter + 1
  kwss_sample <- kwss_coverage(ghu)
  kwss <- ghu[kwss_sample[,2],]
  distances <- matrix(NA, max(kwss$grp_num), 1)
  colnames(distances)[1] <- "TSP distances"
  if (mst == T){distances <- matrix(NA, max(kwss$grp_num), 2)
    colnames(distances) <- c("TSP distances", "MST distances")}
  for (i in sort(unique(kwss$grp_num))){
    kwss_gi <- kwss %>%
      filter(kwss$grp_num == i)
    kwss_gi$idx <- seq_len(nrow(kwss_gi))
    nodes <- c(kwss_gi[, 1], nrow(kwss_gi)+1) #Making space for stopping point at end of list
    locs1 <- cbind(nodes, rbind(sp[i, 1:2], kwss_gi[, 2:3])) #adding the stopping point

    dist_mat <- dist(
      locs1 %>% dplyr::select(x, y),
      method = 'euclidean'
    )*metres

    if (nrow(locs1)>2){

      tsp_prob <- TSP(dist_mat)

      tour <- invisible(suppressWarnings(solve_TSP(
        tsp_prob,
        method = tsp_solver,
        start = 1
      )))
      distances[i,1] <- tour_length(tour)

      if (mst == T){
        nodes_mst <- locs1[,1]
        arcs <- gtools::combinations(max(nodes_mst), 2, nodes_mst)
        d_mst <- round(matrix(dist_mat, length(dist_mat),1),1)
        arcs <- cbind(arcs, d_mst)
        mst_i <- getMinimumSpanningTree(nodes_mst, arcs, algorithm  = "Prim",
                                        show.data  = F, show.graph  = F)
        distances[i,2] <- mst_i$weight
      }
    }else{
      distances[i,1] <- dist_mat*2    #if only one house sampled, walking distance is to the stopping point and back.
    }
  }
  boot_dists[iter,] <- colSums(distances, na.rm = T)/1000 #some entries will be NA, if a group did not appear in the sample.
  coverage[iter,] <- c((sum(kwss$ndogs)/total_dogs), (nrow(kwss)/n)) #% Dogs in first column, % houses in second
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
  saveRDS(boot_dists, file = paste0("location/where/you/want/to/save/walking/distance/output"))}