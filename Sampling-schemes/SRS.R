#SRS
method <- "SRS"
method_full <- "Simple Random Smapling"

#NOTE - this script also contains code to plot a TSP route, if desired. The other sampling scheme scripts
#does not contain this.

srs_coverage <- function(ghu){
  samp_idx <- sample(seq_len(nrow(ghu)), round(0.7*nrow(ghu)))
  samp <- ghu[samp_idx,]
  remainder <- ghu[-samp_idx,]
  return(list("samp" = samp, "remainder" = remainder))
}

#nsim is defined in master script
boot_dists <- matrix(NA, nsim, 1) #Les distances are the walking distance for each iteration
if (mst == T){boot_dists <- matrix(NA, nsim, 2)
              colnames(boot_dists) <- c("TSP distances", "MST distances")}
coverage <- matrix(NA, nsim, 2)
counter <- 0

for (iter in 1:nsim){
  counter <- counter + 1
  srs_func <- srs_coverage(ghu)
  srs <- srs_func$samp
  distances <- matrix(NA, max(srs$grp_num), 1)
  colnames(distances)[1] <- "TSP distances"
  if (mst == T){distances <- matrix(NA, max(srs$grp_num), 2)
                colnames(distances) <- c("TSP distances", "MST distances")}
  for (i in sort(unique(srs$grp_num))){
    srs_gi <- srs %>%
      filter(srs$grp_num == i)
    srs_gi$idx <- seq_len(nrow(srs_gi))
    nodes <- c(srs_gi[, 1], nrow(srs_gi)+1) #Making space for stopping point at end of list
    locs1 <- cbind(nodes, rbind(sp[i, 1:2], srs_gi[, 2:3])) #adding the stopping point
    
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

      #The below code is just preparation to plot the TSP route.

      data1 <- locs1 %>%
        mutate(
          id_order = order(as.integer(tour))
        )

      data1 <- rbind(data1, data1[1,])
      data1[nrow(data1), ncol(data1)] <- nrow(data1)

    }else{
      distances[i,1] <- dist_mat*2    #if only one house sampled, walking distance is to the stopping point and back.
    }

    
    
    # # Plot a map with the data and overlay the optimal path
    # data1 %>%
    #   arrange(id_order) %>%
    #   leaflet() %>%
    #   addTiles() %>%
    #   addCircleMarkers(
    #     ~x,
    #     ~y,
    #     fillColor = ~ifelse(id_order == max(id_order), "yellow", "red"),
    #     fillOpacity = 0.6,
    #     stroke = F
    #   ) %>%
    #   addPolylines(~x, ~y)
    #
    # data1 %>%
    #   arrange(id_order) %>%
    #   leaflet() %>%
    #   addProviderTiles(providers$Esri.WorldImagery) %>%
    #   addCircleMarkers(
    #     ~x,
    #     ~y,
    #     fillColor = ~ifelse(id_order == max(id_order), "yellow", "red"),
    #     fillOpacity = 0.6,
    #     stroke = F
    #   ) %>%
    #   addPolylines(~x, ~y)
    #
    # if(mst == T){
    #   nodes_mst <- locs1[, 1]
    #   arcs <- gtools::combinations(max(nodes_mst), 2, nodes_mst)
    #   d_mst <- round(matrix(dist_mat, length(dist_mat), 1), 1)
    #   arcs <- cbind(arcs, d_mst)
    #   mst_i <- getMinimumSpanningTree(nodes_mst, arcs, algorithm = "Prim",
    #                                   show.data = F, show.graph = F)
    #   distances[i,2] <- mst_i$weight
    # }
    
  }
  boot_dists[iter,] <- colSums(distances, na.rm = T)/1000 #some entries will be NA, if a group did not appear in the sample.
  coverage[iter,] <- c((sum(srs$ndogs)/total_dogs), (nrow(srs)/n)) #% Dogs in first column, % houses in second
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
  nrow(srs)
  nsamp
}

if (save_boot_dists_in_sample_script == T){
  saveRDS(boot_dists, file = paste0("location/where/you/want/to/save/walking/distance/output"))
}