# This script is sourced by the data_preparation.R script.
# Master script should have been run first, to import packages.

houses_knn <- data.frame(cbind(x, y))

n <- nrow(houses_knn)
dis <- as.matrix(dist(houses_knn))


f_x <- c()
for(i in 1:n){
  focal <- sort(dis[, i])[1:k]
  f_x[i] <- 1/focal[k]
}

houses_knn[,3] <- f_x
colnames(houses_knn)[3] <- "f_x"

### Find the k local maximums (clusters)
assignments <- matrix(NA, nrow = n, ncol = 2)
for (i in 1:n){
  d <- data.frame(cbind(dis[, i], f_x, 1:n))
  colnames(d) <- c("dist", "f_x", "index")
  d <- d %>% arrange(dist)
  focal <- d[1:k,]
  lm <- which(focal$f_x == max(focal$f_x))[1] #Sometimes, there can be two observations with he same density.
  lm_index <- focal$index[lm]
  crit <- 0  #criteria
  if (lm_index == i){
    crit <- 1
  }
  while (crit == 0){
    d <- data.frame(cbind(dis[, lm_index], f_x, 1:n))
    colnames(d) <- c("dist", "f_x", "index")
    d <- d %>% arrange(dist)
    focal <- d[1:k,]
    lm <- which(focal$f_x == max(focal$f_x))[1]
    lm_index_new <- focal$index[lm]
    if (lm_index == lm_index_new){
      crit <- 1
    }
    lm_index <- lm_index_new
  }
  assignments[i,] <- c(i, lm_index)
}

sp <- houses_knn[unique(assignments[,2]),]
# The sp variable is then used further in the data_preparation.R script.

if (remove_close_points == T){
  #Remove stopping points that are within a 200m radius of one another
  rm_dist <- 200 #Removal distance, any stopping point within this distance ought to be removed.
  test <- as.matrix(dist(sp[, 1:2]))*metres
  locs <- which(test < rm_dist & test > 0)
  column <- unique(locs%%nrow(test))
  for (index in colnames(test[,column])){
    culprit <- which(rownames(sp) == index) #Find the row index of the stopping point that is within 200m of another.
    sp <- sp[-culprit,]
  }
}