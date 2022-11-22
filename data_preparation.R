###### Data preparation #######

#IMPORTING AND PREPARING THE VILLAGE DATA
file_dd <- "your/data/path" #Once you have the data from Katie Hampson, store it in you desider folder, and source it here.
dd <- openxlsx::read.xlsx(file_dd) # Read the dog data file into your workspace

#For the below script to work, your dd object (dataframe) should be structured as follows:
# Column 1: Latitude
# Column 2: Longitude
# Column 3: Village Name
# Column 4: Number of dogs > 3 months of age
# Column 5: Number of Vaccinated dogs > 3 months of age
# Column 6: Number of dogs < 3 months of age
# Column 7: Number of Vaccinated dogs < 3 months of age

metres <- 111139 #Constant used to convert polar coordinates to metres.
#adjusted_sps = openxlsx::read.xlsx("Adjusted stopping points.xlsx")

#PREPARING VILLAGE
d <- subset(dd, Village==curr_village)
d$total_dogs <- d$`Dogs.(>.3.months)` + d$`Pups.(<.3.months)`
if (exclude_nodoghouses == T){
  d <- d %>% filter(d$total_dogs > 0 )
}
d$Village <- NULL #removing village column (we know we are in Machochwe)
total_adults <- sum(d$`Dogs.(>.3.months)`)
total_pups <- sum(d$`Pups.(<.3.months)`)
total_dogs <- total_adults + total_pups
x <- d[, 2] #Longitude
y <- d[, 1] #Latitude
nvacs <- d$`Vaccinated.dogs.(>.3.months)` + d$`Vaccinated.pups.(<.3.months)`
ndogs <- d$`Dogs.(>.3.months)` + d$`Pups.(<.3.months)`

W <- owin(xrange = c(min(x, na.rm = T), max(x, na.rm = T)), yrange = c(min(y, na.rm = T), max(y, na.rm = T)))
village_pp <- ppp(x, y, window = W)
ch <- convexhull(village_pp) #Creating a convex hull window
ch_sp <- as(ch, "SpatialPolygons") #converting ch to Spatial Polygon for sp package for some sampling schemes.
village <- ppp(x, y, ch)
if (miplot == 1){
  plot(village, main = paste0(curr_village))
}

if (sp_algorithm == "knn"){
  k <- theta
  remove_close_points <- F
  source("path_stop_stopping_point_algorithms/kNN mode-seeking stopping points.R") # Modify this path based on where you stored the SP algorithms

}else if (sp_algorithm == "kmeans"){
  k <- theta
  source("path_stop_stopping_point_algorithms/KMeans stopping points.R") # Modify this path based on where you stored the SP algorithms
}

#Preparing data to determine houses around stopping point
n <- length(x)
houses <- data.frame(idx = 1:n, cbind(x, y), ndogs, nvacs)
houses_and_sps <- data.frame(rbind(houses[, 2:3], sp[, 1:2])) #house locations stacked on top of stopping points

dis2 <- as.matrix(dist(houses_and_sps))[1:n, -(1:n)] #Only the distance of the stopping points to the nearest houses.
gh <- NULL #grouped houses - gh
for (i in seq_len(ncol(dis2))){
  coli <- data.frame(inx = seq(1, nrow(dis2)), dist = dis2[, i]*metres)
  coli <- coli[order(coli$dist),]
  basin <- coli[which(coli$dist < walking_radius),]
  if (miplot == 1){
    plot(houses[,2], houses[, 3])
    points(houses[basin[,1],2], houses[basin[,1], 3], col = "green", pch = 19)
    points(houses_and_sps[n+i,1], houses_and_sps[n+i,2], col="red", pch = 19)
  }
  group <- data.frame(houses[basin[, 1],], grp_num = rep(i, nrow(basin)))
  gh <- rbind(gh, group)
}

ghu <- gh[-which(duplicated(gh$idx) == TRUE),] #ghu - grouped houses UNIQUE (we do not want the same house chosen for more than one stopping point)
rownames(ghu) <- ghu[, 1]
if (dim(ghu)[1] == 0){ghu <- gh} #this code is important! If there are not duplicates, then ghu will be 0, thus we must set ghu to gh then.

if (miplot == 1){
  points(houses[,2], houses[,3])
  plot(houses[,2], houses[,3])
  points(gh[,2], gh[,3], col = gh[,6], pch=19)
  points(ghu[,2], ghu[,3], col = ghu[,6], pch = 19)
  points(houses_and_sps[n:nrow(houses_and_sps),1], houses_and_sps[n:nrow(houses_and_sps),2], col="red", pch = 19)
}
#Here, ghu is the subset of Machochwe which we obtained from a radius around each stopping point.
if (miplot == 1){
  plot(houses[,3], houses[,2])
  points(ghu[,3], ghu[,2], col = "green", pch = 19)
  points(houses_and_sps[n:nrow(houses_and_sps),1], houses_and_sps[n:nrow(houses_and_sps),2], col="red", pch = 19)
}


dogs_avail_prop <- sum(ghu$ndogs)/total_dogs
houses_avail_prop <- nrow(ghu)/n
dogs_avail_prop
houses_avail_prop

if (miplot == 1){
  plot(village)
  points(ghu[,2:3], col = "blue", pch = 19)
}