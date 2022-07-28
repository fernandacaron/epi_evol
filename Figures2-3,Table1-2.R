rm(list = ls())

setwd("Documents/lab/epi_evol")

library(phytools)
library(data.table)
library(picante)
library(caper)
library(geiger)
library(viridis)

## Loading previously saved phylogeny
tr <- read.tree("data/plants_tacted_100_4440sp.tre")

## Loading data
header <- read.csv("data/data_epiphyte_host.csv", header = T, nrow = 1, 
                   check.names = F)
dat <- fread("data/data_epiphyte_host.csv", skip = 1, header = F)
header <- iconv(colnames(header), to = "ASCII", sub = "")
setnames(dat, header)

## Merging lines with the same host name 
dat2 <- aggregate(dat[, 10:length(colnames(dat))], list(dat$Host_name), sum)
colnames(dat2)[1] <- "Host_name" 

## Getting the mean latitude and DBH for all host species
lat <- aggregate(abs(dat$Latitude_corr), list(dat$Host_name), mean)
colnames(lat) <- c("Host_name", "Lat")

dbh <- aggregate(as.numeric(dat$Host_dbh), list(dat$Host_name), mean)
colnames(dbh) <- c("Host_name", "Host_dbh")

## Merging columns with the same epiphytic names
for (i in 1:nrow(dat2)){
	dat2[i, "Dyssochroma_viridiflora"] <- 
									sum(dat2[i, "Dyssochroma_viridiflora"],
									    dat2[i, "Dyssochroma_viridiflora.1"])
	dat2[i, "Campyloneurum_amphostenon"] <- 
									sum(dat2[i, "Campyloneurum_amphostenon"],
									    dat2[i, "Campyloneurum_amphostenon.1"])
	dat2[i, "Lepanthes_VGJ694"] <- sum(dat2[i, "Lepanthes_VGJ694"],
									   dat2[i, "Lepanthes_VGJ694.1"])
}
drop <- c("Dyssochroma_viridiflora.1", "Campyloneurum_amphostenon.1", 
          "Lepanthes_VGJ694.1")
dat2 <- dat2[, -which(colnames(dat2) %in% drop)]

## Replacing all numbers higher than 1 with 1
dat2[, 2:dim(dat2)[2]] <- apply(dat2[, 2:dim(dat2)[2]], 2, as.numeric)

dat2[, 2:dim(dat2)[2]] <- apply(dat2[, 2:dim(dat2)[2]], 2, function(x)
                                ifelse(x>1, 1, x))

## Getting host and epiphytes names
host <- dat2$Host_name
epi <- colnames(dat2)[2:length(colnames(dat2))]

## Making a matrix with epiphytes in columns and hosts in rows
sample <- dat2[, 2:length(colnames(dat2))]
rownames(sample) <- dat2$Host_name

## Calculating SEH for all species, selecting only host plant values and calculating 
## MPD of epiphytes that are present in each host
seh_mpd <- list()
for(i in 1:100) {
	seh <- evol.distinct(tr[[i]], type = "equal.splits")
	seh_host <- seh[seh$Species %in% host, ]
	seh_host[, 3] <- NA

	mpd <- mpd(sample, cophenetic(tr[[i]]))
	names(mpd) <- rownames(sample)

	for(j in 1:nrow(seh_host)){
		seh_host[j, 3] <- mpd[names(mpd) == seh_host[j, 1]]
	}
	colnames(seh_host) <- c("Species", "SEH", "MPD")

	seh_mpd[[i]] <- seh_host
}

## Calculating SEH for each host and the number of epiphyte species of that host
seh_num <- list()
for(i in 1:100) {
	seh <- evol.distinct(tr[[i]], type = "equal.splits")
	seh_host <- seh[seh$Species %in% host, ]
	seh_host[, 3] <- NA

	for(j in 1:nrow(seh_host)){
		seh_host[j, 3] <- ifelse(sum(as.numeric(sample[rownames(sample) == 
		                                        seh_host[j, 1], ])) == 0,
								 NA,
								 sum(as.numeric(sample[rownames(sample) == 
		                                        seh_host[j, 1], ])))
	}
	colnames(seh_host) <- c("Species", "SEH", "Num_spp")

	seh_num[[i]] <- seh_host
	seh_num[[i]] <- seh_num[[i]][seh_num[[i]]$Num_spp > 1, ]
	seh_num[[i]] <- seh_num[[i]][complete.cases(seh_num[[i]][, 1]), ]
}

## Getting host latitude and MPD of the epiphytes on each host
lat_mpd <- list()
	for(i in 1:100) {
		lat_host <- as.data.frame(matrix(nrow = nrow(lat), ncol = 3))
		lat_host[, 2] <- lat$Lat[lat$Host_name %in% host]
		lat_host[, 1] <- lat$Host_name[lat$Host_name %in% host]
		lat_host[, 3] <- NA

	mpd <- mpd(sample, cophenetic(tr[[i]]))
	names(mpd) <- rownames(sample)

	for(j in 1:nrow(lat_host)){
		lat_host[j, 3] <- mpd[names(mpd) == lat_host[j, 1]]
	}
	colnames(lat_host) <- c("Species", "Lat", "MPD")

	lat_mpd[[i]] <- lat_host
}

## Getting DBH from host MPD from its epiphytes
dbh_mpd <- list()
for(i in 1:100) {
	dbh_host <- as.data.frame(matrix(nrow = nrow(dbh), ncol = 3))
	dbh_host[, 2] <- dbh$Host_dbh[dbh$Host_name %in% host]
	dbh_host[, 1] <- dbh$Host_name[dbh$Host_name %in% host]
	dbh_host[, 3] <- NA

	mpd <- mpd(sample, cophenetic(tr[[i]]))
	names(mpd) <- rownames(sample)

	for(j in 1:nrow(dbh_host)){
		dbh_host[j, 3] <- mpd[names(mpd) == dbh_host[j, 1]]
	}
	colnames(dbh_host) <- c("Species", "Host_dbh", "MPD")

	dbh_mpd[[i]] <- dbh_host
}

## Figure 2
cols <- c(rgb(52/255, 35/255, 70/255, 0.7),
          rgb(64/255, 73/255, 142/255, 0.7),
          rgb(53/255, 123/255, 162/255, 0.7),
          rgb(56/255, 170/255, 172/255, 0.7),
          rgb(120/255, 214/255, 174/255, 0.7))

col <- c(rgb(119/255, 119/255, 118/255, 0.5))

pdf("figures/Figure2.pdf", width = 14, height = 5)

layout(matrix(1:3, ncol = 3, byrow = T))

par(mar = c(5, 5, 4, 4))

hist(seh_mpd[[1]]$SEH, col = cols[2], main = "", xlab = "SEH (host)",
     border = mako(7)[3], breaks = 21, cex.lab = 1.9, cex.axis = 1.6)
legend("topright", legend  =  "A", bty  =  'n', cex = 3)

hist(seh_num[[1]]$Num_spp, col = cols[4], main = "",  ylab = "",
     xlab = "Number of epiphyte species per host species", border = mako(7)[5],
     breaks = 23, cex.lab = 1.9, cex.axis = 1.6)
legend("topright", legend  =  "B", bty  =  'n', cex = 3)

hist(seh_mpd[[1]]$MPD, col = cols[5], main = "", xlab = "MPD (epiphytes)",
     border = mako(7)[6], breaks = 32, ylab = "", cex.lab = 1.9, cex.axis = 1.6)
legend("topright", legend  =  "C", bty  =  'n', cex = 3)

dev.off()

## Figure 3
pdf("figures/Figure3.pdf", width = 9, height = 7)

layout(matrix(1:4, ncol = 2, byrow = T))

par(mar = c(4, 6, 1, 3))

plot(seh_mpd[[1]]$MPD ~ seh_mpd[[1]]$SEH, pch = 16, col = col,
     ylab = "MPD (epiphytes)", xlab = "log SEH (host)", log = "x")
legend("topright", legend  =  "A", bty  =  'n', cex = 1.5)

plot(seh_num[[1]]$Num_spp ~ seh_num[[1]]$SEH, pch = 16, col = col,
     ylab = "log Number of epiphytes species\nper host species", 
     xlab = "log SEH (host)", log = "xy")
legend("topright", legend  =  "B", bty  =  'n', cex = 1.5)

plot(lat_mpd[[1]]$MPD ~ lat_mpd[[1]]$Lat, pch = 16, col = col,
     ylab = "MPD (epiphytes)", xlab = "Latitude (host)")
legend("topright", legend  =  "C", bty  =  'n', cex = 1.5)


plot(dbh_mpd[[1]]$MPD ~ dbh_mpd[[1]]$Host_dbh, pch = 16, col = col,
     ylab = "MPD (epiphytes)", xlab = "log DBH (host)", log = "x")
legend("topright", legend  =  "D", bty  =  'n', cex = 1.5)

dev.off()

################

## PGLS and phylogenetic signal analyses

## PGLS (42.43056 mins)

## 1. MPD ~ log(SEH)

res1 <- matrix(ncol = 6, nrow = 100)
colnames(res1) <- c("tree", "slope", "SE_slope", "t_slope", "P_slope",
                    "R_squared")

for (i in 1:100) {
	tr.pruned <- treedata(tr[[i]], sample, warnings = F)$phy
	comp_data <- comparative.data(tr.pruned, seh_mpd[[i]], names.col = Species)
	p1 <- pgls(MPD ~ log(SEH), data = comp_data)

	res1[i, ] <- cbind(i, 
	                  summary(p1)$coefficients[2,1],
	                  summary(p1)$coefficients[2,2],
	                  summary(p1)$coefficients[2,3],
	                  summary(p1)$coefficients[2,4],
	                  summary(p1)$r.squared)
}

## Getting median and range of estimates across topologies
stats1 <- matrix(ncol = 1, nrow = 5)
rownames(stats1)<-c("slope", "SE_slope", "t_slope", "P_slope", "R_squared")

stats1[1, 1] <- paste0(round(median(res1[, 2]), 3), " (", 
                       round(range(res1[, 2])[1], 3), "-",
                       round(range(res1[, 2])[2], 3), ")")
stats1[2, 1] <- paste0(round(median(res1[, 3]), 3), " (", 
                       round(range(res1[, 3])[1], 3), "-",
                       round(range(res1[, 3])[2], 3), ")")
stats1[3, 1] <- paste0(round(median(res1[, 4]), 3), " (", 
                       round(range(res1[, 4])[1], 3), "-",
                       round(range(res1[, 4])[2], 3), ")")
stats1[4, 1] <- paste0(round(median(res1[, 5]), 3), " (", 
                       round(range(res1[, 5])[1], 3), "-",
                       round(range(res1[, 5])[2], 3), ")")
stats1[5, 1] <- paste0(round(median(res1[, 6]), 3), " (", 
                       round(range(res1[, 6])[1], 3), "-",
                       round(range(res1[, 6])[2], 3), ")")

## 2. log(Num_spp) ~ log(SEH)
res2 <- matrix(ncol = 6, nrow = 100)
colnames(res2)<-c("tree", "slope", "SE_slope", "t_slope", "P_slope", 
                     "R_squared")

for (i in 1:100) {
	tr.pruned <- treedata(tr[[i]], sample, warnings = F)$phy
	comp_data <- comparative.data(tr.pruned, seh_num[[i]], names.col = Species)
	p2 <- pgls(log(Num_spp) ~ log(SEH), data = comp_data)

	res2[i, ] <- cbind(i, 
	                  summary(p2)$coefficients[2,1],
	                  summary(p2)$coefficients[2,2],
	                  summary(p2)$coefficients[2,3],
	                  summary(p2)$coefficients[2,4],
	                  summary(p2)$r.squared)
}

stats2 <- matrix(ncol = 1, nrow = 5)
rownames(stats2)<-c("slope", "SE_slope", "t_slope", "P_slope", "R_squared")

stats2[1, 1] <- paste0(round(median(res2[, 2]), 3), " (", 
                       round(range(res2[, 2])[1], 3), "-",
                       round(range(res2[, 2])[2], 3), ")")
stats2[2, 1] <- paste0(round(median(res2[, 3]), 3), " (", 
                       round(range(res2[, 3])[1], 3), "-",
                       round(range(res2[, 3])[2], 3), ")")
stats2[3, 1] <- paste0(round(median(res2[, 4]), 3), " (", 
                       round(range(res2[, 4])[1], 3), "-",
                       round(range(res2[, 4])[2], 3), ")")
stats2[4, 1] <- paste0(round(median(res2[, 5]), 3), " (", 
                       round(range(res2[, 5])[1], 3), "-",
                       round(range(res2[, 5])[2], 3), ")")
stats2[5, 1] <- paste0(round(median(res2[, 6]), 3), " (", 
                       round(range(res2[, 6])[1], 3), "-",
                       round(range(res2[, 6])[2], 3), ")")

## 3. MPD ~ Lat
res3 <- matrix(ncol = 6, nrow = 100)
colnames(res3)<-c("tree", "slope", "SE_slope", "t_slope", "P_slope", 
                     "R_squared")

for (i in 1:100) {
	tr.pruned <- treedata(tr[[i]], sample, warnings = F)$phy
	comp_data <- comparative.data(tr.pruned, lat_mpd[[i]], names.col = Species)
	p3 <- pgls(MPD ~ abs(Lat), data = comp_data)

	res3[i, ] <- cbind(i, 
	                  summary(p3)$coefficients[2,1],
	                  summary(p3)$coefficients[2,2],
	                  summary(p3)$coefficients[2,3],
	                  summary(p3)$coefficients[2,4],
	                  summary(p3)$r.squared)
}

stats3 <- matrix(ncol = 1, nrow = 5)
rownames(stats3)<-c("slope", "SE_slope", "t_slope", "P_slope", "R_squared")

stats3[1, 1] <- paste0(round(median(res3[, 2]), 3), " (", 
                       round(range(res3[, 2])[1], 3), "-",
                       round(range(res3[, 2])[2], 3), ")")
stats3[2, 1] <- paste0(round(median(res3[, 3]), 3), " (", 
                       round(range(res3[, 3])[1], 3), "-",
                       round(range(res3[, 3])[2], 3), ")")
stats3[3, 1] <- paste0(round(median(res3[, 4]), 3), " (", 
                       round(range(res3[, 4])[1], 3), "-",
                       round(range(res3[, 4])[2], 3), ")")
stats3[4, 1] <- paste0(round(median(res3[, 5]), 3), " (", 
                       round(range(res3[, 5])[1], 3), "-",
                       round(range(res3[, 5])[2], 3), ")")
stats3[5, 1] <- paste0(round(median(res3[, 6]), 3), " (", 
                       round(range(res3[, 6])[1], 3), "-",
                       round(range(res3[, 6])[2], 3), ")")

## 4. MPD ~ DBH
res4 <- matrix(ncol = 6, nrow = 100)
colnames(res4)<-c("tree", "slope", "SE_slope", "t_slope", "P_slope", 
                     "R_squared")

for (i in 1:100) {
	tr.pruned <- treedata(tr[[i]], sample, warnings = F)$phy
	comp_data <- comparative.data(tr.pruned, dbh_mpd[[i]], names.col = Species)
	p4 <- pgls(MPD ~ log(Host_dbh), data = comp_data)

	res4[i, ] <- cbind(i, 
	                  summary(p4)$coefficients[2,1],
	                  summary(p4)$coefficients[2,2],
	                  summary(p4)$coefficients[2,3],
	                  summary(p4)$coefficients[2,4],
	                  summary(p4)$r.squared)
}

stats4 <- matrix(ncol = 1, nrow = 5)
rownames(stats4)<-c("slope", "SE_slope", "t_slope", "P_slope", "R_squared")

stats4[1, 1] <- paste0(round(median(res4[, 2]), 3), " (", 
                       round(range(res4[, 2])[1], 3), "-",
                       round(range(res4[, 2])[2], 3), ")")
stats4[2, 1] <- paste0(round(median(res4[, 3]), 3), " (", 
                       round(range(res4[, 3])[1], 3), "-",
                       round(range(res4[, 3])[2], 3), ")")
stats4[3, 1] <- paste0(round(median(res4[, 4]), 3), " (", 
                       round(range(res4[, 4])[1], 3), "-",
                       round(range(res4[, 4])[2], 3), ")")
stats4[4, 1] <- paste0(round(median(res4[, 5]), 3), " (", 
                       round(range(res4[, 5])[1], 3), "-",
                       round(range(res4[, 5])[2], 3), ")")
stats4[5, 1] <- paste0(round(median(res4[, 6]), 3), " (", 
                       round(range(res4[, 6])[1], 3), "-",
                       round(range(res4[, 6])[2], 3), ")")

stats <- cbind(stats1, stats2, stats3, stats4)
colnames(stats) <- c("MPD ~ log(SEH)", "log(Num_spp) ~ log(SEH)",
                     "MPD ~ Latitude", "MPD ~ DBH")

write.csv(stats, "tables/Table1_unformatted.csv")

## Phylogenetic signal

## MPD
physig_mpd <- matrix(ncol = 6, nrow = 100)
for(i in 1:100) {
	mpd <- mpd(sample, cophenetic(tr[[i]]))
	names(mpd) <- rownames(sample)
	mpd <- mpd[complete.cases(mpd)]

	tr.pruned <- treedata(tr[[i]], mpd, warnings = FALSE)$phy
	physig_l <- phylosig(tr.pruned, mpd, method = "lambda", test = TRUE)
	physig_k <- phylosig(tr.pruned, mpd, method = "K", test = TRUE)
	
	physig_mpd[i, ] <- cbind(physig_l$lambda, physig_l$logL, physig_l$logL0, 
	                         physig_l$P, physig_k$K, physig_k$P)
}

## Latitude
physig_lat <- matrix(ncol = 6, nrow = 100)
for(i in 1:100) {
	lat2 <- lat$Lat
	names(lat2) <- lat$Host_name

	tr.pruned <- treedata(tr[[i]], lat2, warnings = FALSE)$phy
	physig_l <- phylosig(tr.pruned, lat2, method = "lambda", test = TRUE)
	physig_k <- phylosig(tr.pruned, lat2, method = "K", test = TRUE)
	
	physig_lat[i, ] <- cbind(physig_l$lambda, physig_l$logL, physig_l$logL0, 
	                         physig_l$P, physig_k$K, physig_k$P)
}

## DBH
physig_dbh <- matrix(ncol = 6, nrow = 100)
for(i in 1:100) {
	dbh2 <- dbh$Host_dbh
	names(dbh2) <- dbh$Host_name
	dbh2 <- dbh2[complete.cases(dbh2)]

	tr.pruned <- treedata(tr[[i]], dbh2, warnings = FALSE)$phy
	physig_l <- phylosig(tr.pruned, dbh2, method = "lambda", test = TRUE)
	physig_k <- phylosig(tr.pruned, dbh2, method = "K", test = TRUE)
	
	physig_dbh[i, ] <- cbind(physig_l$lambda, physig_l$logL, physig_l$logL0, 
	                         physig_l$P, physig_k$K, physig_k$P)
}


physig <- matrix(ncol = 6, nrow = 3)
colnames(physig)<-c("Lambda", "logL", "logL0", "P", "K", "P_K")
rownames(physig)<-c("MPD", "Latitude", "DBH")

physig[1, 1] <- paste0(round(median(physig_mpd[, 1]), 3), " (", 
                       round(range(physig_mpd[, 1])[1], 3), "-",
                       round(range(physig_mpd[, 1])[2], 3), ")")
physig[1, 2] <- paste0(round(median(physig_mpd[, 2]), 3), " (", 
                       round(range(physig_mpd[, 2])[1], 3), "-",
                       round(range(physig_mpd[, 2])[2], 3), ")")
physig[1, 3] <- paste0(round(median(physig_mpd[, 3]), 3), " (", 
                       round(range(physig_mpd[, 3])[1], 3), "-",
                       round(range(physig_mpd[, 3])[2], 3), ")")
physig[1, 4] <- paste0(round(median(physig_mpd[, 4]), 3), " (", 
                       round(range(physig_mpd[, 4])[1], 3), "-",
                       round(range(physig_mpd[, 4])[2], 3), ")")
physig[1, 5] <- paste0(round(median(physig_mpd[, 5]), 3), " (", 
                       round(range(physig_mpd[, 5])[1], 3), "-",
                       round(range(physig_mpd[, 5])[2], 3), ")")
physig[1, 6] <- paste0(round(median(physig_mpd[, 6]), 3), " (", 
                       round(range(physig_mpd[, 6])[1], 3), "-",
                       round(range(physig_mpd[, 6])[2], 3), ")")

physig[2, 1] <- paste0(round(median(physig_lat[, 1]), 3), " (", 
                       round(range(physig_lat[, 1])[1], 3), "-",
                       round(range(physig_lat[, 1])[2], 3), ")")
physig[2, 2] <- paste0(round(median(physig_lat[, 2]), 3), " (", 
                       round(range(physig_lat[, 2])[1], 3), "-",
                       round(range(physig_lat[, 2])[2], 3), ")")
physig[2, 3] <- paste0(round(median(physig_lat[, 3]), 3), " (", 
                       round(range(physig_lat[, 3])[1], 3), "-",
                       round(range(physig_lat[, 3])[2], 3), ")")
physig[2, 4] <- paste0(round(median(physig_lat[, 4]), 3), " (", 
                       round(range(physig_lat[, 4])[1], 3), "-",
                       round(range(physig_lat[, 4])[2], 3), ")")
physig[2, 5] <- paste0(round(median(physig_lat[, 5]), 3), " (", 
                       round(range(physig_lat[, 5])[1], 3), "-",
                       round(range(physig_lat[, 5])[2], 3), ")")
physig[2, 6] <- paste0(round(median(physig_lat[, 6]), 3), " (", 
                       round(range(physig_lat[, 6])[1], 3), "-",
                       round(range(physig_lat[, 6])[2], 3), ")")

physig[3, 1] <- paste0(round(median(physig_dbh[, 1]), 3), " (", 
                       round(range(physig_dbh[, 1])[1], 3), "-",
                       round(range(physig_dbh[, 1])[2], 3), ")")
physig[3, 2] <- paste0(round(median(physig_dbh[, 2]), 3), " (", 
                       round(range(physig_dbh[, 2])[1], 3), "-",
                       round(range(physig_dbh[, 2])[2], 3), ")")
physig[3, 3] <- paste0(round(median(physig_dbh[, 3]), 3), " (", 
                       round(range(physig_dbh[, 3])[1], 3), "-",
                       round(range(physig_dbh[, 3])[2], 3), ")")
physig[3, 4] <- paste0(round(median(physig_dbh[, 4]), 3), " (", 
                       round(range(physig_dbh[, 4])[1], 3), "-",
                       round(range(physig_dbh[, 4])[2], 3), ")")
physig[3, 5] <- paste0(round(median(physig_dbh[, 5]), 3), " (", 
                       round(range(physig_dbh[, 5])[1], 3), "-",
                       round(range(physig_dbh[, 5])[2], 3), ")")
physig[3, 6] <- paste0(round(median(physig_dbh[, 6]), 3), " (", 
                       round(range(physig_dbh[, 6])[1], 3), "-",
                       round(range(physig_dbh[, 6])[2], 3), ")")

write.csv(physig, "tables/Table2_unformatted.csv")
