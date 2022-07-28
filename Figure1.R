rm(list = ls())

setwd("Documents/lab/epi_evol")

library(phytools)
library(viridis)

## Loading the data
tr <- read.tree("data/tact/plants_tacted_100_6562sp.tre")
dat <- read.csv("data/data_epiphyte_host.csv")

## Getting the names of host plants that have data
host <- dat$Host_name
host <- host[complete.cases(host)]

## Getting the names of epiphytes that have data
epi <- colnames(dat)[10:length(colnames(dat))]
for (i in 1:94) {
	epi <- epi[! epi == paste0("NA..", i)]
}
epi <- epi[! epi == "NA."]

## Fixing some names
epi[epi == "Dyssochroma_viridiflora.1"] <- "Dyssochroma_viridiflora"
epi[epi == "Campyloneurum_amphostenon.1"] <- "Campyloneurum_amphostenon"
epi[epi == "Lepanthes_VGJ694.1"] <- "Lepanthes_VGJ694"

## Checking if all host and epiphytes names are in the phylogeny
setdiff(host, tr[[1]]$tip.label)
setdiff(epi, tr[[1]]$tip.label)

## Joining the names
all <- c(epi, host)
all <- names(table(all))

## Pruning the tips of the phylogeny that are not in the data
for (i in 1:100) {
	tr[[i]] <- keep.tip(tr[[i]], all)
}

## Saving the resulting topologies
for (i in 1:100) write.tree(tr[[i]], "data/plants_tacted_100_4440sp.tre", append = TRUE)

## Loading data with families and types (epiphyte or host) of each species
info <- read.csv("data/info_host_epiphyte.csv")

## Getting the type of each species in the phylogeny
types <- info$Type
names(types) <- info$Species
types <- types[names(types) %in% tr[[1]]$tip.label]

## Creating a character map object with a topology
mtree <- make.simmap(tr[[1]], types)

## Setting the colors of each type 
cols <- setNames(c(viridis(20)[17], viridis(20)[11]), c("Host", "Epiphyte"))

pdf("figures/Figure1.pdf", height = 9, width = 8)

plot(mtree, cols, type = "fan", ftype = "off", lwd = 1)

arc.cladelabels(text = "Orchidaceae", mark.node = F, cex = 0.8, 
                col = "black", lwd = 1, ln.offset = 1.05, lab.offset = 1.1, 
                node = findMRCA(mtree, 
                                c(info$Species[info$Family == "Orchidaceae"])))

arc.cladelabels(text = "Bromeliaceae", mark.node = F, cex = 0.8, 
                col = "black", lwd = 1, ln.offset = 1.05, lab.offset = 1.1, 
                node = findMRCA(mtree, 
                                c(info$Species[info$Family == "Bromeliaceae"])))

arc.cladelabels(text = "Ferns", mark.node = F, cex = 0.8, 
                col = "black", lwd = 1, ln.offset = 1.05, lab.offset = 1.1, 
                node = findMRCA(mtree, c(info$Species[info$Seed_plant_Fern == "Fern"])))

add.simmap.legend(colors = cols, prompt = FALSE, x = 0.9*par()$usr[1],
                  y = 0.8*par()$usr[3], fsize = 0.8)

dev.off()

