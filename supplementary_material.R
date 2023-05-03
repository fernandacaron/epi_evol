#Demonstration of how incomplete taxon sampling affects estimates of MPD
library(phytools)
library(picante)
set.seed(12345)

tr1<-pbtree(n=100)
tr2<-drop.tip(tr1, sample(tr1$tip.label, 0.8*Ntip(tr1)))
data<-matrix(rbinom(2000,1,0.5), ncol=20)
colnames(data)<-tr2$tip.label
plot(mpd(data, cophenetic(tr1)), mpd(data, cophenetic(tr2)))


#Demonstration of how incomplete taxon sampling affects estimates of SEH
library(phytools)
library(picante)
set.seed(12345)

tr1<-pbtree(n=10000)
tr2<-drop.tip(tr1, sample(tr1$tip.label, 0.8*Ntip(tr1)))

x1<-evol.distinct(tr1)
x2<-evol.distinct(tr2)

x3 <- x1[x1$Species %in% x2$Species, ]

cor.test(x2$w, x3$w, method="spearman")
	# Spearman's rank correlation rho

# data:  x2$w and x3$w
# S = 710282731, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
      # rho 
# 0.4672878 