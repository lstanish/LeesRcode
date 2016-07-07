setwd("/Users/stanish/Documents/R")   # for Macs

library(vegan)
source("transposeOTU.r")      ## contains function t.otu()
source("Fill.Taxonomy.r")     ## ## contains function fill.taxonomy()
library(scatterplot3d)
source("big.cor.r")
library(ggplot2)

otu <- read.csv("OhioRiver_OTU_table_noCPs_UCsParsed_19350.csv")
otu <- read.csv("Ross_OTUtable_raw.csv")
env <- read.csv(file="OhioRiver_Metadata.csv", row.names=1)
chem <- read.csv("OhioRiver_Chemistry.csv", row.names=1)

alldat <- cbind(env,chem)
alldat.num <- alldat[,c(13:50)]

# Fill in empty taxonomy elements [Lee's function]
otu4 <- fill.taxonomy(otu)
# Transpose OTU table [Lee's function]
otu5 <- t.otu(otu4, s=2, e=68)
# Create taxonomy matrix 
taxon.key <- otu4[,c(1,69:75)]


#Transform OTU table
otu.pct <- matrix(0,nrow(otu5),ncol(otu5))
for (k in 1:nrow(otu5)){
  otu.pct[k,]<-otu5[k,]/sum(otu5[k,])
}


### Acidithiobacillus 
acidi <- which(taxon.key$Genus=="Acidithiobacillus")
acidi.pct <- otu.pct[,acidi]*100

### Firmicutes
firms <- which(taxon.key$Phylum=="Firmicutes")
firm.tot <- vector(mode="numeric", length=nrow(otu.pct))
for(i in 1:length(firm.tot)) {
  firm.tot[i] <- sum(otu.pct[i,firms])
}

#### Proteobacteria
prots <- which(taxon.key$Class=="Alphaproteobacteria")
prots.tot <- vector(mode="numeric", length=nrow(otu.pct))
for(i in 1:length(prots.tot)) {
  prots.tot[i] <- sum(otu.pct[i,prots])
}


### Create objects to hold sums for a and b proteos
aprots <- prots.tot
bprots <- prots.tot
Free.Cl <- alldat$Free.Cl
State <- alldat$State
plotdat <- data.frame(cbind(Free.Cl,"Tot.Cl"=alldat$Tot.Cl,
                            "Fl"=alldat$F.mgL, "Fluoridation"=alldat$Fluoridation,
                            "aprots"=aprots, "firms"=firm.tot,"acidi"=acidi))

cors.3 <- cor(alldat.num,firm.tot)

ggplot(plotdat, aes(x=Free.Cl,y=acidi.pct, colour=State) ) + 
  geom_point(size=3) +
   labs(x=expression("Free Cl"^"-"), y="Relative abundance") + 
  theme(panel.grid.minor.x=element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        legend.key=element_blank(),
        legend.background=element_rect(colour="black"),
        axis.text=element_text(colour="black"),
        axis.ticks=element_line(colour="black"),
        legend.position=c(0.1, 0.8) )



ggplot(plotdat, aes(Free.Cl, firm.tot ) ) +
  geom_point(size = 3) +
  stat_smooth() + labs(x=expression("Free Cl"^"-"))

ggplot(plotdat, aes(Free.Cl, firm.tot ) ) +
  geom_point(size = 3) +
  stat_smooth() + labs(x=expression("Free Cl"^"-"))
