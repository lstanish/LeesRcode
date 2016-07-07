######### Ordination analysis of otu tables using non-metric MDS #####

setwd("/Users/stanish/Documents/R")   # for Macs
### Generified script to import an OTU table and prepare for downstream analysis
### LFS June 22 2016

library(vegan)
library(ggplot2)
library(grid)


### Load required function ###
t.otu <- function(x, s, e) {
  ## Written by Lee Stanish
  ## transpose OTU table for ordination analysis##
  ## x=OTU table; s=starting column; e=ending column
  x <- as.data.frame(x) 
  x1 <- t(x[,c(s:e)])      # transpose OTU counts
  colnames(x1) <- x[,1]     #replace OTU ids as column names 
  x1
}
### End required function ##


######    Start of code    ######

setwd("/Users/stanish/Documents/R")  # set to your dir
otu <- read.csv("OhioRiver_OTU_table_noCPs_UCsParsed_19350.csv")  # set to your file


## transpose OTU table using function t.otu 
otu.t <- t.otu(x=otu,s=2,e=68)  # adjust numbers as needed for your OTU table
taxonomy <- otu$Taxonomy   # check that #columns in otu.t=length of taxonomy

diversity(otu.t)
specnumber(otu.t)


###########           Analysis code (for later)       ####################

###Rarefy OTU table###
### A) Calculate total # seqs for each sample
tototu <- apply(otu.t, MARGIN=1, FUN=sum)
min(tototu)			# rarefy to this value

### B) Rarefy
oturare <- rrarefy(otu.t, min(tototu))

### C) Generate rarefaction curves (test)
rarecurve(otu.t, step=50,sample=min(tototu))

### D) Remove OTUs with 0 abundance
otumax <- apply(oturare, 2, max)
otu.1 <- oturare[, otumax >= 1]
taxonomy.1 <- droplevels(taxonomy[otumax >= 1])

### E) Create key for OTU taxonomies and OTU
otu.export <- as.matrix(t(otu.1))
OTU <- 1:nrow(otu.export)

### F) Finalize OTU table
otu.final <- cbind(OTU,otu.export)

### G) Export OTU table
write.csv(otu.final, paste("OhioRiver_OTU_table_noCPs_UCsParsed_",min(tototu),".csv", sep=""))

#Transform rarefied OTU table
otu.pct <- decostand(otu.1, method="total")
otu.hell <- decostand(otu.1, method="hellinger")  # Hellinger transformation


######## nMDS  #####
mds <- metaMDS(otu.hell, k=3, distance="bray", trymax=20, autotransform=FALSE)

# Shepard plot to compare actual dissimilarities to computed dissimilarities
diss <- vegdist(otu.hell, distance="bray")
stressplot(mds, diss) ## pretty nice fit

mds.sc <- scores(mds, display="sites", choices=c(1:3))
mds.spec <- scores(mds, display="species", choices=c(1:3))


########     Plot results (not generified)    ########

###All samples#
plot(mds, type="p", choices=c(1,3), cex.lab=1.5,cex.axis=1.5,font=1.5,font.axis=2,font.lab=1.5)

mds.sc <- data.frame(mds.sc)


ggplot(mds.sc, aes(NMDS1, NMDS2)) + 
  geom_point(data=mds.sc,aes(NMDS1, NMDS2, size=alldat$Tot.Cl, shape=alldat$ChlType), colour="darkgrey") +          
  geom_segment(data=fit.df, aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow=arrow(length=unit(0.2,"cm")), colour="black", inherit_aes=FALSE) +
  geom_text(data=fit.df, aes(x=NMDS1*1.1, y=NMDS2*1.1, label=species), size=4.2)+ coord_fixed() +
  guides(shape=guide_legend("Disinfectant"), 
         size=guide_legend(expression(paste("Total Cl-", ", mg ", l^-1, sep=""))) ) + 
  theme_bw() +
  theme(axis.title=element_text(size=14), 
        axis.text=element_text(size=11),
        legend.key=element_blank(),
        legend.key=element_rect(colour="white"), 
        legend.text=element_text(size=11), 
        legend.title=element_text(size=11))
#+ ggtitle("NMDS of bacterial communities (3D stress=0.12)")


######## Select more abundant, 5%, OTUs to plot ######
spec1max <- apply(otu.1,2,max)
mdsspec5 <- mds$species[spec1max>968,]   ##nMDS species coordinates for OTUs

spec5IDs <- taxon.key1[spec1max>968,]

spec5IDs$Taxonomy <- factor(spec5IDs$Taxonomy)   ## correct levels for factors
spec5IDs$Domain <- factor(spec5IDs$Domain)
spec5IDs$Phylum <- factor(spec5IDs$Phylum)
spec5IDs$Class <- factor(spec5IDs$Class)
spec5IDs$Order <- factor(spec5IDs$Order)
spec5IDs$Family <- factor(spec5IDs$Family)
spec5IDs$Genus <- factor(spec5IDs$Genus)

### Replace proteo genera with proteo subphylum ###
pal <- c("#76EEC6", "#0000EE", "#8A2BE2", "#8B2323", "#5F9EA0", 
         "#458B00", "#EE7621","#EE6A50", "#6495ED", "#EEAD0E",
         "#8B008B", "#FF7F00", "#8B0000", "#9BCD9B", "#2F4F4F")
phy <- which(spec5IDs$Phylum!="Proteobacteria")
cla <- which(spec5IDs$Phylum=="Proteobacteria")
plot.tax <- vector(mode="character", length=nrow(spec5IDs))
plot.tax[phy] <- as.vector(spec5IDs$Phylum[phy])
plot.tax[cla] <- as.vector(spec5IDs$Class[cla])
plot.tax <- factor(plot.tax)
col.int <- as.integer(plot.tax)
plot.cols <- pal[col.int]

mdsspec5 <- data.frame(mdsspec5)
mds.sc <- data.frame(mds.sc)

ggplot(aes(NMDS1, NMDS2), data=mds.sc) + 
  geom_point(colour="darkgrey", shape=21, size=0.9) +
  

### GGPLOT WITH PROTEO SUBPHYLA NAMES ###
  q <- ggplot(data=mdsspec5, aes(x=MDS1, y=MDS2))
r <- q + geom_point(size=4, pch=21) +
  aes(fill=plot.tax)
r + guides(size=FALSE) +
  guides(fill=guide_legend(title="Phylum")) +
  labs(x="NMDS1", y="NMDS2", size=20) +
  theme_bw()+ 
  theme(axis.title=element_text(size=14), 
        axis.text=element_text(size=11),
        legend.key=element_rect(colour="white"), 
        legend.text=element_text(size=11), 
        legend.title=element_text(size=11)) 



##### CODE FOR REMOVING GRID LINES: KEEP CODE#####
#theme(panel.grid.major=element_blank(), 
#  panel.grid.minor=element_blank(),
#  legend.key=element_rect(colour="white")) + 	  


######## Change plotting colors for taxa labels #####
colors <- c("black","red","green3","blue","magenta","darkgrey","cyan")
palette(colors)


##save to pdf file##
pdf(file="OR_3DnMDS_Taxa10percent.pdf")
plot(mds.sc[,1], mds.sc[,2], pch=20,col="gray", cex=0.8, xlim=c(-1.5,1.5),
     ylim=c(-1,1.3))
text(mdsspec5[,1], mdsspec5[,2], labels=spec5IDs$Class, cex=0.9, 
     col=as.integer(spec5IDs$Phylum))
mtext("Bacterial taxa >= 10% abundance classified to family level")
legend("topleft", bty="n", legend=c(levels(spec5IDs$Phylum)), pt.bg=c(1:7),
       pch=22)
dev.off()



#########  -----------MULTIVARIATE ANALYSES ------------- #########
#Step 1: remove groups with only one sample
to.rm <- c(which(alldat$Municipality=="Pikeville"),
           which(alldat$Municipality=="Morgantown") )
alldat.rm <- alldat[-to.rm,]
otu.hell.rm <- otu.hell[-to.rm,]

diss <- vegdist(otu.hell.rm, method="euclidean")
x <- betadisper(diss,group=alldat.rm$Municipality)
permutest(x)  
anova(x)
boxplot(x)
TukeyHSD(x)
#Removing those 2 samples greatly improves stability using either BC or Euc dist

##### PERMANOVA analysis for group membership    #####
adonis(diss~State, data=alldat.rm)

#########  MANTEL TESTS  #######
envdis <- vegdist(alldat.rm$Free.Cl,method="euclidean")
mantel(envdis,diss,strata=alldat.rm$Municipality)

####  ANOSIM for group membership   ####
mod <- anosim(otu.hell[-52,], grouping=env$ChlType[-52],distance="bray",
              strata=env$State[-52])