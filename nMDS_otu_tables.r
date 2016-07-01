######### Ordination analysis of otu tables using non-metric MDS #####

setwd("yourdir") 

library(vegan)
library(ggplot2)
library(grid)
source("Fill.Taxonomy.r")    ## contains functions fill.taxonomy() and t.otu()
library(scatterplot3d)
source("big.cor1.r")

otu <- read.csv("your_otu_table.csv")
env <- read.csv(file="your_metadata_table.csv", row.names=1)


## transpose OTU table, if needed ##
use function t.otu(s=startingCol, e=endingCol)
otu2 <- t.otu(otu,s=yourStart,e=yourEnd)

#Remove samples with min tots < yourMin.

oturm <- vector("numeric",0)
for(i in 1:nrow(otu2)) {
	if (sum(otu2[i,])<1000) {
		oturm <- c(oturm,i)
	}
	otu2a <- otu2[-oturm,]
}
 env2a <- env[-oturm,]


###Rarefy OTU table###
	### A) Calculate total # seqs for each sample
  tototu <- apply(otu2a, MARGIN=1, FUN=sum)
  target.rare <- min(tototu)			# rarefy to this value

	### B) Rarefy
  oturare <- rrarefy(otu2a, min(tototu))
  
  ### C) Remove OTUs with 0 abundance
  otumax <- apply(oturare, 2, max)
  otu3 <- oturare[, otumax >= 1]
  
  ### D) Generate rarefaction curves (test)
  rarecurve(otu3, step=50,sample=min(tototu))
  
  ### E) Create key for OTU taxonomies and OTU
  otu.export <- as.matrix(t(otu3))
  OTU <- 1:nrow(otu.export)
  
  ### F) Finalize OTU table
  otu.final <- cbind(OTU,otu.export)
  
  ### G) Export OTU table
  write.csv(otu.final, paste("fileName",target.rare,".csv", sep=""))
   

### The following code can be run independently of above code, but requires the 
### file created in step G

#After re-arranging column names of taxonomies: import rarefied OTU 
# table. OTU IDs are now numbers, not taxonomies
otu.tmp <- read.csv("your_OTU_table.csv")

# Fill in empty taxonomy elements using fill.taxonomy.R()
otu4<- fill.taxonomy(otu)

# Transpose OTU table using t.otu()
otu5 <- t.otu(otu4, s=2, e=68)

# Create taxonomy matrix 
taxon.key <- otu4[,c(1,69:75)]

#OPTIONAL: remove singletons from OTU table of entire dataset
otumax <- apply(otu5, 2, max)   # or just otu2
#otumax
otu3 <- otu5[, otumax > 1]      
# create new taxon key
taxon.key1 <- taxon.key[otumax>1,]

#Transform rarefied OTU table
otu.pct <- decostand(otu3, method="total")
otu.hell <- decostand(otu3, method="hellinger")  # Hellinger transformation


      ######## nMDS  #####
mds <- metaMDS(otu.hell, k=3, distance="bray", trymax=20, autotransform=FALSE)

# Shepard plot to compare actual dissimilarities to computed dissimilarities
diss <- vegdist(otu.hell, distance="bray")
stressplot(mds, diss) ## check fit

# Obtain species and site scores
mds.sc <- scores(mds, display="sites", choices=c(1:3))
mds.spec <- scores(mds, display="species", choices=c(1:3))

## Plot results
  ###All samples#
plot(mds, type="none", choices=c(1,3), cex.lab=1.5,cex.axis=1.5,font=1.5,font.axis=2,font.lab=1.5)

points(mds.sc[,1], mds.sc[,2], pch=19, col="black")
points(mds.sc[,1], mds.sc[,2], pch=1,col="black",cex=1)
text(mds.sc[,1], mds.sc[,2],cex=.7, col="black")
mtext("nMDS of microbial communities, min x seqs/samp, rarefied", side=3, line=1, cex=1)
ordihull(mds, InsertYourFactorHere, lty=2)

legend("bottomleft", bty="n", legend="YourLegend", cex=0.6)


## PLOT WITH FITTED ENVTAL VARS  ##

### OPTIONAL: remove covarying environmental variables prior to fitting
big.cor(env, 0.6)
env.sub <- env[,-c(YourVarsToRemove)]
fit <- envfit(mds.sc[,1:2], env.sub,permutations=1000,
                   strata=yourOptionalStrata)
fit


library(grid)
fit.df <- data.frame(fit$vectors$arrows*sqrt(fit$vectors$r))

fit.df$species <- rnms

mds.sc <- data.frame(mds.sc)


ggplot(mds.sc, aes(NMDS1, NMDS2)) + 
    geom_point(data=mds.sc,aes(NMDS1, NMDS2, size=YourSizeParam, shape=yourShapeParam), colour="darkgrey") +          
    geom_segment(data=fit.df, aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow=arrow(length=unit(0.2,"cm")), colour="black", inherit_aes=FALSE) +
    geom_text(data=fit.df, aes(x=NMDS1*1.1, y=NMDS2*1.1, label=species), size=4.2)+ coord_fixed() +
    guides(shape=guide_legend("YourVariable"), 
           size=guide_legend(expression(paste("YourLegend", sep=""))) ) + 
    theme_bw() +
    theme(axis.title=element_text(size=14), 
          axis.text=element_text(size=11),
          legend.key=element_blank(),
          legend.key=element_rect(colour="white"), 
          legend.text=element_text(size=11), 
          legend.title=element_text(size=11))
    #+ ggtitle("YourTitle (3D stress=xxx)")


######## Plot on ordination results OTUs with > 5% abundance ######
spec1max <- apply(otu3,2,max)
mdsspec5 <- mds$species[spec1max>968,]   ##nMDS species coordinates for OTUs

spec5IDs <- taxon.key1[spec1max>968,]

  spec5IDs$Taxonomy <- factor(spec5IDs$Taxonomy)   ## correct levels for factors
  spec5IDs$Domain <- factor(spec5IDs$Domain)
  spec5IDs$Phylum <- factor(spec5IDs$Phylum)
  spec5IDs$Class <- factor(spec5IDs$Class)
  spec5IDs$Order <- factor(spec5IDs$Order)
  spec5IDs$Family <- factor(spec5IDs$Family)
  spec5IDs$Genus <- factor(spec5IDs$Genus)

### OPTIONAL: Replace proteo generic with proteo subphylum ###
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
  

##### OPTIONAL: CODE FOR REMOVING GRID LINES #####
    #theme(panel.grid.major=element_blank(), 
    #  panel.grid.minor=element_blank(),
    #  legend.key=element_rect(colour="white")) + 	  
    
    
######## Change plotting colors for taxa labels #####
colors <- c("black","red","green3","blue","magenta","darkgrey","cyan")
palette(colors)

##save to pdf file##
pdf(file="YourFile.pdf")
plot(mds.sc[,1], mds.sc[,2], pch=20,col="gray", cex=0.8, xlim=c(-1.5,1.5),
     ylim=c(-1,1.3))
text(mdsspec5[,1], mdsspec5[,2], labels=spec5IDs$Class, cex=0.9, 
     col=as.integer(spec5IDs$Phylum))
mtext("Bacterial taxa >= 10% abundance classified to family level")
legend("topleft", bty="n", legend=c(levels(spec5IDs$Phylum)), pt.bg=c(1:7),
       pch=22)
dev.off()


###### ---- CREATE BARPLOT OF DOMINANT TAXA --- ########
domtax <- otu.pct[,spec1max>YourMax]*100
domtax.t <- t(domtax)
taxname <- spec5IDs$Genus
labs <- alldat$NewID

mp <- barplot(domtax.t,col=rainbow(18),axisnames=FALSE, axes=FALSE,	 
	border=NA, axis.lty=1, ylim=c(0,100), ylab="% abundance")
text(mp,par("usr")[3], labels=labs, srt=65, adj=c(1.1,1.1), xpd=TRUE,     
	cex=0.4)
axis(side=2, lwd=0.5, tck= -0.02)
legend(-4,-8, xpd=TRUE, inset=0.02, bty="n", legend=taxname, 	  pch=22,pt.bg=rainbow(18), cex=0.7, ncol=4)


#########  -----------MULTIVARIATE ANALYSES ------------- #########

diss <- vegdist(otu.hell, method="YourMethod")
x <- betadisper(diss,group=YourGroup)
permutest(x)  
anova(x)
boxplot(x)
TukeyHSD(x)
#Note that removing groups that are represented by a single sample will greatly improve stability using either BC or Euc dist

##### PERMANOVA analysis for group membership    #####
adonis(diss~YourVar, data=env)

#########  MANTEL TESTS  #######
envdis <- vegdist(YourVar,method="YourMethod")
mantel(envdis,diss,strata=YourOptionalStrata)

####  ANOSIM for group membership   ####
mod <- anosim(otu.hell, grouping=YourGrouping,distance="YourDistanceMetric",
              strata=YourOptionalStrata)
