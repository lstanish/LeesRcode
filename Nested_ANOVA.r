###### Nested ANOVA ########

setwd("/Users/stanish/Documents/R")   # for Macssetwd("C:\\Users\\Lee Stanish\\Documents\\R")

library(vegan)
library(ggplot2)
source("Fill.Taxonomy.r")
source("inksplot3.r")

spec <- read.csv("OhioRiver_OTU_table_noCPs_19350.csv")
env <- read.csv("OhioRiver_Metadata.csv", row.names=1)
chem <- read.csv("OhioRiver_Chemistry.csv", row.names=1)


## Locate samples from Munis with no reps prior to analyses
beta <- vector(length=17)
for (i in 1:17) {
  beta[i] <- (length(which(env$Muni.Numeric==i)))
}
noreps <- which(beta<2)
reps <- which(beta>1)
beta1 <- vector(length=length(noreps))
for(i in 1:length(noreps)) {
  beta1[i] <- which(env$Muni.Numeric==noreps[i], arr.ind=TRUE)
}

# Remove samples from singleton munis
env.rm <- env[-beta1,]
spec.rm <- spec[,-(beta1+1)]


######Ohio River alone: convert OTU table to percents ####
spec1 <- fill.taxonomy(spec.rm)
taxon.key <- spec1[,(ncol(spec1)-6):ncol(spec1)]
spec.t <- t.otu(spec.rm, s=2,e=66)
#spec.pa <- decostand(spec.t, method="pa", margin=1)

spec.pct <- matrix(data=NA, nrow=nrow(spec.t),ncol=ncol(spec.t))
for(i in 1:ncol(spec.t)) {
  for(j in 1:nrow(spec.t)) {
    spec.pct[j,i] <- spec.t[j,i]/sum(spec.t[j,])
  }
}




############### RESUME ANALYSIS WITH CHOSEN PERCENT OTU TABLE #####

ntax <- ncol(spec.pct)

#### Loop for obtaining Pvalues from nested ANOVAS ###
muni <- env.rm$Municipality
temp.pvals <- matrix(data=NA, nrow=ntax, ncol=3)
for(i in 1:ntax) {
  mod <- aov(spec.pct[,i]~env.rm$ChlType+muni%in%env.rm$ChlType)
  mod.sum <- summary(mod)
  temp.pvals[i,] <- mod.sum[[1]]$'Pr(>F)'
}

#### Loop for obtaining P/A sums for each taxon by chlorine type ##
#pa.sums <- matrix(NA, nrow=ntax, ncol=2, dimnames=list(taxon.key$Genus, 
#                                             c("Chloramine", "Chlorine")))
#for(j in 1:ntax) {
#  pa.sums[j,] <- tapply(spec.pa[,j], env$ChlType, sum)
#}

#library(dplyr)
#pa.sums.df <- data.frame(pa.sums)
#pa.sums.df <- mutate(pa.sums.df, chloramine.norm = Chloramine/9, 
#             chlorine.norm = Chlorine/58, ratio = chlorine.norm/chloramine.norm)



## Option #1: Select taxon columns where Pval among group is at least 100x greater than
### within group Pval (eg Chltype way more sig than Muni)
sigdiff <- which(temp.pvals[,1]/temp.pvals[,2]<0.01)
sig.tax <- temp.pvals[sigdiff,]
spec.sig <- spec.pct[,sigdiff]
taxname <- factor(taxon.key$Genus[sigdiff])
colnames(spec.sig) <- taxon.key$Taxonomy[sigdiff]
rownames(spec.sig) <- rownames(env.rm)
#write.csv(spec.sig,"OR_sigtax_nestANOVA.csv")


## Option #2: Select taxon columns where Pval among group is at least 100x greater than
### within group Pval (eg Chltype way more sig than Muni) *INCLUDING* MLE1-12
### (326) and Mycobacterium (41)
sigdiff <- c(which(temp.pvals[,1]/temp.pvals[,2]<0.01),326,41)
sig.tax <- temp.pvals[sigdiff,]
spec.sig <- spec.pct.rm[,sigdiff]
taxname <- factor(taxon.key$Genus[sigdiff])
colnames(spec.sig) <- taxon.key$Taxonomy[sigdiff]
rownames(spec.sig) <- rownames(env.rm)

###### ---- CREATE BARPLOT OF DOMINANT TAXA --- ########
chl <- which(env.rm$ChlType=="Chlorine")
chm <- which(env.rm$ChlType=="Chloramine")

mx <- apply(spec.sig,2,max)
domtax <- spec.sig[,mx>0.01]*100
domnames <- factor(taxname[mx>0.01])
domtax.t <- t(domtax)

domtax.sort <- rbind(domtax[chl,],domtax[chm,])
domtax.names <- row.names(domtax.t)
plotIDs <- as.character(env.rm$newID)
dom.labs <- c(plotIDs[chl],plotIDs[chm])
rownames(domtax.sort) <- dom.labs
domtax.tax <- colnames(domtax.sort)
taxonomy <- gsub("Bacteria/|Proteobacteria/|Rhizobiales/","",domtax.tax) 
taxonomy1 <- gsub("Nitrospira/Nitrospirales/Nitrospiraceae/|Bradyrhizobiaceae/","",taxonomy)
taxonomy2 <- gsub("Methylobacteriaceae/|Rhodospirillales/Rhodospirillaceae/","",taxonomy1)
taxonomy3 <- gsub("Sphingomonadales/Sphingomonadaceae/|Burkholderiales/Comamonadaceae/","",taxonomy2)
taxonomy4 <- gsub("Hydrogenophilales/Hydrogenophilaceae/|Nitrosomonadaceae/","",taxonomy3)
taxonomy5 <- gsub("Nitrosomonadales/|Gallionellaceae/|Rhodocyclales/","",taxonomy4)
taxonomy6 <- gsub("Sphingobacteriia/Sphingobacteriales/Chitinophagaceae/","",taxonomy5)
taxonomy7 <- gsub("Actinobacteria/Corynebacteriales/Mycobacteriaceae/","",taxonomy6)

colnames(domtax.sort) <- taxonomy7

## Option: Normalize percentages by excluding MLE1-12 and Mycos
domtax.sort1 <- data.frame(domtax.sort)
domtax.sort1 <- mutate(domtax.sort1,normsum=100-(Cyanobacteria.MLE1.12+ Actinobacteria.Mycobacterium) )
domtax.sort2 <- domtax.sort1[,1:(ncol(domtax.sort1)-3)]
domtax.pct <- matrix(data=NA, nrow=nrow(domtax.sort2),ncol=ncol(domtax.sort2))
for(i in 1:ncol(domtax.sort2)) {
  for(j in 1:nrow(domtax.sort2)) {
    domtax.pct[j,i] <- (domtax.sort2[j,i]/domtax.sort1$normsum[j])*100
  }
}
colnames(domtax.pct) <- taxonomy7
rownames(domtax.pct) <- dom.labs


### Nov 26 2015: fix a.tax and b.tax so they work for non-transposed table###
a.tax <- c(1:2,4:7,9:13,19)
b.tax <- c(8,14:17,3,18)
domtax.a <- domtax.sort[a.tax,]
domtax.b <- domtax.sort[b.tax,]


### Inksplot ###
# Non-normalized percentages #
pdf("OR_inksplot_PLoS_revised.pdf", width=8.5, height=5, bg="white")
inksplot(domtax.sort)
dev.off()


# Normalized with respect to MLE1-12 and Mycobacteria
pdf("OR_inksplot_Fig_norm.pdf", width=8.5, height=5, bg="white")
inksplot(domtax.pct)
dev.off()



### "A" PLOT ####
library(RColorBrewer)
colors <- brewer.pal(12,"Set3")
pdf("OR_barplot_Fig5A_colortest.pdf", width=6, height=4.3, bg="white")
mp.a <- barplot(height=domtax.t[a.tax,],col=colors,axisnames=FALSE,    
                axes=FALSE, border=NA, axis.lty=0, ylim=c(0,90), ylab="% abundance")
text(mp.a,par("usr")[3], labels=dom.labs, srt=65, adj=c(1.1,1.1), xpd=TRUE,     
     cex=0.4)
axis(side=2, lwd=0.5, tck= -0.02)
legend(-5,-9, xpd=TRUE, inset=0.02, bty="n", legend=domnames[a.tax], 
       pch=22,col=colors, pt.bg=colors, cex=0.7, ncol=4)
dev.off()


### "B" PLOT ###
library(gplots)
colors <- brewer.pal(7,"Set2")

pdf("OR_barplot_Fig5B_colortest.pdf", width=6, height=4.3, bg="white")
mp.b <- barplot(domtax.t[b.tax,],col=colors,axisnames=FALSE, axes=FALSE,   
                border=NA, axis.lty=0, ylim=c(0,12), ylab="% abundance")
text(mp.b,par("usr")[3], labels=dom.labs, srt=65, adj=c(1.1,1.1), xpd=TRUE,     
     cex=0.4)
axis(side=2, lwd=0.5, tck= -0.02)
legend(-2,-1.2, xpd=TRUE, inset=0.02, bty="n", legend=domnames[b.tax], 
       pch=22,col=colors, pt.bg=colors, cex=0.7, ncol=4)
dev.off()
