setwd("/Users/stanish/Documents/R")

library(vegan)
library(ggplot2)
source("Fill.Taxonomy.r")
library(BiodiversityR)

spec <- read.csv("OhioRiver_OTU_table_noCPs_19350.csv")
env <- read.csv("OhioRiver_Metadata.csv", row.names=1)
chem <- read.csv("OhioRiver_Chemistry.csv", row.names=1)

spec1 <- fill.taxonomy(spec)
taxon.key <- spec1[,(ncol(spec1)-6):ncol(spec1)]
spec.t <- t.otu(spec, s=2,e=68)

spec2 <- spec.t
colnames(spec2) <- taxon.key$Genus

xr <- rankabundance(spec2)
rankabunplot(xr,scale="proportion",type="h",lwd=3)

