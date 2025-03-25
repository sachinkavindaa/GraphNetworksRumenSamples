library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(readxl)
library(ggplot2)
library(gridExtra)

data <- read_excel("ORFMatrix_Heritable.xlsx")
data <- as.data.frame(data)
rownames(data) <- data[, 1]  # Replace "1" with the correct column index
data <- data[, -1] 

gsg <- goodSamplesGenes(data, verbose =3)
summary(gsg)
gsg$allOK


data_samples <- t(data)

sampleTree <- hclust(dist(data_samples), method = "average")
plot(sampleTree)



