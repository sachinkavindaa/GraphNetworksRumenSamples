library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(readxl)
library(ggplot2)
library(gridExtra)

data <- read_excel("ORF_matrix_all_heritable-2.xlsx")
rownames(data) <- data$AnimalID
data$AnimalID <- NULL  # Remove the AnimalID column

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

data_matrix <- as.matrix(data)

norm.counts <- t((data_matrix)) 

power <- c(c(1:10), seq(from = 12, to =50 , by = 2))

sft <- pickSoftThreshold(norm.counts,
                  powerVector = power,
                  networkType = "signed",
                  verbose = 5)

sft.data <- sft$fitIndices

a1 <- ggplot(sft.data, aes(Power,SFT.R.sq,label = power))+
  geom_point() +
  geom_text(nudge_y = 0.1)+
  geom_hline(yintercept = 0.8, color = 'red' ) +
  labs(x = 'Power', y = 'Scale free topolog model fit, sigend R^2') + 
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k. , label = Power )) +
  geom_point() +
  geom_text(nudge_y = 0.1 ) +
  labs(x = 'Power' ,  y= 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

# Find the smallest power where R² ≥ 0.8
optimal_power <- sft.data$Power[
  which(sft.data$SFT.R.sq >= 0.8)[1]  # First power meeting the threshold
]

# If no power reaches R² ≥ 0.8, use the power closest to 0.8
if (is.na(optimal_power)) {
  optimal_power <- sft.data$Power[which.max(sft.data$SFT.R.sq)]
}

print(paste("Optimal power:", optimal_power))

norm.counts[] <- sapply(norm.counts,as.numeric)

soft_power <- 12
temp_cor <- cor 
cor <- WGCNA::cor

# Check dimensions (genes should be rows, samples as columns)
dim(norm.counts)  # Expected: [genes] x [samples], e.g., 717 x 100

# Check class (must be a numeric matrix)
class(norm.counts)  # Should return "matrix" or "data.frame"

bwnet <- blockwiseModules(norm.counts,
                 maxBlockSize = 20000,
                 TOMType = "signed",
                 power = soft_power,
                 mergeCutHeight = 0.25,
                 numericLabels = FALSE,
                 randomSeed = 1234,
                 verbose = 3)

cor <- temp_cor

module_eigengenes <- bwnet
head(module_eigengenes)

table(bwnet$colors)


plotDendroAndColors(bwnet$dendrograms[[1]],cbind(bwnet$unmergedColors,bwnet$colors),
                          c("unmerged","merged"),
                          dendrapply = FALSE,
                          addGuide = TRUE,
                          hang = 0.03,
                          guidehang = 0.05)


# Original data orientation (samples x genes)
dim(data)  # Should be [41 samples] x [717 genes]

# Transpose to genes x samples
norm.counts <- t(as.matrix(data))
dim(norm.counts)  # Should be [717 genes] x [41 samples]

# Check for zero-variance genes
gene_vars <- apply(norm.counts, 1, var, na.rm = TRUE)
sum(gene_vars == 0)  # Number of genes to remove

# Keep only genes with variance > 0
norm.counts_clean <- norm.counts[gene_vars > 0, ]
dim(norm.counts_clean)  # e.g., 700 genes x 41 samples

bwnet <- blockwiseModules(
  norm.counts_clean,  # Genes = rows, samples = columns
  power = 12,
  maxBlockSize = 20000,  # Process all genes in one block
  networkType = "signed",
  TOMType = "signed",
  numericLabels = FALSE,
  verbose = 3
  
  plotDendroAndColors(
    dendro = bwnet$dendrograms[[1]], 
    colors = cbind(bwnet$unmergedColors, bwnet$colors),
    groupLabels = c("Unmerged", "Merged"),
    dendroLabels = FALSE,
    hang = 0.03,
    main = "Gene Clustering and Module Colors"
  )


  # Check number of modules before/after merging
  table(bwnet$unmergedColors)  # Unmerged modules
  table(bwnet$colors)          # Merged modules

  plotDendroAndColors(
    dendro = bwnet$dendrograms[[1]], 
    colors = cbind(bwnet$unmergedColors, bwnet$colors),
    groupLabels = c("Unmerged", "Merged"),  # Must match number of color columns
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = "Gene Clustering and Module Assignment"
  )

  

