########################
#    WGCNA analysis    #
########################

# Some part of the code were taken from tutorial available online.

library(WGCNA)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(fastDummies)
library(gridExtra)
library(mlr3measures)

conditions.ad.l4 <- condition_table[condition_table$stage=="adult" | condition_table$stage=="larva f",]
keep <- row.names(ad.l4.tpm[rowSums(ad.l4.tpm>1)>=2,]) 
datacleaned <- t.ad.l4.ctf[,colnames(t.ad.l4.ctf) %in% keep]
asin.t.ad.l4.ctf <- asinh(datacleaned)                  

#### 1. Choose the soft power ----

### compute R^2 for all exponents
power = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(asin.t.ad.l4.ctf, powerVector = power, networkType = "signed", verbose = 5) #if not specified, the network type used is unsigned

### plot the results to determine the soft threshold that better fits a free scale topology

sft.data <- sft$fitIndices

# visualization to pick power
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.08) + 
  geom_hline(yintercept = 0.80, color = 'red') +
  labs(x=NULL, y=NULL) +
  #labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 500) +
  #geom_hline(yintercept = 100, color = 'red') +
  labs(x='Power', y=NULL) +
  #labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

#tiff("choice_power_signed_CTF.tiff", res=1000, width=8000, height=6424)
grid.arrange(a1, a2, nrow = 2)
#dev.off()

#### 2. Generate the network ----

soft_power <- 9

# performes automatic network construction
bwnet <- blockwiseModules(asin.t.ad.l4.ctf,
                          TOMType = "signed",
                          power = soft_power,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors underneath
#tiff("dendrogram_CTF.tiff", res=300, width=3000, height=1980)
plotDendroAndColors(
  bwnet$dendrograms[[1]],
  bwnet$colors[bwnet$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05)
#dev.off()

# Module Eigengene
module_eigengenes <- bwnet$MEs
head(module_eigengenes)

## Save data 
module_df <- data.frame(
  gene_id = names(bwnet$colors),
  colors = bwnet$colors
)

#### 3. Module-trait (tissue) association ----

# Define numbers of genes and samples
nGenes = ncol(asin.t.ad.l4.ctf)
nSamples = nrow(asin.t.ad.l4.ctf)
bin.tissues <- fastDummies::dummy_cols(conditions.ad.l4$tissue1, remove_selected_columns=TRUE) #conditions or tissue1
rownames(bin.tissues) <- rownames(conditions.ad.l4)
names(bin.tissues) <- sub('.data_', '', names(bin.tissues))
module.trait.correlation = cor(MEs0, bin.tissues, method = "pearson") #one continuous and one categorical variable
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples) #calculate the p-value associated with the correlation
  
Signif.star <- symnum(module.trait.Pvalue, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                 symbols = c("***", "**", "*", ".", " "))

textMatrixstar = paste(signif(module.trait.correlation, 2), " ", Signif.star, sep = "");
dim(textMatrixstar) = dim(module.trait.correlation)

#tiff("Heatmap_corr_modules_child_tissues_CTF.tiff", res=300, width=3300, height=3100) 
par(mar = c(6, 10, 2.5, 0.5)) #padding bottom, left, up, right                        
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = module.trait.correlation,
               xLabels = names(bin.tissues),
               yLabels = names(MEs0),
               ySymbols = names(MEs0),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrixstar,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#dev.off()
