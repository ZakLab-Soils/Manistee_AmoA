library(vegan)
library(tidyverse)

phy.aoa.pruned <- readRDS("Phyloseq_AOA_Pruned.rds")

#Transforma data with decostand and make a bray-curtis distance matrix
phy.aoa.hel <- decostand(phyloseq::otu_table(phy.aoa.pruned), method = "hellinger")
phy.aoa.bc <- vegdist(phy.aoa.hel, method = "bray")

#metaMDS for community data with k=3 (scree plot not shown showed best fit at 3).
phy.aoa.nmds <- metaMDS(phy.aoa.bc, autotransform = FALSE, distance = "bray", maxtry = 100, k=3)

#Mantel test for continuous variables
sample.data.aoa <- data.frame(sample_data(phy.aoa.pruned))
dist.soilph.aoa <- vegdist(sample.data.aoa$SOIL.PH, method = "euclidean")
dist.nmin.aoa <- vegdist(sample.data.aoa$N.MIN, method = "euclidean")
dist.ammpost.aoa <- vegdist(sample.data.aoa$AMM.CORRECTED.POST, method = "euclidean")
mantel.soilph.aoa <- mantel(phy.aoa.bc, dist.soilph.aoa, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel.nmin.aoa <- mantel(phy.aoa.bc, dist.nmin.aoa, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel.ammpost.aoa <- mantel(phy.aoa.bc, dist.ammpost.aoa, method = "spearman", permutations = 9999, na.rm = TRUE)

#adonis (permanvoa) works on discrete variables, so STAND or create groupings based on environmental variables (A vs B or A,B,C) to test.
adonis2(phy.aoa.bc ~STAND , data =sample.data.aoa)

#Create NMDS Plot
nmds.data.aoa <- as.data.frame(scores(phy.aoa.nmds))
nmds.data.aoa$sample <- rownames(nmds.data.aoa)
nmds.data.aoa$STAND <- as.factor(as.integer(gsub("Stand_", "", sample.data.aoa$STAND))
                                 
#Coloring can be based on environmental factor but removing scale_fill_manual option just defaults
nmds.colors <- c("darkblue", "coral3", "dodgerblue4", "dodgerblue2", "pink2", "coral4", "purple")
nmds.plot.aoa <- ggplot(data = nmds.data.aoa, aes(x = NMDS1, y = NMDS2, fill = STAND)) + geom_point(aes(), colour = "black", pch = 21, size = 3) + xlab(paste0("NMDS1")) + ylab(paste0("NMDS2")) + stat_ellipse(aes(x = NMDS1, y= NMDS2, color = STAND), level = 0.5) + scale_fill_manual(values = nmds.colors) + coord_fixed() + theme_linedraw(base_size = 18) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(face = "bold")) + theme(text=element_text(size = 20)) + ggtitle("NMDS of ASV amoA of AOA")
