library(vegan)
library(tidyverse)
set.seed(496571)

phy.aob.pruned <- readRDS("Phyloseq_AOB_Pruned.rds")

#Transforma data with decostand and make a bray-curtis distance matrix
phy.aob.hel <- decostand(phyloseq::otu_table(phy.aob.pruned), method = "hellinger")
phy.aob.bc <- vegdist(phy.aob.hel, method = "bray")

#metaMDS for community data with k=3 (scree plot not shown showed best fit at 3).
phy.aob.nmds <- metaMDS(phy.aob.bc, autotransform = FALSE, distance = "bray", maxtry = 100, k=3)

#Mantel test for continuous variables
sample.data.aob <- data.frame(sample_data(phy.aob.pruned))
dist.soilph.aob <- vegdist(sample.data.aob$SOIL.PH, method = "euclidean")
dist.nmin.aob <- vegdist(sample.data.aob$N.MIN, method = "euclidean")
dist.ammpost.aob <- vegdist(sample.data.aob$AMM.CORRECTED.POST, method = "euclidean")
mantel.soilph.aob <- mantel(phy.aob.bc, dist.soilph.aob, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel.nmin.aob <- mantel(phy.aob.bc, dist.nmin.aob, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel.ammpost.aob <- mantel(phy.aob.bc, dist.ammpost.aob, method = "spearman", permutations = 9999, na.rm = TRUE)

#adonis (permanvoa) works on discrete variables, so STAND or create groupings based on environmental variables (A vs B or A,B,C) to test.
adonis.aob <- adonis2(phy.aob.bc ~STAND , data =sample.data.aob)

myresults.mantel.adonis.aob <- list(mantel.soilph.aob, mantel.nmin.aob, mantel.ammpost.aob, adonis.aob)
capture.output(myresults.mantel.adonis.aob, file = "Mantel_Adonis_AOB_Results.txt")

#Create NMDS Plot
nmds.data.aob <- as.data.frame(scores(phy.aob.nmds))
nmds.data.aob$sample <- rownames(nmds.data.aob)
nmds.data.aob$STAND <- as.factor(as.integer(gsub("Stand_", "", sample.data.aob$STAND)))
                                 
#Coloring can be based on environmental factor but removing scale_fill_manual option just defaults
nmds.colors.aob <- c("darkblue", "coral3","coral", "dodgerblue3", "navy", "pink2", "pink", "dodgerblue")
nmds.plot.aob <- ggplot(data = nmds.data.aob, aes(x = NMDS1, y = NMDS2, fill = STAND)) + geom_point(aes(), colour = "black", pch = 21, size = 3) + xlab(paste0("NMDS1")) + ylab(paste0("NMDS2")) + stat_ellipse(aes(x = NMDS1, y= NMDS2, color = STAND), level = 0.5) + scale_fill_manual(values = nmds.colors.aob) + coord_fixed() + theme_linedraw(base_size = 18) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(face = "bold")) + theme(text=element_text(size = 20)) + ggtitle("NMDS of ASV amoA of AOB")
ggsave("NMDS_AOB.pdf", width = 8.5, height = 10)
