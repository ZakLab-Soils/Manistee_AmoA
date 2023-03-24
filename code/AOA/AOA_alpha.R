##Alpha diversity
library(phyloseq)
library(vegan)
library(ggplot2)

phy.aoa.pruned <- readRDS("Phyloseq_AOA_Pruned.rds")

phy.aoa.rare <- rarefy_even_depth(phy.aoa.pruned, rngseed =49657120, sample.size=min(sample_sums(phy.aoa.pruned)), replace=F)
sample_data(phy.aoa.rare)$STAND.CLEAN <- as.factor(gsub("Stand_*", "", fixed=FALSE, sample_data(phy.aoa.rare)$STAND))

plot.shannon.data <- plot_richness(phy.aoa.rare, x="STAND.CLEAN", measures = "Shannon", color = "STAND.CLEAN")$data
plot.shannon.data$STAND.CLEAN <- factor(plot.chao.data$STAND.CLEAN, levels = c("58", "6", "7", "41", "24", "100", "22"))

plot.chao.data <- plot_richness(phy.aoa.rare, x="STAND.CLEAN", measures = "Chao1", color = "STAND.CLEAN")$data
plot.chao.data$STAND.CLEAN <- factor(plot.chao.data$STAND.CLEAN, levels = c("58", "6", "7", "41", "24", "100", "22"))

boxplot.colors <- c("coral4", "coral3", "coral2", "coral", "dodgerblue", "dodgerblue3", "dodgerblue4" )

shannon.plot <- ggplot(data=plot.shannon.data, aes(x=STAND.CLEAN, y=value, fill = STAND.CLEAN)) + geom_boxplot(alpha=0.7) + xlab(paste0("STAND")) + ylab(paste0(" Diversity")) + guides(fill=guide_legend(title="Stand"))+scale_fill_manual(values = boxplot.colors)+theme(text=element_text(size = 18)) + ggtitle("Shannon Diversity of amoA in AOA")

chao.plot <- ggplot(data=plot.chao.data, aes(x=STAND.CLEAN, y=value, fill = STAND.CLEAN)) + geom_boxplot(alpha=0.7) + xlab(paste0("STAND")) + ylab(paste0("Chao1 Estimation")) + guides(fill=guide_legend(title="Stand"))+scale_fill_manual(values = boxplot.colors)+theme(text=element_text(size = 18)) + ggtitle("Chao1 of amoA in AOA)

##Rarefaction by Stand

sample_data(phy.aoa.rare)$STAND.CLEAN <- as.factor(gsub("Stand_*", "", fixed=FALSE, sample_data(phy.aoa.rare)$STAND))
phy.aoa.rare.merged <- merge_samples(phy.aoa.rare, "STAND.CLEAN")

rarefaction.aoa <- rarecurve(as(otu_table(phy.aoa.rare.merged), "matrix"), step = 100, col=boxplot.colors)
names(rarefaction.aoa) <- sample_names(phy.aoa.rare.merged)

protox <- mapply(FUN = function(x, y) {
    mydf <- as.data.frame(x)
    colnames(mydf) <- "value"
    mydf$SampleID <- y
    mydf$subsample <- attr(x, "Subsample")
    mydf
}, x = rarefaction.aoa, y = as.list(names(rarefaction.aoa)), SIMPLIFY = FALSE)
xy <- do.call(rbind, protox)
rownames(xy) <- NULL

Rarefaction.aoa.ggplot <- ggplot(xy, aes(x=subsample, y=value, group=SampleID)) + geom_line(aes(color=SampleID), size = 0.8) + scale_color_manual(values=boxplot.colors) + theme_bw() +ggtitle("Rarefaction Curves by Stand") + theme(plot.title =element_text(hjust = 0.5)) + xlab("Subsampled Reads") + ylab("Observed ASVs") + theme(plot.title = element_text(size = 20),legend.key.size = unit(.8,"cm"),legend.text = element_text(size = 16), legend.title=element_text(size = 18), axis.text = element_text(size = 10), axis.title = element_text(size = 16)) + guides(color=guide_legend(title = "Stand"))
