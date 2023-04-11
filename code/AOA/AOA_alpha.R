##Alpha diversity
library(phyloseq)
library(vegan)
library(ggplot2)

phy.aoa.pruned <- readRDS("Phyloseq_AOA_Pruned.rds")

phy.aoa.rare <- rarefy_even_depth(phy.aoa.pruned, rngseed =49657120, sample.size=min(sample_sums(phy.aoa.pruned)), replace=F)
sample_data(phy.aoa.rare)$STAND.CLEAN <- as.factor(gsub("Stand_*", "", fixed=FALSE, sample_data(phy.aoa.rare)$STAND))

plot.shannon.data.aoa <- plot_richness(phy.aoa.rare, x="STAND.CLEAN", measures = "Shannon", color = "STAND.CLEAN")$data
plot.shannon.data.aoa$STAND.CLEAN <- factor(plot.shannon.data$STAND.CLEAN, levels = c("58", "7", "41", "100", "24", "22", "6"))

plot.chao.data.aoa <- plot_richness(phy.aoa.rare, x="STAND.CLEAN", measures = "Chao1", color = "STAND.CLEAN")$data
plot.chao.data.aoa$STAND.CLEAN <- factor(plot.chao.data$STAND.CLEAN, levels = c("58", "7", "41", "100", "24", "22", "6"))

boxplot.colors.aoa <- c("coral4", "coral3", "pink2", "purple", "dodgerblue2", "dodgerblue4", "darkblue")

shannon.plot.aoa <- ggplot(data=plot.shannon.data.aoa, aes(x=STAND.CLEAN, y=value, fill = STAND.CLEAN)) + geom_boxplot(alpha=0.7) + xlab(paste0("STAND")) + ylab(paste0(" Diversity")) + guides(fill=guide_legend(title="Stand"))+scale_fill_manual(values = boxplot.colors.aoa)+theme(text=element_text(size = 18)) + theme(panel.background = element_blank()) + ggtitle("Shannon Diversity of amoA in AOA")
ggsave("Shannon_boxplot_AOA.pdf", width = 8.5, height = 10)

chao.plot.aoa <- ggplot(data=plot.chao.data.aoa, aes(x=STAND.CLEAN, y=value, fill = STAND.CLEAN)) + geom_boxplot(alpha=0.7) + xlab(paste0("STAND")) + ylab(paste0("Chao1 Estimation")) + guides(fill=guide_legend(title="Stand"))+scale_fill_manual(values = boxplot.colors.aoa)+theme(text=element_text(size = 18)) + theme(panel.background = element_blank()) + ggtitle("Chao1 of amoA in AOA")
ggsave("Chao_boxplot_AOA.pdf", width = 8.5, height = 10)

##Rarefaction by Stand
sample_data(phy.aoa.rare)$STAND.CLEAN <- as.factor(gsub("Stand_*", "", fixed=FALSE, sample_data(phy.aoa.rare)$STAND))
sample_data(phy.aoa.rare)$STAND.CLEAN <- factor(sample_data(phy.aoa.rare)$STAND.CLEAN, c("41", "7", "58", "100", "6", "24", "22")) 
phy.aoa.rare.merged <- merge_samples(phy.aoa.rare, "STAND.CLEAN")
#To average the rarefied data that is merged by the number of plots in each stand
stand.counts.aoa <- c(5, 5, 5, 5, 6, 6, 6)

rarefaction.aoa <- rarecurve(as(round(otu_table(phy.aoa.rare.merged)/stand.counts.aoa), "matrix"), step = 100, col = boxplot.colors
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
xy$SampleID <- factor(xy$SampleID, levels=c("41", "7", "58", "100", "6", "24","22"))

Rarefaction.aoa.ggplot <- ggplot(xy, aes(x=subsample, y=value, group=SampleID)) + geom_line(aes(color=SampleID), size = 0.8) + scale_color_manual(values=boxplot.colors) + theme_bw() +ggtitle("Rarefaction Curves by Stand") + theme(plot.title =element_text(hjust = 0.5)) + xlab("Subsampled Reads") + ylab("Observed ASVs") + theme(plot.title = element_text(size = 20),legend.key.size = unit(.8,"cm"),legend.text = element_text(size = 16), legend.title=element_text(size = 18), axis.text = element_text(size = 10), axis.title = element_text(size = 16)) + guides(color=guide_legend(title = "Stand"))
ggsave("Rarefaction_AOA_Stand_Averages.pdf")
