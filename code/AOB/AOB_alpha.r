#Alpha diversity
library(phyloseq)
library(vegan)
library(ggplot2)

phy.aob.pruned <- readRDS("Phyloseq_AOB_Pruned.rds")

phy.aob.rare <- rarefy_even_depth(phy.aob.pruned, rngseed =496571, sample.size=min(sample_sums(phy.aob.pruned)), replace=F)
sample_data(phy.aob.rare)$STAND.CLEAN <- as.factor(gsub("Stand_*", "", fixed=FALSE, sample_data(phy.aob.rare)$STAND))

plot.shannon.data.aob <- plot_richness(phy.aob.rare, x="STAND.CLEAN", measures = "Shannon", color = "STAND.CLEAN")$data
plot.shannon.data.aob$STAND.CLEAN <- factor(plot.shannon.data.aob$STAND.CLEAN, levels = c("58", "7", "9", "41", "100", "24", "22", "6"))

plot.chao.data.aob <- plot_richness(phy.aob.rare, x="STAND.CLEAN", measures = "Chao1", color = "STAND.CLEAN")$data
plot.chao.data.aob$STAND.CLEAN <- factor(plot.chao.data.aob$STAND.CLEAN, levels = c("58", "7", "9", "41", "100", "24", "22", "6"))

boxplot.colors.aob <- c("coral4", "coral3", "coral", "pink2", "dodgerblue", "dodgerblue2", "dodgerblue4", "darkblue")

shannon.plot.aob <- ggplot(data=plot.shannon.data.aob, aes(x=STAND.CLEAN, y=value, fill = STAND.CLEAN)) + geom_boxplot(alpha=0.7) + xlab(paste0("STAND")) + ylab(paste0(" Diversity")) + guides(fill=guide_legend(title="Stand"))+scale_fill_manual(values = boxplot.colors.aob)+theme(text=element_text(size = 18)) + theme(panel.background = element_blank()) + ggtitle("Shannon Diversity of amoA in AOB")
ggsave("Shannon_boxplot_AOB.pdf", width = 8.5, height = 10)

chao.plot.aob <- ggplot(data=plot.chao.data.aob, aes(x=STAND.CLEAN, y=value, fill = STAND.CLEAN)) + geom_boxplot(alpha=0.7) + xlab(paste0("STAND")) + ylab(paste0("Chao1 Estimation")) + guides(fill=guide_legend(title="Stand"))+scale_fill_manual(values = boxplot.colors.aob)+theme(text=element_text(size = 18)) + theme(panel.background = element_blank()) + ggtitle("Chao1 of amoA in AOB")
ggsave("Chao_boxplot_AOB.pdf", width = 8.5, height = 10)

##Rarefaction by Stand
sample_data(phy.aob.rare)$STAND.CLEAN <- as.factor(gsub("Stand_*", "", fixed=FALSE, sample_data(phy.aob.rare)$STAND))
sample_data(phy.aob.rare)$STAND.CLEAN <- factor(sample_data(phy.aob.rare)$STAND.CLEAN, levels = c("41", "7", "9", "58", "100", "6", "24", "22")) 
phy.aob.rare.merged <- merge_samples(phy.aob.rare, "STAND.CLEAN")
#To average the rarefied data that is merged by the number of plots in each stand
stand.counts.aob <- c(4, 6, 4, 4, 6, 5, 6, 4)

rarefaction.aob <- rarecurve(as(round(otu_table(phy.aob.rare.merged)/stand.counts.aob), "matrix"), step = 100, col = boxplot.colors.aob)
names(rarefaction.aob) <- sample_names(phy.aob.rare.merged)

protox <- mapply(FUN = function(x, y) {
    mydf <- as.data.frame(x)
    colnames(mydf) <- "value"
    mydf$SampleID <- y
    mydf$subsample <- attr(x, "Subsample")
    mydf
}, x = rarefaction.aob, y = as.list(names(rarefaction.aob)), SIMPLIFY = FALSE)
xy <- do.call(rbind, protox)
rownames(xy) <- NULL
xy$SampleID <- factor(xy$SampleID, levels=c("41", "7", "9", "58", "100", "6", "24","22"))

Rarefaction.aob.ggplot <- ggplot(xy, aes(x=subsample, y=value, group=SampleID)) + geom_line(aes(color=SampleID), size = 0.8) + scale_color_manual(values=boxplot.colors.aob) + theme_bw() +ggtitle("Rarefaction Curves by Stand") + theme(plot.title =element_text(hjust = 0.5)) + xlab("Subsampled Reads") + ylab("Observed ASVs") + theme(plot.title = element_text(size = 20),legend.key.size = unit(.8,"cm"),legend.text = element_text(size = 16), legend.title=element_text(size = 18), axis.text = element_text(size = 10), axis.title = element_text(size = 16)) + guides(color=guide_legend(title = "Stand"))
ggsave("Rarefaction_AOB_Stand_Averages.pdf", width=8.5, height = 10)
