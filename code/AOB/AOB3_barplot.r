#Create and analyze barplot for Cluster groupings
library(tidyverse)
library(ggplot2)

#phy.aob.pruned <- readRDS("Phyloseq_AOB_Pruned.rds")
Cl.all.aob.tbl <- readRDS("CL_ALL_AOB_Table.rds")

CL1.aob.tbl <-  Cl.all.aob.tbl %>% group_by(CL1, STAND, PLOT) %>% summarize(COUNT = sum(COUNT))%>% ungroup()%>% select(CL1, STAND, PLOT, COUNT)
CL1.aob.tmp.tbl <- CL1.aob.tbl  %>% group_by(CL1) %>% mutate(CL1.TOTAL = sum(COUNT))%>% ungroup() %>% filter(CL1.TOTAL > 0) %>% group_by(PLOT) %>% mutate(PLOT.COUNT = sum(COUNT))%>% ungroup() %>% group_by(CL1, PLOT) %>% mutate(PROP = sum(COUNT)/PLOT.COUNT) %>% ungroup() %>% mutate(HELLINGER = sqrt(PROP))

#At the 0.1 pairwise clustering, all AOB ASVs collapse into 6 clusters
CL1.stand.aob.mean <- aggregate(PROP~CL1+STAND, CL1.aob.tmp.tbl, mean)
CL1.stand.aob.mean$STAND <- gsub("Stand_*", "", fixed=FALSE, CL1.stand.aob.mean$STAND)
CL1.stand.aob.mean <- CL1.stand.aob.mean[order(CL1.stand.aob.mean$STAND), ]
CL1.stand.aob.mean$STAND <- factor(CL1.stand.aob.mean$STAND, levels = c("58", "9", "7", "41", "100", "6", "24", "22"))

CL1.aob.colors <- c("purple", "mediumorchid2" , "mediumpurple", "purple4", "navy", "lavender")
CL1.stand.aob.mean.plot <- ggplot(CL1.stand.aob.mean, aes(fill=factor(CL1), y=PROP, x=STAND)) + geom_bar(position = "fill", stat = "identity") + ylab(paste0("Relative Abundance")) + xlab(paste0("Stand")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) +scale_fill_manual(values = CL1.aob.colors) + guides(fill=guide_legend(title="Cluster\n0.1"))

ggsave("Barplot_Cluster_1_AOB.pdf", width = 8.5, height = 10)

#Sideways barplot with Lineages separated out
boxplot.colors.aob  <- c("coral4", "coral3", "pink2", "pink", "purple", "dodgerblue2", "dodgerblue4", "darkblue")
CL1.stand.aob.mean.group.plot <- ggplot(CL1.stand.aob.mean, aes(fill=STAND, y=factor(CL1), x=PROP)) + geom_bar(position = position_dodge(), stat = "identity") + xlab(paste0("Relative Abundance")) + ylab(paste0("Cluster")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +scale_fill_manual(values = boxplot.colors.aob) + ggtitle("Relative Cluster Abundance in each Stand") + theme(plot.title = element_text(size = 20, hjust=0.5), axis.text=element_text(size = 18), axis.title = element_text(size =16), legend.text=element_text(size=14), legend.position = c(0.8,0.5))

ggsave("Barplot_Cluster_Groups_1_AOB.pdf", width = 8.5, height = 10)
