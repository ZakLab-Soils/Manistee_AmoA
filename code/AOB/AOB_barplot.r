#Create and analyze barplot for Cluster groupings
library(tidyverse)
library(ggplot2)


phy.aob.pruned <- readRDS("Phyloseq_AOB_Pruned.rds")

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

ggsave("Barplot_Cluster_P1_AOB.pdf", width = 8.5, height = 10)

#Sideways barplot with TAX2 Clades separated out
boxplot.colors.aob  <- c("coral4", "coral3", "pink2", "pink", "purple", "dodgerblue2", "dodgerblue4", "darkblue")
CL1.stand.aob.mean.group.plot <- ggplot(CL1.stand.aob.mean, aes(fill=STAND, y=factor(CL1), x=PROP)) + geom_bar(position = position_dodge(), stat = "identity") + xlab(paste0("Relative Abundance")) + ylab(paste0("Cluster")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +scale_fill_manual(values = boxplot.colors.aob) + ggtitle("Relative Cluster Abundance in each Stand") + theme(plot.title = element_text(size = 20, hjust=0.5), axis.text=element_text(size = 18), axis.title = element_text(size =16), legend.text=element_text(size=14), legend.position = c(0.8,0.5))

ggsave("Barplot_Cluster_Groups_P1_AOB.pdf", width = 8.5, height = 10)

#Examining Cluster 0.05
CL05.aob.tbl <-  Cl.all.aob.tbl %>% group_by(CL05, STAND, PLOT) %>% summarize(COUNT = sum(COUNT))%>% ungroup()%>% select(CL05, STAND, PLOT, COUNT)

CL05.aob.tmp.tbl <- CL05.aob.tbl  %>% group_by(CL05) %>% mutate(CL05.TOTAL = sum(COUNT))%>% ungroup() %>% filter(CL05.TOTAL > 0) %>% group_by(PLOT) %>% mutate(PLOT.COUNT = sum(COUNT))%>% ungroup() %>% group_by(CL05, PLOT) %>% mutate(PROP = sum(COUNT)/PLOT.COUNT) %>% ungroup() %>% mutate(HELLINGER = sqrt(PROP))

#At the 0.05 pairwise clustering, all AOB ASVs collapse into 14 clusters
CL05.stand.aob.mean <- aggregate(PROP~CL05+STAND, CL05.aob.tmp.tbl, mean)
CL05.stand.aob.mean$STAND <- gsub("Stand_*", "", fixed=FALSE, CL05.stand.aob.mean$STAND)
CL05.stand.aob.mean <- CL05.stand.aob.mean[order(CL05.stand.aob.mean$STAND), ]
CL05.stand.aob.mean$STAND <- factor(CL05.stand.aob.mean$STAND, levels = c("58", "9", "7", "41", "100", "6", "24", "22"))

CL05.aob.colors <- c("orchid2", "orchid", "orchid3", "mediumorchid1", "mediumorchid2", "mediumorchid3", "mediumorchid", "mediumpurple1", "mediumpurple2", "mediumpurple3", "mediumpurple", "purple", "purple3", "lavender")

CL05.stand.aob.mean.plot <- ggplot(CL05.stand.aob.mean, aes(fill=factor(CL05), y=PROP, x=STAND)) + geom_bar(position = "fill", stat = "identity") + ylab(paste0("Relative Abundance")) + xlab(paste0("Stand")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) +scale_fill_manual(values = CL05.aob.colors) + guides(fill=guide_legend(title="Cluster\n0.05"))

#Examining Cluster 0.03
CL03.aob.tbl <-  Cl.all.aob.tbl %>% group_by(CL03, STAND, PLOT) %>% summarize(COUNT = sum(COUNT))%>% ungroup()%>% select(CL03, STAND, PLOT, COUNT)

CL03.aob.tmp.tbl <- CL03.aob.tbl  %>% group_by(CL03) %>% mutate(CL03.TOTAL = sum(COUNT))%>% ungroup() %>% filter(CL03.TOTAL > 0) %>% group_by(PLOT) %>% mutate(PLOT.COUNT = sum(COUNT))%>% ungroup() %>% group_by(CL03, PLOT) %>% mutate(PROP = sum(COUNT)/PLOT.COUNT) %>% ungroup() %>% mutate(HELLINGER = sqrt(PROP))

#At the 0.03 pairwise clustering, all AOB ASVs collapse into 14 clusters
CL03.stand.aob.mean <- aggregate(PROP~CL03+STAND, CL03.aob.tmp.tbl, mean)
CL03.stand.aob.mean$STAND <- gsub("Stand_*", "", fixed=FALSE, CL03.stand.aob.mean$STAND)
CL03.stand.aob.mean <- CL03.stand.aob.mean[order(CL03.stand.aob.mean$STAND), ]
CL03.stand.aob.mean$STAND <- factor(CL03.stand.aob.mean$STAND, levels = c("58", "9", "7", "41", "100", "6", "24", "22"))

#24 different clusters
CL03.aob.colors <- c("lavender", "pink", "lightpink", "pink2", "lightpink2", "orchid2", "orchid", "orchid3", "magenta3", "mediumorchid1", "mediumorchid2", "mediumorchid3", "mediumorchid","orchid4","mediumorchid4", "magenta4","mediumpurple1", "mediumpurple2", "mediumpurple3", "mediumpurple", "purple", "purple2" , "purple3", "purple4", "navy")

CL03.stand.aob.mean.plot <- ggplot(CL03.stand.aob.mean, aes(fill=factor(CL03), y=PROP, x=STAND)) + geom_bar(position = "fill", stat = "identity") + ylab(paste0("Relative Abundance")) + xlab(paste0("Stand")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) +scale_fill_manual(values = CL03.aob.colors) + guides(fill=guide_legend(title="Cluster\n0.03"))






###STILL EDITING - NEED TO DECIDE WHICH LEVEL TO USE

#Test for significance for each Clade with Kruskal test and Pairwise-Wilcox
TAX2.plot.stand.aoa.mean <- aggregate(PROP~TAX2+PLOT+STAND, TAX2.aoa.tmp.tbl, mean)

tax2.nsbeta.aoa <- TAX2.plot.stand.aoa.mean[TAX2.plot.stand.aoa.mean$TAX2 == "NS-Beta",]
tax2.nsdelta.aoa <- TAX2.plot.stand.aoa.mean[TAX2.plot.stand.aoa.mean$TAX2 == "NS-Delta",]
tax2.nsalpha.aoa <- TAX2.plot.stand.aoa.mean[TAX2.plot.stand.aoa.mean$TAX2 == "NS-Alpha",]
tax2.nsgamma.aoa <- TAX2.plot.stand.aoa.mean[TAX2.plot.stand.aoa.mean$TAX2 == "NS-Gamma",]
tax2.ntalpha.aoa <- TAX2.plot.stand.aoa.mean[TAX2.plot.stand.aoa.mean$TAX2 == "NT-Alpha",]

ktnsa.aoa <- kruskal.test(PROP~STAND, data = tax2.nsalpha.aoa)
ktnta.aoa <- kruskal.test(PROP~STAND, data = tax2.ntalpha.aoa)
ktnsg.aoa <- kruskal.test(PROP~STAND, data = tax2.nsgamma.aoa)
ktnsd.aoa <- kruskal.test(PROP~STAND, data = tax2.nsdelta.aoa)
ktnsb.aoa <- kruskal.test(PROP~STAND, data = tax2.nsbeta.aoa)

ptnsa.aoa <- pairwise.wilcox.test(tax2.nsalpha.aoa$PROP, g=tax2.nsalpha.aoa$STAND, p.adjust.method = "BH")
ptnsg.aoa  <- pairwise.wilcox.test(tax2.nsgamma.aoa$PROP, g=tax2.nsgamma.aoa$STAND, p.adjust.method = "BH")
ptnsd.aoa  <- pairwise.wilcox.test(tax2.nsdelta.aoa$PROP, g=tax2.nsdelta.aoa$STAND, p.adjust.method = "BH")
ptnsb.aoa  <- pairwise.wilcox.test(tax2.nsbeta.aoa$PROP, g=tax2.nsbeta.aoa$STAND, p.adjust.method = "BH")
ptnta.aoa  <- pairwise.wilcox.test(tax2.ntalpha.aoa$PROP, g=tax2.ntalpha.aoa$STAND, p.adjust.method = "BH")

myresults <- list(ktnsa.aoa, ktnsb.aoa, ktnsd.aoa, ktnsg.aoa, ktnta.aoa, ptnsa.aoa, ptnsb.aoa, ptnsd.aoa, ptnsg.aoa, ptnta.aoa)
capture.output(myresults, file = "Kruskal_Wilcox_Results_Cluster1.txt")
