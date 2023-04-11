#Create and analyze barplot for TAX groups
library(tidyverse)
library(ggplot2)

phy.aoa.pruned <- readRDS("Phyloseq_AOA_Pruned.rds")

phy.aoa.asv.df <- data.frame(otu_table(phy.aoa.pruned)) %>% as_tibble(., rownames="SAMPLE.ID")
phy.aoa.samdata.df <- data.frame(sample_data(phy.aoa.pruned)) %>%  as_tibble(., rownames = "SAMPLE.ID")
phy.aoa.tax.df <- data.frame(tax_table(phy.aoa.pruned)) %>% as_tibble(., rownames = "ASV") %>% rename_all(toupper) %>% pivot_longer(-ASV, names_to = "RANKING", values_to = "TAXON") %>% pivot_wider(id_cols = ASV, names_from = "RANKING", values_from = TAXON)
names(phy.aoa.tax.df)[] <- c("ASV", "TAX1", "TAX2", "TAX3", "TAX4", "TAX5", "TAX6", "TAX7", "TAX8", "TAX9", "TAX10")


ASV.aoa.tbl <- phy.aoa.asv.df %>% 
pivot_longer(-SAMPLE.ID, names_to = "ASV", values_to = "COUNT") %>%
inner_join(phy.aoa.samdata.df, ., by = "SAMPLE.ID") %>% 
group_by(PLOT, STAND, ASV) %>%
summarise(COUNT = sum(COUNT)) %>%
ungroup() %>%
group_by(ASV) %>%
mutate(ASV.TOTAL = sum(COUNT)) %>%
ungroup() %>%
filter(ASV.TOTAL > 0) %>%
select(-ASV.TOTAL)

TAX.all.aoa.tbl <- ASV.aoa.tbl %>% inner_join(phy.aoa.tax.df, ., by = "ASV")

#Will use these again for TITAN analysis later
saveRDS(ASV.aoa.tbl, file="ASV_AOA_Table.rds")
saveRDS(TAX.all.aoa.tbl, file = “TAX_ALL_AOA_Table.rds)

TAX2.aoa.tbl <-  TAX.all.aoa.tbl %>% group_by(TAX2, STAND, PLOT) %>% summarize(COUNT = sum(COUNT))%>% ungroup()%>% select(TAX2, STAND, PLOT, COUNT)

TAX2.aoa.tmp.tbl <- TAX2.aoa.tbl  %>% group_by(TAX2) %>% mutate(TAX2.TOTAL = sum(COUNT))%>% ungroup() %>% filter(TAX2.TOTAL > 0) %>% group_by(PLOT) %>% mutate(PLOT.COUNT = sum(COUNT))%>% ungroup() %>% group_by(TAX2, PLOT) %>% mutate(PROP = sum(COUNT)/PLOT.COUNT) %>% ungroup() %>% mutate(HELLINGER = sqrt(PROP))




#For all the TAX2 groups available
#TAX2.stand.aoa.mean <- aggregate(PROP~TAX2+STAND, TAX2.aoa.tmp.tbl, mean)
#TAX2.stand.aoa.mean$STAND <- gsub("Stand_*", "", fixed=FALSE, #TAX2.stand.aoa.mean$STAND)
#TAX2.stand.aoa.mean <- TAX2.stand.aoa.mean[order(TAX2.stand.aoa.mean$STAND), ]
#TAX2.stand.aoa.mean$STAND <- factor(TAX2.stand.aoa.mean$STAND, levels = c("58", "7", "41", "100", "24", "22", "6"))
#TAX2.aoa.colors <- c("purple", "mediumorchid2", "mediumpurple2", "purple3", "lavender", "navy")
#TAX2.stand.aoa.mean.plot <- ggplot(TAX2.stand.aoa.mean, aes(fill=TAX2, y=PROP, x=STAND)) + geom_bar(position = "fill", stat = "identity") + ylab(paste0("Relative Abundance")) + xlab(paste0("Stand")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) +scale_fill_manual(values = TAX2.aoa.colors)

#Removed NS-OTU1 since <0.05% of data
TAX2.stand.aoa.mean.lessNSotu1 <- aggregate(PROP~TAX2+STAND, TAX2.aoa.tmp.tbl[TAX2.aoa.tmp.tbl$TAX2 != "NS_OTU1",], mean)
TAX2.stand.aoa.mean.lessNSotu1$STAND <- gsub("Stand_*", "", fixed=FALSE, TAX2.stand.aoa.mean.lessNSotu1$STAND)
TAX2.stand.aoa.mean.lessNSotu1$STAND <- factor(TAX2.stand.aoa.mean.lessNSotu1$STAND, levels = c("58", "7", "41", "100", "24", "22", "6"))

TAX2.aoa.colors2 <- c("purple", "mediumorchid2", "mediumpurple2", "purple3","navy")

TAX2.stand.aoa.mean.lessNSotu1.plot <- ggplot(TAX2.stand.aoa.mean.lessNSotu1, aes(fill=TAX2, y=PROP, x=STAND)) + geom_bar(position = "fill", stat = "identity") + ylab(paste0("Relative Abundance")) + xlab(paste0("Stand")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) +scale_fill_manual(values = TAX2.aoa.colors2)

ggsave("Barplot_TAXON2_AOA.pdf", width = 8.5, height = 10)
