#These commands to ASV.all.tbl also in titan code
#ps.ASV.tbl <- as.data.frame(otu_table(phy.Rs)) %>% as_tibble(., rownames = "SAMPLE.ID")
#ps.sample.tbl <- as.data.frame(sample_data(phy.Rs)) %>% as_tibble(.)

#ASV.all.tbl <- ps.ASV.tbl %>%
#pivot_longer(-SAMPLE.ID, names_to = "ASV", values_to = "COUNT") %>%
#inner_join(ps.sample.tbl, ., by = "SAMPLE.ID") %>%
#filter(TYPE == "environmental") %>%
#group_by(PLOT, STAND, ASV) %>%
#summarise(COUNT = sum(COUNT)) %>%
#ungroup() %>%
#group_by(ASV) %>%
#mutate(ASV.TOTAL = sum(COUNT)) %>%
#ungroup() %>%
#filter(ASV.TOTAL > 0) %>%
#select(-ASV.TOTAL

TAX2.all.tbl <- ASV.all.tbl %>% inner_join(ps.tax.tbl, ., by = "ASV") %>% group_by(TAX2, STAND, PLOT) %>% summarize(COUNT = sum(COUNT))%>% ungroup()%>% select(TAX2, STAND, PLOT, COUNT)

TAX2.all.tmp.tbl <- TAX2.all.tbl  %>% group_by(TAX2) %>% mutate(TAX2.TOTAL = sum(COUNT))%>% ungroup() %>% filter(TAX2.TOTAL > 0) %>% group_by(PLOT) %>% mutate(PLOT.COUNT = sum(COUNT))%>% ungroup() %>% group_by(TAX2, PLOT) %>% mutate(PROP = sum(COUNT)/PLOT.COUNT) %>% ungroup() %>% mutate(HELLINGER = sqrt(PROP))

TAX2.stand.mean <- aggregate(PROP~TAX2+STAND, TAX2.all.tmp.tbl, mean)
TAX2.stand.mean$STAND <- gsub("Stand_*", "", fixed=FALSE, TAX2.stand.mean$STAND)
TAX2.stand.mean <- TAX2.stand.mean[order(TAX2.stand.mean$STAND), ]

TAX2.colors <- c("purple", "mediumorchid2", "mediumpurple2", "purple3", "lavender", "navy")
TAX2.stand.mean$STAND <- factor(TAX2.stand.mean$STAND, levels = c("58", "7", "6", "41", "24", "100", "22"))
ggplot(TAX2.stand.mean, aes(fill=TAX2, y=PROP, x=STAND)) + geom_bar(position = "fill", stat = "identity") + ylab(paste0("Relative Abundance")) + xlab(paste0("Stand")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) +scale_fill_manual(values = TAX2.colors)

#Removed NS-OTU1 since <0.05% of data
TAX2.stand.mean.nonsout1 <- aggregate(PROP~TAX2+STAND, TAX2.all.tmp.tbl[TAX2.all.tmp.tbl$TAX2 != "NS_OTU1",], mean)
TAX2.stand.mean.nonsout1$STAND <- gsub("Stand_*", "", fixed=FALSE, TAX2.stand.mean.nonsout1$STAND)
TAX2.stand.mean.nonsout1$STAND <- factor(TAX2.stand.mean.nonsout1$STAND, levels = c("58", "7", "6", "41", "24", "100", "22"))
ggplot(TAX2.stand.mean.nonsout1, aes(fill=TAX2, y=PROP, x=STAND)) + geom_bar(position = "fill", stat = "identity") + ylab(paste0("Relative Abundance")) + xlab(paste0("Stand")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) +scale_fill_manual(values = TAX2.colors.nonsotu1)
