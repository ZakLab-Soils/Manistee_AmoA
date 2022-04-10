
TAX2.all.tbl <- ASV.all.tbl %>% inner_join(ps.tax.tbl, ., by = "ASV") %>% group_by(TAX2, STAND, PLOT) %>% summarize(COUNT = sum(COUNT))%>% ungroup()%>% select(TAX2, STAND, PLOT, COUNT)

TAX2.all.tmp.tbl <- TAX2.all.tbl  %>% group_by(TAX2) %>% mutate(TAX2.TOTAL = sum(COUNT))%>% ungroup() %>% filter(TAX2.TOTAL > 0) %>% group_by(PLOT) %>% mutate(PLOT.COUNT = sum(COUNT))%>% ungroup() %>% group_by(TAX2, PLOT) %>% mutate(PROP = sum(COUNT)/PLOT.COUNT) %>% ungroup() %>% mutate(HELLINGER = sqrt(PROP))

TAX2.stand.mean <- aggregate(PROP~TAX2+STAND, TAX2.all.tmp.tbl, mean)
TAX2.stand.mean$STAND <- gsub("Stand_*", "", fixed=FALSE, TAX2.stand.mean$STAND)
TAX2.stand.mean <- TAX2.stand.mean[order(TAX2.stand.mean$STAND), ]

TAX2.colors <- c("purple", "mediumorchid2", "mediumpurple2", "purple3", "lavender", "navy")

ggplot(TAX2.stand.mean, aes(fill=TAX2, y=PROP, x=STAND)) + geom_bar(position = "fill", stat = "identity") + ylab(paste0("Relative Abundance")) + xlab(paste0("Stand")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) +scale_fill_manual(values = TAX2.colors)
