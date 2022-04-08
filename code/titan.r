
ps.tax.tbl <- as.data.frame(tax_table(phy.Rs)) %>% as_tibble(., rownames = "ASV") %>% rename_all(toupper) %>% pivot_longer(-ASV, names_to = "RANKING", values_to = "TAXON") %>% pivot_wider(id_cols = ASV, names_from = "RANKING", values_from = TAXON)

ps.ASV.tbl <- as.data.frame(otu_table(phy.Rs)) %>% as_tibble(., rownames = "SAMPLE.ID")

ps.sample.tbl <- as.data.frame(sample_data(phy.Rs)) %>% as_tibble(.)

ASV.all.tbl <- ps.ASV.tbl %>%
pivot_longer(-SAMPLE.ID, names_to = "ASV", values_to = "COUNT") %>%
inner_join(ps.sample.tbl, ., by = "SAMPLE.ID") %>%
filter(TYPE == "environmental") %>%
group_by(PLOT, STAND, ASV) %>%
summarise(COUNT = sum(COUNT)) %>%
ungroup() %>%
group_by(ASV) %>%
mutate(ASV.TOTAL = sum(COUNT)) %>%
ungroup() %>%
filter(ASV.TOTAL > 0) %>%
select(-ASV.TOTAL)

names(ps.tax.tbl)[] <- c("ASV", "TAX1", "TAX2", "TAX3", "TAX4", "TAX5", "TAX6", "TAX7", "TAX8", "TAX9", "TAX10")

TAX10.all.tbl <- ASV.all.tbl %>% inner_join(ps.tax.tbl, ., by = "ASV") %>% group_by(TAX10, STAND, PLOT) %>% summarize(COUNT = sum(COUNT))%>% ungroup()%>% select(TAX10, STAND, PLOT, COUNT)

TAX10.all.perc.tbl <- TAX10.all.tbl  %>% group_by(TAX10) %>% summarise(TAX10.COUNT = sum(COUNT)) %>% ungroup() %>% mutate(TOTAL.COUNT = sum(TAX10.COUNT)) %>% mutate(TAX10.PERC = 100 *(TAX10.COUNT/TOTAL.COUNT))%>% select(TAX10, TAX10.PERC)

TAX10.all.tmp.tbl <- TAX10.all.tbl  %>% group_by(TAX10) %>% mutate(TAX10.TOTAL = sum(COUNT))%>% ungroup() %>% filter(TAX10.TOTAL > 0) %>% group_by(PLOT) %>% mutate(PLOT.COUNT = sum(COUNT))%>% ungroup() %>% group_by(TAX10, PLOT) %>% mutate(PROP = sum(COUNT)/PLOT.COUNT) %>% ungroup() %>% mutate(HELLINGER = sqrt(PROP))

TAX10.all.trim.tbl <- TAX10.all.tbl %>% mutate(PA = ifelse(COUNT > 0,1,0))%>% group_by(TAX10) %>% inner_join(., TAX10.all.perc.tbl, by = "TAX10") %>% mutate(OCCURRENCE = sum(PA))%>% ungroup() %>% filter(TAX10.PERC >= 0.1 & OCCURRENCE >= 3)

TAX10.all.trim.Occ5.tbl <- TAX10.all.tbl %>% mutate(PA = ifelse(COUNT > 0,1,0))%>% group_by(TAX10) %>% inner_join(., TAX10.all.perc.tbl, by = "TAX10") %>% mutate(OCCURRENCE = sum(PA))%>% ungroup() %>% filter(TAX10.PERC >= 0.1 & OCCURRENCE >= 5)

env.N.min.df <- soil.compiled.data.tbl %>% select(PLOT, N.MIN) %>% arrange(PLOT) %>% column_to_rownames(var = "PLOT") %>% as.data.frame(.)

env.ph.df <- soil.ph.tbl %>% select(PLOT, SOIL.PH) %>% arrange(PLOT) %>% column_to_rownames(var = "PLOT") %>% as.data.frame(.)

aoa.plots <- unique(ASV.all.tbl$PLOT)

env.N.min.df2 <- env.N.min.df %>% filter(row.names(env.N.min.df) %in% aoa.plots)

env.ph.df2 <- env.ph.df %>% filter(row.names(env.ph.df) %in% aoa.plots)

phy.Rs.hel <- decostand(phyloseq::otu_table(phy.Rs), method = "hellinger")

TAX10.all.trim.titan.in.df <- TAX10.all.trim.tbl %>% select(TAX10) %>% distinct() %>% inner_join(., TAX10.all.tmp.tbl, by = "TAX10") %>% select(PLOT, TAX10, HELLINGER) %>% pivot_wider(id_cols = PLOT, names_from = "TAX10", values_from = "HELLINGER") %>% arrange(PLOT) %>% column_to_rownames(var = "PLOT") %>% as.data.frame(.)

TAX10.all.trim.titan <- titan(env.N.min.df2, TAX10.all.trim.titan.in.df, minSplt =5, numPerm = 1000, boot = TRUE, nBoot = 1000, imax = FALSE, ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95, ncpus =4, memory = TRUE)

TAX10.group.titan.tbl <- TAX10.TAX2.tbl %>% mutate(TAX2.GROUP = ifelse(TAX2 == "NT-Alpha", "NT-Alpha", NA)) %>% mutate(TAX2.GROUP = ifelse(TAX2 == "NS-Gamma", "NS-Gamma", TAX2.GROUP)) %>% mutate(TAX2.GROUP = ifelse(TAX2 == "NS-Delta", "NS-Delta", TAX2.GROUP)) %>% mutate(TAX2.GROUP = ifelse(TAX2 == "NS-Beta", "NS-Beta", TAX2.GROUP)) %>% mutate(TAX2.GROUP = ifelse(TAX2 == "NS-Alpha", "NS-Alpha", TAX2.GROUP)) %>% select(TAX10, TAX2.GROUP) %>% distinct()

TAX10.titan.response.tbl <- as_tibble(TAX10.al.trim.titan$sppmax, rownames = NA) %>% rownames_to_column(var = "TAX10") %>% select (TAX10, zenv.cp, '5%', '95%', z.median, filter, maxgrp) %>% rename(ZENV.CP = "zenv.cp", LCI = "5%", UCI = "95%", Z.MEDIAN = "z.median", GROUP = "filter", MAXGRP = "maxgrp") %>% mutate(Z.MEDIAN = ifelse(MAXGRP == 1, -1 * Z.MEDIAN, Z.MEDIAN)) %>% inner_join(.,TAX10.group.titan.tbl, by = "TAX10") %>% mutate(TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NT-Alpha", "NT-Alpha", NA), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Gamma", "NS-Gamma", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Delta", "NS-Delta", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Beta", "NS-Beta", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Alpha", "NS-Alpha", TAX2.GROUP.LABEL)) %>% mutate(TAX2.GROUP.LABEL = factor(TAX2.GROUP.LABEL, levels = c("NT-Alpha", "NS-Gamma", "NS-Delta", "NS-Beta", "NS-Alpha"))) %>% arrange(TAX2.GROUP.LABEL, Z.MEDIAN) %>% mutate(TAX10 = factor(TAX10, levels = TAX10)) %>% mutate(SIG = ifelse(GROUP >0, "Significant", "Not Significant")) %>% mutate(SIG = factor(SIG, levels = c("Significant", "Not Significant")))

TAX10.fig.input.tbl <- as_tibble(TAX10.al.trim.titan$sppmax, rownames = "TAX10") %>% inner_join(., TAX10.group.titan.tbl, by = "TAX10") %>% dplyr::rename(RESPONSE.GROUP = "filter")  %>% select(TAX10, TAX2.GROUP, RESPONSE.GROUP) %>% distinct()  %>% group_by(TAX2.GROUP)%>% summarise(POSITIVE.RESPONSE = length(TAX2.GROUP[RESPONSE.GROUP == 2]), NEGATIVE.RESPONSE = length(TAX2.GROUP[RESPONSE.GROUP ==1]), NO.RESPONSE = length(TAX2.GROUP[RESPONSE.GROUP ==0]), TOTAL.TAX10 = length(TAX2.GROUP)) %>% ungroup() %>% mutate(TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NT-Alpha", "NT-Alpha", NA), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Gamma", "NS-Gamma", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Delta", "NS-Delta", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Beta", "NS-Beta", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Alpha", "NS-Alpha", TAX2.GROUP.LABEL)) %>% mutate(TAX2.GROUP.LABEL = factor(TAX2.GROUP.LABEL, levels = c("NT-Alpha", "NS-Gamma", "NS-Delta", "NS-Beta", "NS-Alpha")))

TAX10.titan.plot <- ggplot() +   geom_bar(data = TAX10.fig.input.tbl, aes(x = TAX2.GROUP.LABEL, y = PROPORTION, fill = RESPONSE),colour = NA, stat = "identity", alpha = 0.75) +   geom_hline(yintercept = 0,linetype = 2,size = 0.5,colour = "black") + geom_text(data = TAX10.fig.input.tbl, aes(x = TAX2.GROUP.LABEL, y = LABEL.POSITION, colour = RESPONSE), label = TAX10.fig.input.tbl$SIG.LABEL, size = 2.5) + scale_colour_viridis(name = "Association with\ninorganic N", labels = c("Negative", "Positive"), discrete = TRUE, begin = 0.2, end = 0.8, option = "plasma", guide = FALSE) + scale_fill_viridis(name = "Association with\ninorgaic N", discrete = TRUE, labels = c("Negative", "Positive"), begin = 0.2, end = 0.8, option = "plasma") +   scale_y_continuous(limits = c(-1.2, 1.2), expand = expansion(mult = c(0.1, 0.1)), breaks = c(-1, -0.5, 0, 0.5, 1),labels = c("1.0", "0.5", "0", "0.5", "1.0")) + coord_flip() + scale_x_discrete(limits = rev(levels(TAX10.fig.input.tbl$TAX2.GROUP.LABEL))) + labs(x = NULL, y = "Proportion of ASV in Clade") + theme(axis.line = element_line(colour = "black", size = 0.5),panel.background = element_rect(fill = NA), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title = element_text(size = 12, colour = "black", hjust = -0.35),axis.title.x = element_text(size = 7, colour = "black"), axis.title.y = element_text(size = 7, colour = "black"), axis.text.x = element_text(size = 6, colour = "black"), axis.text.y = element_text(size = 6, colour = "black"), legend.title = element_text(size = 7, colour = "black"), legend.text = element_text(size = 6, colour = "black"), legend.key = element_blank())

TAX10.titan.response.plot <- ggplot()+geom_bar(data=TAX10.titan.response.tbl, aes(x=TAX10, y=Z.MEDIAN, fill = TAX2.GROUP.LABEL, alpha = SIG), stat = "identity", width = 0.75) + scale_fill_viridis_d(name = NULL, option = "plasma", begin =0.1, end = 0.9) + scale_alpha_manual(name = NULL, values = c(1, 0.3)) + geom_hline(yintercept = 0, linetype = 2, size = 0.5, colour = "black") + labs(x="ASV Taxon", y = "Median Z-score\n(response to inorganic N availabilty)") + coord_flip() + scale_x_discrete(limits=rev(levels(TAX10.titan.response.tbl$TAX10))) + theme(axis.line = element_line(colour = "black", size = 0.5), panel.background = element_rect(fill = NA), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title = element_text(size = 12, colour = "black", hjust = -0.5),axis.title.x = element_text(size = 7, colour = "black"), axis.title.y = element_text(size = 7, colour = "black"), axis.text.x = element_text(size = 6, colour = "black"), axis.text.y = element_text(size = 6, colour = "black"), legend.title = element_text(size = 7, colour = "black"), legend.text = element_text(size = 6, colour = "black"), legend.key = element_blank())

TAX10.all.trim.ph.titan <- titan(env.ph.df2, TAX10.all.trim.titan.in.df, minSplt =5, numPerm = 1000, boot = TRUE, nBoot = 1000, imax = FALSE, ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95, ncpus =4, memory = TRUE)
