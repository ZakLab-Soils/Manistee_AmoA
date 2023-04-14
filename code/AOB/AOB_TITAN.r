library(tidyverse)
library(TITAN2)
library(ggplot2)
library(virdisLite)

phy.aob.pruned <- readRDS("Phyloseq_AOB_Pruned.rds")
ASV.aob.tbl <- readRDS("ASV_AOB_Table.rds")
CL.all.aob.tbl <- readRDS("CL_ALL_AOB_Table.rds")

###STILL WORKING TO GET DATA INTO CORRECT FORMAT

#Titan

CL01.tbl <- CL.all.aob.tbl %>% group_by(CL01, STAND, PLOT) %>% summarise(COUNT = sum(COUNT)) %>% ungroup() %>% select(CL01, STAND, PLOT, COUNT)
CL03.tbl <- CL.all.aob.tbl %>% group_by(CL03, STAND, PLOT) %>% summarise(COUNT = sum(COUNT)) %>% ungroup() %>% select(CL03, STAND, PLOT, COUNT)
CL1.tbl  <- CL.all.aob.tbl %>% group_by(CL1, STAND, PLOT) %>% summarise(COUNT = sum(COUNT)) %>% ungroup() %>% select(CL1, STAND, PLOT, COUNT)

CL01.perc.tbl <- CL01.tbl  %>% group_by(CL01) %>% summarise(CL01.COUNT = sum(COUNT)) %>% ungroup() %>% mutate(TOTAL.COUNT = sum(CL01.COUNT)) %>% mutate(CL01.PERC = 100 *(CL01.COUNT/TOTAL.COUNT))%>% select(CL01, CL01.PERC)
CL03.perc.tbl <- CL03.tbl  %>% group_by(CL03) %>% summarise(CL03.COUNT = sum(COUNT)) %>% ungroup() %>% mutate(TOTAL.COUNT = sum(CL03.COUNT)) %>% mutate(CL03.PERC = 100 *(CL03.COUNT/TOTAL.COUNT))%>% select(CL03, CL03.PERC)

CL01.aob.tmp.tbl <- CL01.tbl  %>% group_by(CL01) %>% mutate(CL01.TOTAL = sum(COUNT))%>% ungroup() %>% filter(CL01.TOTAL > 0) %>% group_by(PLOT) %>% mutate(PLOT.COUNT = sum(COUNT))%>% ungroup() %>% group_by(CL01, PLOT) %>% mutate(PROP = sum(COUNT)/PLOT.COUNT) %>% ungroup() %>% mutate(HELLINGER = sqrt(PROP))
CL03.aob.tmp.tbl <- CL03.tbl  %>% group_by(CL03) %>% mutate(CL03.TOTAL = sum(COUNT))%>% ungroup() %>% filter(CL03.TOTAL > 0) %>% group_by(PLOT) %>% mutate(PLOT.COUNT = sum(COUNT))%>% ungroup() %>% group_by(CL03, PLOT) %>% mutate(PROP = sum(COUNT)/PLOT.COUNT) %>% ungroup() %>% mutate(HELLINGER = sqrt(PROP))

CL01.aob.trim.tbl <- CL01.tbl %>% mutate(PA = ifelse(COUNT > 0,1,0))%>% group_by(CL01) %>% inner_join(., CL01.perc.tbl, by = "CL01") %>% mutate(OCCURRENCE = sum(PA))%>% ungroup() %>% filter(CL01.PERC >= 0.1 & OCCURRENCE >= 3)
CL03.aob.trim.tbl <- CL03.tbl %>% mutate(PA = ifelse(COUNT > 0,1,0))%>% group_by(CL03) %>% inner_join(., CL03.perc.tbl, by = "CL03") %>% mutate(OCCURRENCE = sum(PA))%>% ungroup() %>% filter(CL03.PERC >= 0.1 & OCCURRENCE >= 3)

CL01.CL1.tbl <- CL.all.aob.tbl %>% group_by(CL01,CL1, STAND, PLOT) %>% summarize(COUNT = sum(COUNT))%>% ungroup()%>% select(CL01, CL1, STAND, PLOT, COUNT)
CL03.CL1.tbl <- CL.all.aob.tbl %>% group_by(CL03,CL1, STAND, PLOT) %>% summarize(COUNT = sum(COUNT))%>% ungroup()%>% select(CL03, CL1, STAND, PLOT, COUNT)

#Dropped Cluster 1-6 (CL 0.1, number 6)
CL01.group.titan.aob.tbl <- CL01.CL1.tbl %>% mutate(CL1.GROUP = ifelse(CL1 == "1", "Cluster 1-1", NA)) %>% mutate(CL1.GROUP = ifelse(CL1 == "2", "Cluster 1-2", CL1.GROUP)) %>% mutate(CL1.GROUP = ifelse(CL1 == "3", "Cluster 1-3", CL1.GROUP)) %>% mutate(CL1.GROUP = ifelse(CL1 == "4", "Cluster 1-4", CL1.GROUP)) %>% mutate(CL1.GROUP = ifelse(CL1 == "5", "Cluster 1-5", CL1.GROUP)) %>% select(CL01, CL1.GROUP) %>% distinct()

CL03.group.titan.aob.tbl <- CL03.CL1.tbl %>% mutate(CL1.GROUP = ifelse(CL1 == "1", "Cluster 1-1", NA)) %>% mutate(CL1.GROUP = ifelse(CL1 == "2", "Cluster 1-2", CL1.GROUP)) %>% mutate(CL1.GROUP = ifelse(CL1 == "3", "Cluster 1-3", CL1.GROUP)) %>% mutate(CL1.GROUP = ifelse(CL1 == "4", "Cluster 1-4", CL1.GROUP)) %>% mutate(CL1.GROUP = ifelse(CL1 == "5", "Cluster 1-5", CL1.GROUP)) %>% select(CL03, CL1.GROUP) %>% distinct()

phy.aob.samdata.df <- data.frame(sample_data(phy.aob.pruned)) %>%  as_tibble(., rownames = "SAMPLE.ID")


env.nmin.aob.df <- phy.aob.samdata.df %>% select(PLOT, N.MIN) %>% arrange(PLOT) %>% column_to_rownames(var = "PLOT") %>% as.data.frame(.)
env.ph.aob.df <- phy.aob.samdata.df %>% select(PLOT, SOIL.PH) %>% arrange(PLOT) %>% column_to_rownames(var = "PLOT") %>% as.data.frame(.)
env.ammpost.aob.df <- phy.aob.samdata.df %>% select(PLOT, AMM.CORRECTED.POST) %>% arrange(PLOT) %>% column_to_rownames(var = "PLOT") %>% as.data.frame(.)

phy.aob.hel <- decostand(phyloseq::otu_table(phy.aob.pruned), method = "hellinger")

CL01.aob.trim.titan.in.df <- CL01.aob.trim.tbl %>% select(CL01) %>% distinct() %>% inner_join(., CL01.aob.tmp.tbl, by = "CL01") %>% select(PLOT, CL01, HELLINGER) %>% pivot_wider(id_cols = PLOT, names_from = "CL01", values_from = "HELLINGER") %>% arrange(PLOT) %>% column_to_rownames(var = "PLOT") %>% as.data.frame(.)
CL03.aob.trim.titan.in.df <- CL03.aob.trim.tbl %>% select(CL03) %>% distinct() %>% inner_join(., CL03.aob.tmp.tbl, by = "CL03") %>% select(PLOT, CL03, HELLINGER) %>% pivot_wider(id_cols = PLOT, names_from = "CL03", values_from = "HELLINGER") %>% arrange(PLOT) %>% column_to_rownames(var = "PLOT") %>% as.data.frame(.)

CL01.aob.trim.titan.nmin <- titan(env.nmin.aob.df, CL01.aob.trim.titan.in.df, minSplt =5, numPerm = 1000, boot = TRUE, nBoot = 1000, imax = FALSE, ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95, ncpus =4, memory = FALSE)
CL03.aob.trim.titan.nmin <- titan(env.nmin.aob.df, CL03.aob.trim.titan.in.df, minSplt =5, numPerm = 1000, boot = TRUE, nBoot = 1000, imax = FALSE, ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95, ncpus =4, memory = FALSE)

CL01.aob.trim.titan.ph <- titan(env.ph.aob.df, CL01.aob.trim.titan.in.df, minSplt =5, numPerm = 1000, boot = TRUE, nBoot = 1000, imax = FALSE, ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95, ncpus =4, memory = FALSE)
CL03.aob.trim.titan.ph <- titan(env.ph.aob.df, CL03.aob.trim.titan.in.df, minSplt =5, numPerm = 1000, boot = TRUE, nBoot = 1000, imax = FALSE, ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95, ncpus =4, memory = FALSE)

CL01.aob.trim.titan.ammpost <- titan(env.ammpost.aob.df, CL01.aob.trim.titan.in.df, minSplt =5, numPerm = 1000, boot = TRUE, nBoot = 1000, imax = FALSE, ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95, ncpus =4, memory = FALSE)
CL03.aob.trim.titan.ammpost <- titan(env.ammpost.aob.df, CL03.aob.trim.titan.in.df, minSplt =5, numPerm = 1000, boot = TRUE, nBoot = 1000, imax = FALSE, ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95, ncpus =4, memory = FALSE)


#Maybe some of the TITAN plots
#plot_sumz_density(CL01.aob.trim.titan.nmin, ribbon=FALSE, points = TRUE)





TAX10.titan.response.aoa.nmin.tbl <- as_tibble(TAX10.aoa.trim.titan.nmin$sppmax, rownames = NA) %>% rownames_to_column(var = "TAX10") %>% select (TAX10, zenv.cp, '5%', '95%', z.median, filter, maxgrp) %>% rename(ZENV.CP = "zenv.cp", LCI = "5%", UCI = "95%", Z.MEDIAN = "z.median", GROUP = "filter", MAXGRP = "maxgrp") %>% mutate(Z.MEDIAN = ifelse(MAXGRP == 1, -1 * Z.MEDIAN, Z.MEDIAN)) %>% inner_join(.,TAX10.group.titan.tbl, by = "TAX10") %>% mutate(TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NT-Alpha", "NT-Alpha", NA), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Gamma", "NS-Gamma", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Delta", "NS-Delta", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Beta", "NS-Beta", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Alpha", "NS-Alpha", TAX2.GROUP.LABEL)) %>% mutate(TAX2.GROUP.LABEL = factor(TAX2.GROUP.LABEL, levels = c("NT-Alpha", "NS-Gamma", "NS-Delta", "NS-Beta", "NS-Alpha"))) %>% arrange(TAX2.GROUP.LABEL, Z.MEDIAN) %>% mutate(TAX10 = factor(TAX10, levels = TAX10)) %>% mutate(SIG = ifelse(GROUP >0, "Significant", "Not Significant")) %>% mutate(SIG = factor(SIG, levels = c("Significant", "Not Significant")))

TAX10.titan.response.aoa.ph.tbl <- as_tibble(TAX10.aoa.trim.titan.ph$sppmax, rownames = NA) %>% rownames_to_column(var = "TAX10") %>% select (TAX10, zenv.cp, '5%', '95%', z.median, filter, maxgrp) %>% rename(ZENV.CP = "zenv.cp", LCI = "5%", UCI = "95%", Z.MEDIAN = "z.median", GROUP = "filter", MAXGRP = "maxgrp") %>% mutate(Z.MEDIAN = ifelse(MAXGRP == 1, -1 * Z.MEDIAN, Z.MEDIAN)) %>% inner_join(.,TAX10.group.titan.tbl, by = "TAX10") %>% mutate(TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NT-Alpha", "NT-Alpha", NA), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Gamma", "NS-Gamma", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Delta", "NS-Delta", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Beta", "NS-Beta", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Alpha", "NS-Alpha", TAX2.GROUP.LABEL)) %>% mutate(TAX2.GROUP.LABEL = factor(TAX2.GROUP.LABEL, levels = c("NT-Alpha", "NS-Gamma", "NS-Delta", "NS-Beta", "NS-Alpha"))) %>% arrange(TAX2.GROUP.LABEL, Z.MEDIAN) %>% mutate(TAX10 = factor(TAX10, levels = TAX10)) %>% mutate(SIG = ifelse(GROUP >0, "Significant", "Not Significant")) %>% mutate(SIG = factor(SIG, levels = c("Significant", "Not Significant")))

TAX10.titan.response.aoa.ammpost.tbl <- as_tibble(TAX10.aoa.trim.titan.ammpost$sppmax, rownames = NA) %>% rownames_to_column(var = "TAX10") %>% select (TAX10, zenv.cp, '5%', '95%', z.median, filter, maxgrp) %>% rename(ZENV.CP = "zenv.cp", LCI = "5%", UCI = "95%", Z.MEDIAN = "z.median", GROUP = "filter", MAXGRP = "maxgrp") %>% mutate(Z.MEDIAN = ifelse(MAXGRP == 1, -1 * Z.MEDIAN, Z.MEDIAN)) %>% inner_join(.,TAX10.group.titan.tbl, by = "TAX10") %>% mutate(TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NT-Alpha", "NT-Alpha", NA), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Gamma", "NS-Gamma", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Delta", "NS-Delta", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Beta", "NS-Beta", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Alpha", "NS-Alpha", TAX2.GROUP.LABEL)) %>% mutate(TAX2.GROUP.LABEL = factor(TAX2.GROUP.LABEL, levels = c("NT-Alpha", "NS-Gamma", "NS-Delta", "NS-Beta", "NS-Alpha"))) %>% arrange(TAX2.GROUP.LABEL, Z.MEDIAN) %>% mutate(TAX10 = factor(TAX10, levels = TAX10)) %>% mutate(SIG = ifelse(GROUP >0, "Significant", "Not Significant")) %>% mutate(SIG = factor(SIG, levels = c("Significant", "Not Significant")))

TAX10.titan.response.aoa.nmin.plot <- ggplot()+geom_bar(data=TAX10.titan.response.aoa.nmin.tbl, aes(x=TAX10, y=Z.MEDIAN, fill = TAX2.GROUP.LABEL, alpha = SIG), stat = "identity", width = 0.75) + scale_fill_viridis_d(name = NULL, option = "plasma", begin =0.1, end = 0.9) + scale_alpha_manual(name = NULL, values = c(1, 0.3)) + geom_hline(yintercept = 0, linetype = 2, size = 0.5, colour = "black") + labs(x="ASV Taxon", y = "Median Z-score\n(response to inorganic N availabilty)") + coord_flip() + scale_x_discrete(limits=rev(levels(TAX10.titan.response.aoa.nmin.tbl$TAX10))) + theme(axis.line = element_line(colour = "black", size = 0.5), panel.background = element_rect(fill = NA), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title = element_text(size = 12, colour = "black", hjust = -0.5),axis.title.x = element_text(size = 7, colour = "black"), axis.title.y = element_text(size = 7, colour = "black"), axis.text.x = element_text(size = 6, colour = "black"), axis.text.y = element_text(size = 6, colour = "black"), legend.title = element_text(size = 7, colour = "black"), legend.text = element_text(size = 6, colour = "black"), legend.key = element_blank())

TAX10.titan.response.aoa.ph.plot <- ggplot()+geom_bar(data=TAX10.titan.response.aoa.ph.tbl, aes(x=TAX10, y=Z.MEDIAN, fill = TAX2.GROUP.LABEL, alpha = SIG), stat = "identity", width = 0.75) + scale_fill_viridis_d(name = NULL, option = "plasma", begin =0.1, end = 0.9) + scale_alpha_manual(name = NULL, values = c(1, 0.3)) + geom_hline(yintercept = 0, linetype = 2, size = 0.5, colour = "black") + labs(x="ASV Taxon", y = "Median Z-score\n(response to inorganic N availabilty)") + coord_flip() + scale_x_discrete(limits=rev(levels(TAX10.titan.response.aoa.ph.tbl$TAX10))) + theme(axis.line = element_line(colour = "black", size = 0.5), panel.background = element_rect(fill = NA), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title = element_text(size = 12, colour = "black", hjust = -0.5),axis.title.x = element_text(size = 7, colour = "black"), axis.title.y = element_text(size = 7, colour = "black"), axis.text.x = element_text(size = 6, colour = "black"), axis.text.y = element_text(size = 6, colour = "black"), legend.title = element_text(size = 7, colour = "black"), legend.text = element_text(size = 6, colour = "black"), legend.key = element_blank())

TAX10.titan.response.aoa.ammpost.plot <- ggplot()+geom_bar(data=TAX10.titan.response.aoa.ammpost.tbl, aes(x=TAX10, y=Z.MEDIAN, fill = TAX2.GROUP.LABEL, alpha = SIG), stat = "identity", width = 0.75) + scale_fill_viridis_d(name = NULL, option = "plasma", begin =0.1, end = 0.9) + scale_alpha_manual(name = NULL, values = c(1, 0.3)) + geom_hline(yintercept = 0, linetype = 2, size = 0.5, colour = "black") + labs(x="ASV Taxon", y = "Median Z-score\n(response to inorganic N availabilty)") + coord_flip() + scale_x_discrete(limits=rev(levels(TAX10.titan.response.aoa.ammpost.tbl$TAX10))) + theme(axis.line = element_line(colour = "black", size = 0.5), panel.background = element_rect(fill = NA), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title = element_text(size = 12, colour = "black", hjust = -0.5),axis.title.x = element_text(size = 7, colour = "black"), axis.title.y = element_text(size = 7, colour = "black"), axis.text.x = element_text(size = 6, colour = "black"), axis.text.y = element_text(size = 6, colour = "black"), legend.title = element_text(size = 7, colour = "black"), legend.text = element_text(size = 6, colour = "black"), legend.key = element_blank())

TAX10.fig.input.aoa.nmin.tbl <- as_tibble(TAX10.aoa.trim.titan.nmin$sppmax, rownames = "TAX10") %>% inner_join(., TAX10.group.titan.aoa.tbl, by = "TAX10") %>% dplyr::rename(RESPONSE.GROUP = "filter")  %>% select(TAX10, TAX2.GROUP, RESPONSE.GROUP) %>% distinct()  %>% group_by(TAX2.GROUP)%>% summarise(POSITIVE.RESPONSE = length(TAX2.GROUP[RESPONSE.GROUP == 2]), NEGATIVE.RESPONSE = length(TAX2.GROUP[RESPONSE.GROUP ==1]), NO.RESPONSE = length(TAX2.GROUP[RESPONSE.GROUP ==0]), TOTAL.TAX10 = length(TAX2.GROUP)) %>% ungroup() %>% mutate(TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NT-Alpha", "NT-Alpha", NA), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Gamma", "NS-Gamma", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Delta", "NS-Delta", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Beta", "NS-Beta", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Alpha", "NS-Alpha", TAX2.GROUP.LABEL)) %>% mutate(TAX2.GROUP.LABEL = factor(TAX2.GROUP.LABEL, levels = c("NT-Alpha", "NS-Gamma", "NS-Delta", "NS-Beta", "NS-Alpha"))) %>% pivot_longer(cols=POSITIVE.RESPONSE : NEGATIVE.RESPONSE, names_to = "RESPONSE", values_to = "COUNT") %>% mutate(PROPORTION = COUNT/TOTAL.TAX10) %>% mutate(PROPORTION = ifelse(RESPONSE == "NEGATIVE.RESPONSE", -1 *PROPORTION, PROPORTION)) %>% mutate(SIG.LABEL = paste("(", COUNT, "/", TOTAL.TAX10, ")", sep="")) %>% mutate(LABEL.POSITION = ifelse(RESPONSE =="POSITIVE.RESPONSE", 0.2 + PROPORTION, PROPORTION - 0.2))

TAX10.fig.input.aoa.ph.tbl <- as_tibble(TAX10.aoa.trim.titan.ph$sppmax, rownames = "TAX10") %>% inner_join(., TAX10.group.titan.aoa.tbl, by = "TAX10") %>% dplyr::rename(RESPONSE.GROUP = "filter")  %>% select(TAX10, TAX2.GROUP, RESPONSE.GROUP) %>% distinct()  %>% group_by(TAX2.GROUP)%>% summarise(POSITIVE.RESPONSE = length(TAX2.GROUP[RESPONSE.GROUP == 2]), NEGATIVE.RESPONSE = length(TAX2.GROUP[RESPONSE.GROUP ==1]), NO.RESPONSE = length(TAX2.GROUP[RESPONSE.GROUP ==0]), TOTAL.TAX10 = length(TAX2.GROUP)) %>% ungroup() %>% mutate(TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NT-Alpha", "NT-Alpha", NA), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Gamma", "NS-Gamma", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Delta", "NS-Delta", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Beta", "NS-Beta", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Alpha", "NS-Alpha", TAX2.GROUP.LABEL)) %>% mutate(TAX2.GROUP.LABEL = factor(TAX2.GROUP.LABEL, levels = c("NT-Alpha", "NS-Gamma", "NS-Delta", "NS-Beta", "NS-Alpha"))) %>% pivot_longer(cols=POSITIVE.RESPONSE : NEGATIVE.RESPONSE, names_to = "RESPONSE", values_to = "COUNT") %>% mutate(PROPORTION = COUNT/TOTAL.TAX10) %>% mutate(PROPORTION = ifelse(RESPONSE == "NEGATIVE.RESPONSE", -1 *PROPORTION, PROPORTION)) %>% mutate(SIG.LABEL = paste("(", COUNT, "/", TOTAL.TAX10, ")", sep="")) %>% mutate(LABEL.POSITION = ifelse(RESPONSE =="POSITIVE.RESPONSE", 0.2 + PROPORTION, PROPORTION - 0.2))

TAX10.fig.input.aoa.ammpost.tbl <- as_tibble(TAX10.aoa.trim.titan.ammpost$sppmax, rownames = "TAX10") %>% inner_join(., TAX10.group.titan.aoa.tbl, by = "TAX10") %>% dplyr::rename(RESPONSE.GROUP = "filter")  %>% select(TAX10, TAX2.GROUP, RESPONSE.GROUP) %>% distinct()  %>% group_by(TAX2.GROUP)%>% summarise(POSITIVE.RESPONSE = length(TAX2.GROUP[RESPONSE.GROUP == 2]), NEGATIVE.RESPONSE = length(TAX2.GROUP[RESPONSE.GROUP ==1]), NO.RESPONSE = length(TAX2.GROUP[RESPONSE.GROUP ==0]), TOTAL.TAX10 = length(TAX2.GROUP)) %>% ungroup() %>% mutate(TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NT-Alpha", "NT-Alpha", NA), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Gamma", "NS-Gamma", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Delta", "NS-Delta", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Beta", "NS-Beta", TAX2.GROUP.LABEL), TAX2.GROUP.LABEL = ifelse(TAX2.GROUP == "NS-Alpha", "NS-Alpha", TAX2.GROUP.LABEL)) %>% mutate(TAX2.GROUP.LABEL = factor(TAX2.GROUP.LABEL, levels = c("NT-Alpha", "NS-Gamma", "NS-Delta", "NS-Beta", "NS-Alpha"))) %>% pivot_longer(cols=POSITIVE.RESPONSE : NEGATIVE.RESPONSE, names_to = "RESPONSE", values_to = "COUNT") %>% mutate(PROPORTION = COUNT/TOTAL.TAX10) %>% mutate(PROPORTION = ifelse(RESPONSE == "NEGATIVE.RESPONSE", -1 *PROPORTION, PROPORTION)) %>% mutate(SIG.LABEL = paste("(", COUNT, "/", TOTAL.TAX10, ")", sep="")) %>% mutate(LABEL.POSITION = ifelse(RESPONSE =="POSITIVE.RESPONSE", 0.2 + PROPORTION, PROPORTION - 0.2))

TAX10.titan.aoa.nmin.plot <- ggplot() +   geom_bar(data = TAX10.fig.input.aoa.nmin.tbl, aes(x = TAX2.GROUP.LABEL, y = PROPORTION, fill = RESPONSE),colour = NA, stat = "identity", alpha = 0.75) +   geom_hline(yintercept = 0,linetype = 2,size = 0.5,colour = "black") + geom_text(data = TAX10.fig.input.aoa.nmin.tbl, aes(x = TAX2.GROUP.LABEL, y = LABEL.POSITION, colour = RESPONSE), label = TAX10.fig.input.aoa.nmin.tbl$SIG.LABEL, size = 2.5) + scale_colour_viridis_d(name = "Association with\ninorganic N", labels = c("Negative", "Positive"), begin = 0.2, end = 0.8, option = "plasma", guide = FALSE) + scale_fill_viridis_d(name = "Association with\ninorgaic N", labels = c("Negative", "Positive"), begin = 0.2, end = 0.8, option = "plasma") +   scale_y_continuous(limits = c(-1.2, 1.2), expand = expansion(mult = c(0.1, 0.1)), breaks = c(-1, -0.5, 0, 0.5, 1),labels = c("1.0", "0.5", "0", "0.5", "1.0")) + coord_flip() + scale_x_discrete(limits = rev(levels(TAX10.fig.input.aoa.nmin.tbl$TAX2.GROUP.LABEL))) + labs(x = NULL, y = "Proportion of ASV in Clade") + theme(axis.line = element_line(colour = "black", size = 0.5),panel.background = element_rect(fill = NA), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title = element_text(size = 12, colour = "black", hjust = -0.35),axis.title.x = element_text(size = 7, colour = "black"), axis.title.y = element_text(size = 7, colour = "black"), axis.text.x = element_text(size = 6, colour = "black"), axis.text.y = element_text(size = 6, colour = "black"), legend.title = element_text(size = 7, colour = "black"), legend.text = element_text(size = 6, colour = "black"), legend.key = element_blank())




#Dataframe and adding plot and stand for TITAN adding up
merged.clusterized.aa.seqtab.03.t <- as.data.frame(seqtab.multisamp.nochim %>% t %>% rowsum(clusterized.aa$cluster_0_03) %>% t)
merged.clusterized.aa.seqtab.03.t$PLOT <- sample_data(ps)$PLOT
merged.clusterized.aa.seqtab.03.t$STAND <- sample_data(ps)$STAND

#Trying to manipulate the tables I have to match to what I need for TITAN. I actually already have them all summed by cluster so some of the code for AOA (or ITS) is unnecessary but I do need to get % and all that

#TAX10.all.perc.tbl <- TAX10.all.tbl  %>% group_by(TAX10) %>% summarise(TAX10.COUNT = sum(COUNT)) %>% ungroup() %>% mutate(TOTAL.COUNT = sum(TAX10.COUNT)) %>% mutate(TAX10.PERC = 100 *(TAX10.COUNT/TOTAL.COUNT))%>% select(TAX10, TAX10.PERC)

#sum(merged.clusterized.aa.seqtab.03.t[1:25]) ##gives total read sum
colSums(merged.clusterized.aa.seqtab.03.t[1:25])

cluster.03.perc <- 100 * colSums(merged.clusterized.aa.seqtab.03.t[1:25])/sum(merged.clusterized.aa.seqtab.03.t[1:25])

ps.ASV.tbl <- as.data.frame(otu_table(ps)) %>% as_tibble(., rownames="SAMPLE.ID")
ps.sample.tbl <- as.data.frame(sample_data(ps)) %>% as_tibble(.)
ps.ASV.sample.tbl <- ps.ASV.tbl %>% pivot_longer(-SAMPLE.ID, names_to = "ASV", values_to = "COUNT") %>% inner_join(ps.sample.tbl, ., by="SAMPLE.ID") %>% filter(TYPE == "environmental") %>% group_by(PLOT, STAND, ASV) %>% summarise(COUNT=sum(COUNT)) %>% ungroup() %>% group_by(ASV) %>% mutate(ASV.TOTAL = sum(COUNT)) %>% ungroup() %>% filter(ASV.TOTAL >0) %>% select(-ASV.TOTAL)

> clusterized.aa.pared <- clusterized.aa[2:11]

#IF I want to rename the the “taxa” clusters so they display what level clustering they are.
paste0("Level1_", clusterized.aa.pared$cluster_0_01)

ps.tax.tbl <- as.data.frame(clusterized.aa.pared) %>% as_tibble(., rownames="ASV") %>% rename_all(toupper) %>% pivot_longer(-ASV, names_to = "RANKING", values_to="TAXON") %>% pivot_wider(id_cols = ASV, names_from="RANKING", values_from = TAXON)

#ASV naming issues - not sure why it worked before but I am editing this in the future ~wooo
clusterized.aa.pared$ASV <- taxa_names(phyloseq.multsamp.aob)

Cluster003.all.tbl <- ps.ASV.sample.tbl %>% inner_join(ps.tax.tbl, ., by="ASV") %>% group_by(CLUSTER_0_03, STAND, PLOT) %>% summarise(COUNT = sum(COUNT)) %>% ungroup() %>% select(CLUSTER_0_03, STAND, PLOT, COUNT)

#Repeated for 0.03 and 0.01. Redo what I did before but with this new tbl to make it the same code from AOA.

Cluster003.all.perc.tbl <- Cluster003.all.tbl  %>% group_by(CLUSTER_0_03) %>% summarise(CLUSTER003.COUNT = sum(COUNT)) %>% ungroup() %>% mutate(TOTAL.COUNT = sum(CLUSTER003.COUNT)) %>% mutate(CLUSTER003.PERC = 100 *(CLUSTER003.COUNT/TOTAL.COUNT))%>% select(CLUSTER_0_03, CLUSTER003.PERC)
