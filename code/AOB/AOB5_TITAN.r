library(tidyverse)
library(TITAN2)
library(ggplot2)
library(virdisLite)

phy.aob.pruned <- readRDS("Phyloseq_AOB_Pruned.rds")
ASV.aob.tbl <- readRDS("ASV_AOB_Table.rds")
CL.all.aob.tbl <- readRDS("CL_ALL_AOB_Table.rds")

#Titan
CL01.tbl <- CL.all.aob.tbl %>% group_by(CL01, STAND, PLOT) %>% summarise(COUNT = sum(COUNT)) %>% ungroup() %>% select(CL01, STAND, PLOT, COUNT)
CL1.tbl  <- CL.all.aob.tbl %>% group_by(CL1, STAND, PLOT) %>% summarise(COUNT = sum(COUNT)) %>% ungroup() %>% select(CL1, STAND, PLOT, COUNT)

CL01.perc.tbl <- CL01.tbl  %>% group_by(CL01) %>% summarise(CL01.COUNT = sum(COUNT)) %>% ungroup() %>% mutate(TOTAL.COUNT = sum(CL01.COUNT)) %>% mutate(CL01.PERC = 100 *(CL01.COUNT/TOTAL.COUNT))%>% select(CL01, CL01.PERC)
CL01.aob.tmp.tbl <- CL01.tbl  %>% group_by(CL01) %>% mutate(CL01.TOTAL = sum(COUNT))%>% ungroup() %>% filter(CL01.TOTAL > 0) %>% group_by(PLOT) %>% mutate(PLOT.COUNT = sum(COUNT))%>% ungroup() %>% group_by(CL01, PLOT) %>% mutate(PROP = sum(COUNT)/PLOT.COUNT) %>% ungroup() %>% mutate(HELLINGER = sqrt(PROP))
CL01.aob.trim.tbl <- CL01.tbl %>% mutate(PA = ifelse(COUNT > 0,1,0))%>% group_by(CL01) %>% inner_join(., CL01.perc.tbl, by = "CL01") %>% mutate(OCCURRENCE = sum(PA))%>% ungroup() %>% filter(CL01.PERC >= 0.1 & OCCURRENCE >= 3)
CL01.CL1.tbl <- CL.all.aob.tbl %>% group_by(CL01,CL1, STAND, PLOT) %>% summarize(COUNT = sum(COUNT))%>% ungroup()%>% select(CL01, CL1, STAND, PLOT, COUNT)

#Dropped Cluster 1-6 (CL 0.1, number 6)
CL01.group.titan.aob.tbl <- CL01.CL1.tbl %>% mutate(CL1.GROUP = ifelse(CL1 == "1", "Cluster 1-1", NA)) %>% mutate(CL1.GROUP = ifelse(CL1 == "2", "Cluster 1-2", CL1.GROUP)) %>% mutate(CL1.GROUP = ifelse(CL1 == "3", "Cluster 1-3", CL1.GROUP)) %>% mutate(CL1.GROUP = ifelse(CL1 == "4", "Cluster 1-4", CL1.GROUP)) %>% mutate(CL1.GROUP = ifelse(CL1 == "5", "Cluster 1-5", CL1.GROUP)) %>% select(CL01, CL1.GROUP) %>% distinct()

phy.aob.samdata.df <- data.frame(sample_data(phy.aob.pruned)) %>%  as_tibble(., rownames = "SAMPLE.ID")
env.nmin.aob.df <- phy.aob.samdata.df %>% select(PLOT, N.MIN) %>% arrange(PLOT) %>% column_to_rownames(var = "PLOT") %>% as.data.frame(.)
env.ph.aob.df <- phy.aob.samdata.df %>% select(PLOT, SOIL.PH) %>% arrange(PLOT) %>% column_to_rownames(var = "PLOT") %>% as.data.frame(.)
env.ammpost.aob.df <- phy.aob.samdata.df %>% select(PLOT, AMM.CORRECTED.POST) %>% arrange(PLOT) %>% column_to_rownames(var = "PLOT") %>% as.data.frame(.)

phy.aob.hel <- decostand(phyloseq::otu_table(phy.aob.pruned), method = "hellinger")

CL01.aob.trim.titan.in.df <- CL01.aob.trim.tbl %>% select(CL01) %>% distinct() %>% inner_join(., CL01.aob.tmp.tbl, by = "CL01") %>% select(PLOT, CL01, HELLINGER) %>% pivot_wider(id_cols = PLOT, names_from = "CL01", values_from = "HELLINGER") %>% arrange(PLOT) %>% column_to_rownames(var = "PLOT") %>% as.data.frame(.)
CL01.aob.trim.titan.nmin <- titan(env.nmin.aob.df, CL01.aob.trim.titan.in.df, minSplt =5, numPerm = 1000, boot = TRUE, nBoot = 1000, imax = FALSE, ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95, ncpus =4, memory = FALSE)
CL01.aob.trim.titan.ph <- titan(env.ph.aob.df, CL01.aob.trim.titan.in.df, minSplt =5, numPerm = 1000, boot = TRUE, nBoot = 1000, imax = FALSE, ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95, ncpus =4, memory = FALSE)
CL01.aob.trim.titan.ammpost <- titan(env.ammpost.aob.df, CL01.aob.trim.titan.in.df, minSplt =5, numPerm = 1000, boot = TRUE, nBoot = 1000, imax = FALSE, ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95, ncpus =4, memory = FALSE)
CL01.group.titan.aob.tbl$CL01 <- as.character(CL01.group.titan.aob.tbl$CL01)
CL01.titan.response.aob.nmin.tbl <- as_tibble(CL01.aob.trim.titan.nmin$sppmax, rownames = NA) %>% rownames_to_column(var = "CL01") %>% select (CL01, zenv.cp, '5%', '95%', z.median, filter, maxgrp) %>% rename(ZENV.CP = "zenv.cp", LCI = "5%", UCI = "95%", Z.MEDIAN = "z.median", GROUP = "filter", MAXGRP = "maxgrp") %>% mutate(Z.MEDIAN = ifelse(MAXGRP == 1, -1 * Z.MEDIAN, Z.MEDIAN)) %>% inner_join(.,CL01.group.titan.aob.tbl, by = "CL01") %>% mutate(CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-1", "Cluster 1-1", NA), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-2", "Cluster 1-2", CL1.GROUP.LABEL), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-3", "Cluster 1-3", CL1.GROUP.LABEL), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-4", "Cluster 1-4", CL1.GROUP.LABEL), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-5", "Cluster 1-5", CL1.GROUP.LABEL)) %>% mutate(CL1.GROUP.LABEL = factor(CL1.GROUP.LABEL, levels = c("Cluster 1-1", "Cluster 1-2", "Cluster 1-3", "Cluster 1-4", "Cluster 1-5"))) %>% arrange(CL1.GROUP.LABEL, Z.MEDIAN) %>% mutate(CL01 = factor(CL01, levels = unique(CL01))) %>% mutate(SIG = ifelse(GROUP >0, "Significant", "Not Significant")) %>% mutate(SIG = factor(SIG, levels = c("Significant", "Not Significant")))
CL01.titan.response.aob.ph.tbl <- as_tibble(CL01.aob.trim.titan.ph$sppmax, rownames = NA) %>% rownames_to_column(var = "CL01") %>% select (CL01, zenv.cp, '5%', '95%', z.median, filter, maxgrp) %>% rename(ZENV.CP = "zenv.cp", LCI = "5%", UCI = "95%", Z.MEDIAN = "z.median", GROUP = "filter", MAXGRP = "maxgrp") %>% mutate(Z.MEDIAN = ifelse(MAXGRP == 1, -1 * Z.MEDIAN, Z.MEDIAN)) %>% inner_join(.,CL01.group.titan.aob.tbl, by = "CL01") %>% mutate(CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-1", "Cluster 1-1", NA), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-2", "Cluster 1-2", CL1.GROUP.LABEL), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-3", "Cluster 1-3", CL1.GROUP.LABEL), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-4", "Cluster 1-4", CL1.GROUP.LABEL), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-5", "Cluster 1-5", CL1.GROUP.LABEL)) %>% mutate(CL1.GROUP.LABEL = factor(CL1.GROUP.LABEL, levels = c("Cluster 1-1", "Cluster 1-2", "Cluster 1-3", "Cluster 1-4", "Cluster 1-5"))) %>% arrange(CL1.GROUP.LABEL, Z.MEDIAN) %>% mutate(CL01 = factor(CL01, levels = unique(CL01))) %>% mutate(SIG = ifelse(GROUP >0, "Significant", "Not Significant")) %>% mutate(SIG = factor(SIG, levels = c("Significant", "Not Significant")))
CL01.titan.response.aob.ammpost.tbl <- as_tibble(CL01.aob.trim.titan.ammpost$sppmax, rownames = NA) %>% rownames_to_column(var = "CL01") %>% select (CL01, zenv.cp, '5%', '95%', z.median, filter, maxgrp) %>% rename(ZENV.CP = "zenv.cp", LCI = "5%", UCI = "95%", Z.MEDIAN = "z.median", GROUP = "filter", MAXGRP = "maxgrp") %>% mutate(Z.MEDIAN = ifelse(MAXGRP == 1, -1 * Z.MEDIAN, Z.MEDIAN)) %>% inner_join(.,CL01.group.titan.aob.tbl, by = "CL01") %>% mutate(CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-1", "Cluster 1-1", NA), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-2", "Cluster 1-2", CL1.GROUP.LABEL), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-3", "Cluster 1-3", CL1.GROUP.LABEL), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-4", "Cluster 1-4", CL1.GROUP.LABEL), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-5", "Cluster 1-5", CL1.GROUP.LABEL)) %>% mutate(CL1.GROUP.LABEL = factor(CL1.GROUP.LABEL, levels = c("Cluster 1-1", "Cluster 1-2", "Cluster 1-3", "Cluster 1-4", "Cluster 1-5"))) %>% arrange(CL1.GROUP.LABEL, Z.MEDIAN) %>% mutate(CL01 = factor(CL01, levels = unique(CL01))) %>% mutate(SIG = ifelse(GROUP >0, "Significant", "Not Significant")) %>% mutate(SIG = factor(SIG, levels = c("Significant", "Not Significant")))

CL01.titan.response.aob.nmin.plot <- ggplot()+geom_bar(data=CL01.titan.response.aob.nmin.tbl, aes(x=CL01, y=Z.MEDIAN, fill = CL1.GROUP.LABEL, alpha = SIG), stat = "identity", width = 0.75) + scale_fill_viridis_d(name = NULL, option = "plasma", begin =0.1, end = 0.9) + scale_alpha_manual(name = NULL, values = c(1, 0.3)) + geom_hline(yintercept = 0, linetype = 2, size = 0.5, colour = "black") + labs(x="ASV Taxon", y = "Median Z-score\n(response to inorganic N availabilty)") + coord_flip() + scale_x_discrete(limits=rev(levels(CL01.titan.response.aob.nmin.tbl$CL01))) + theme(axis.line = element_line(colour = "black", size = 0.5), panel.background = element_rect(fill = NA), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title = element_text(size = 12, colour = "black", hjust = -0.5),axis.title.x = element_text(size = 7, colour = "black"), axis.title.y = element_text(size = 7, colour = "black"), axis.text.x = element_text(size = 6, colour = "black"), axis.text.y = element_text(size = 6, colour = "black"), legend.title = element_text(size = 7, colour = "black"), legend.text = element_text(size = 6, colour = "black"), legend.key = element_blank())
CL01.titan.response.aob.ph.plot <- ggplot()+geom_bar(data=CL01.titan.response.aob.ph.tbl, aes(x=CL01, y=Z.MEDIAN, fill = CL1.GROUP.LABEL, alpha = SIG), stat = "identity", width = 0.75) + scale_fill_viridis_d(name = NULL, option = "plasma", begin =0.1, end = 0.9) + scale_alpha_manual(name = NULL, values = c(1, 0.3)) + geom_hline(yintercept = 0, linetype = 2, size = 0.5, colour = "black") + labs(x="ASV Taxon", y = "Median Z-score\n(response to Soil pH)") + coord_flip() + scale_x_discrete(limits=rev(levels(CL01.titan.response.aob.ph.tbl$CL01))) + theme(axis.line = element_line(colour = "black", size = 0.5), panel.background = element_rect(fill = NA), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title = element_text(size = 12, colour = "black", hjust = -0.5),axis.title.x = element_text(size = 7, colour = "black"), axis.title.y = element_text(size = 7, colour = "black"), axis.text.x = element_text(size = 6, colour = "black"), axis.text.y = element_text(size = 6, colour = "black"), legend.title = element_text(size = 7, colour = "black"), legend.text = element_text(size = 6, colour = "black"), legend.key = element_blank())
CL01.titan.response.aob.ammpost.plot <- ggplot()+geom_bar(data=CL01.titan.response.aob.ammpost.tbl, aes(x=CL01, y=Z.MEDIAN, fill = CL1.GROUP.LABEL, alpha = SIG), stat = "identity", width = 0.75) + scale_fill_viridis_d(name = NULL, option = "plasma", begin =0.1, end = 0.9) + scale_alpha_manual(name = NULL, values = c(1, 0.3)) + geom_hline(yintercept = 0, linetype = 2, size = 0.5, colour = "black") + labs(x="ASV Taxon", y = "Median Z-score\n(response to Soil pH)") + coord_flip() + scale_x_discrete(limits=rev(levels(CL01.titan.response.aob.ammpost.tbl$CL01))) + theme(axis.line = element_line(colour = "black", size = 0.5), panel.background = element_rect(fill = NA), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title = element_text(size = 12, colour = "black", hjust = -0.5),axis.title.x = element_text(size = 7, colour = "black"), axis.title.y = element_text(size = 7, colour = "black"), axis.text.x = element_text(size = 6, colour = "black"), axis.text.y = element_text(size = 6, colour = "black"), legend.title = element_text(size = 7, colour = "black"), legend.text = element_text(size = 6, colour = "black"), legend.key = element_blank())

CL01.fig.input.aob.nmin.tbl <- as_tibble(CL01.aob.trim.titan.nmin$sppmax, rownames = "CL01") %>% inner_join(., CL01.group.titan.aob.tbl, by = "CL01") %>% dplyr::rename(RESPONSE.GROUP = "filter")  %>% select(CL01, CL1.GROUP, RESPONSE.GROUP) %>% distinct()  %>% group_by(CL1.GROUP)%>% summarise(POSITIVE.RESPONSE = length(CL1.GROUP[RESPONSE.GROUP == 2]), NEGATIVE.RESPONSE = length(CL1.GROUP[RESPONSE.GROUP ==1]), NO.RESPONSE = length(CL1.GROUP[RESPONSE.GROUP ==0]), TOTAL.CL01 = length(CL1.GROUP)) %>% ungroup() %>% mutate(CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-1", "Cluster 1-1", NA), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-2", "Cluster 1-2", CL1.GROUP.LABEL), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-3", "Cluster 1-3", CL1.GROUP.LABEL), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-4", "Cluster 1-4", CL1.GROUP.LABEL), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-5", "Cluster 1-5", CL1.GROUP.LABEL)) %>% mutate(CL1.GROUP.LABEL = factor(CL1.GROUP.LABEL, levels = c("Cluster 1-1", "Cluster 1-2", "Cluster 1-3", "Cluster 1-4", "Cluster 1-5"))) %>% pivot_longer(cols=POSITIVE.RESPONSE : NEGATIVE.RESPONSE, names_to = "RESPONSE", values_to = "COUNT") %>% mutate(PROPORTION = COUNT/TOTAL.CL01) %>% mutate(PROPORTION = ifelse(RESPONSE == "NEGATIVE.RESPONSE", -1 *PROPORTION, PROPORTION)) %>% mutate(SIG.LABEL = paste("(", COUNT, "/", TOTAL.CL01, ")", sep="")) %>% mutate(LABEL.POSITION = ifelse(RESPONSE =="POSITIVE.RESPONSE", 0.2 + PROPORTION, PROPORTION - 0.2))
CL01.fig.input.aob.ph.tbl <- as_tibble(CL01.aob.trim.titan.ph$sppmax, rownames = "CL01") %>% inner_join(., CL01.group.titan.aob.tbl, by = "CL01") %>% dplyr::rename(RESPONSE.GROUP = "filter")  %>% select(CL01, CL1.GROUP, RESPONSE.GROUP) %>% distinct()  %>% group_by(CL1.GROUP)%>% summarise(POSITIVE.RESPONSE = length(CL1.GROUP[RESPONSE.GROUP == 2]), NEGATIVE.RESPONSE = length(CL1.GROUP[RESPONSE.GROUP ==1]), NO.RESPONSE = length(CL1.GROUP[RESPONSE.GROUP ==0]), TOTAL.CL01 = length(CL1.GROUP)) %>% ungroup() %>% mutate(CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-1", "Cluster 1-1", NA), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-2", "Cluster 1-2", CL1.GROUP.LABEL), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-3", "Cluster 1-3", CL1.GROUP.LABEL), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-4", "Cluster 1-4", CL1.GROUP.LABEL), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-5", "Cluster 1-5", CL1.GROUP.LABEL)) %>% mutate(CL1.GROUP.LABEL = factor(CL1.GROUP.LABEL, levels = c("Cluster 1-1", "Cluster 1-2", "Cluster 1-3", "Cluster 1-4", "Cluster 1-5"))) %>% pivot_longer(cols=POSITIVE.RESPONSE : NEGATIVE.RESPONSE, names_to = "RESPONSE", values_to = "COUNT") %>% mutate(PROPORTION = COUNT/TOTAL.CL01) %>% mutate(PROPORTION = ifelse(RESPONSE == "NEGATIVE.RESPONSE", -1 *PROPORTION, PROPORTION)) %>% mutate(SIG.LABEL = paste("(", COUNT, "/", TOTAL.CL01, ")", sep="")) %>% mutate(LABEL.POSITION = ifelse(RESPONSE =="POSITIVE.RESPONSE", 0.2 + PROPORTION, PROPORTION - 0.2))
CL01.fig.input.aob.ammpost.tbl <- as_tibble(CL01.aob.trim.titan.ammpost$sppmax, rownames = "CL01") %>% inner_join(., CL01.group.titan.aob.tbl, by = "CL01") %>% dplyr::rename(RESPONSE.GROUP = "filter")  %>% select(CL01, CL1.GROUP, RESPONSE.GROUP) %>% distinct()  %>% group_by(CL1.GROUP)%>% summarise(POSITIVE.RESPONSE = length(CL1.GROUP[RESPONSE.GROUP == 2]), NEGATIVE.RESPONSE = length(CL1.GROUP[RESPONSE.GROUP ==1]), NO.RESPONSE = length(CL1.GROUP[RESPONSE.GROUP ==0]), TOTAL.CL01 = length(CL1.GROUP)) %>% ungroup() %>% mutate(CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-1", "Cluster 1-1", NA), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-2", "Cluster 1-2", CL1.GROUP.LABEL), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-3", "Cluster 1-3", CL1.GROUP.LABEL), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-4", "Cluster 1-4", CL1.GROUP.LABEL), CL1.GROUP.LABEL = ifelse(CL1.GROUP == "Cluster 1-5", "Cluster 1-5", CL1.GROUP.LABEL)) %>% mutate(CL1.GROUP.LABEL = factor(CL1.GROUP.LABEL, levels = c("Cluster 1-1", "Cluster 1-2", "Cluster 1-3", "Cluster 1-4", "Cluster 1-5"))) %>% pivot_longer(cols=POSITIVE.RESPONSE : NEGATIVE.RESPONSE, names_to = "RESPONSE", values_to = "COUNT") %>% mutate(PROPORTION = COUNT/TOTAL.CL01) %>% mutate(PROPORTION = ifelse(RESPONSE == "NEGATIVE.RESPONSE", -1 *PROPORTION, PROPORTION)) %>% mutate(SIG.LABEL = paste("(", COUNT, "/", TOTAL.CL01, ")", sep="")) %>% mutate(LABEL.POSITION = ifelse(RESPONSE =="POSITIVE.RESPONSE", 0.2 + PROPORTION, PROPORTION - 0.2))

CL01.titan.aob.nmin.plot <- ggplot() +   geom_bar(data = CL01.fig.input.aob.nmin.tbl, aes(x = CL1.GROUP.LABEL, y = PROPORTION, fill = RESPONSE),colour = NA, stat = "identity", alpha = 0.75) +   geom_hline(yintercept = 0,linetype = 2,size = 0.5,colour = "black") + geom_text(data = CL01.fig.input.aob.nmin.tbl, aes(x = CL1.GROUP.LABEL, y = LABEL.POSITION, colour = RESPONSE), label = CL01.fig.input.aob.nmin.tbl$SIG.LABEL, size = 2.5) + scale_colour_viridis_d(name = "Association with\ninorganic N", labels = c("Negative", "Positive"), begin = 0.2, end = 0.8, option = "plasma", guide = FALSE) + scale_fill_viridis_d(name = "Association with\ninorgaic N", labels = c("Negative", "Positive"), begin = 0.2, end = 0.8, option = "plasma") +   scale_y_continuous(limits = c(-1.2, 1.2), expand = expansion(mult = c(0.1, 0.1)), breaks = c(-1, -0.5, 0, 0.5, 1),labels = c("1.0", "0.5", "0", "0.5", "1.0")) + coord_flip() + scale_x_discrete(limits = rev(levels(CL01.fig.input.aob.nmin.tbl$CL1.GROUP.LABEL))) + labs(x = NULL, y = "Proportion of ASV in 0.1 Cluster") + theme(axis.line = element_line(colour = "black", size = 0.5),panel.background = element_rect(fill = NA), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title = element_text(size = 12, colour = "black", hjust = -0.35),axis.title.x = element_text(size = 7, colour = "black"), axis.title.y = element_text(size = 7, colour = "black"), axis.text.x = element_text(size = 6, colour = "black"), axis.text.y = element_text(size = 6, colour = "black"), legend.title = element_text(size = 7, colour = "black"), legend.text = element_text(size = 6, colour = "black"), legend.key = element_blank())
ggsave("Titan_Nmin_AOB.pdf", width=8.5, height=10)

CL01.titan.aob.ph.plot <- ggplot() +   geom_bar(data = CL01.fig.input.aob.ph.tbl, aes(x = CL1.GROUP.LABEL, y = PROPORTION, fill = RESPONSE),colour = NA, stat = "identity", alpha = 0.75) +   geom_hline(yintercept = 0,linetype = 2,size = 0.5,colour = "black") + geom_text(data = CL01.fig.input.aob.ph.tbl, aes(x = CL1.GROUP.LABEL, y = LABEL.POSITION, colour = RESPONSE), label = CL01.fig.input.aob.ph.tbl$SIG.LABEL, size = 2.5) + scale_colour_viridis_d(name = "Association with\nSoil pH", labels = c("Negative", "Positive"), begin = 0.2, end = 0.8, option = "plasma", guide = FALSE) + scale_fill_viridis_d(name = "Association with\nSoil pH", labels = c("Negative", "Positive"), begin = 0.2, end = 0.8, option = "plasma") +   scale_y_continuous(limits = c(-1.2, 1.2), expand = expansion(mult = c(0.1, 0.1)), breaks = c(-1, -0.5, 0, 0.5, 1),labels = c("1.0", "0.5", "0", "0.5", "1.0")) + coord_flip() + scale_x_discrete(limits = rev(levels(CL01.fig.input.aob.ph.tbl$CL1.GROUP.LABEL))) + labs(x = NULL, y = "Proportion of ASV in 0.1 Cluster") + theme(axis.line = element_line(colour = "black", size = 0.5),panel.background = element_rect(fill = NA), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title = element_text(size = 12, colour = "black", hjust = -0.35),axis.title.x = element_text(size = 7, colour = "black"), axis.title.y = element_text(size = 7, colour = "black"), axis.text.x = element_text(size = 6, colour = "black"), axis.text.y = element_text(size = 6, colour = "black"), legend.title = element_text(size = 7, colour = "black"), legend.text = element_text(size = 6, colour = "black"), legend.key = element_blank())
ggsave("Titan_Soil_pH_AOB.pdf", width=8.5, height=10)

CL01.titan.aob.ammpost.plot <- ggplot() +   geom_bar(data = CL01.fig.input.aob.ammpost.tbl, aes(x = CL1.GROUP.LABEL, y = PROPORTION, fill = RESPONSE),colour = NA, stat = "identity", alpha = 0.75) +   geom_hline(yintercept = 0,linetype = 2,size = 0.5,colour = "black") + geom_text(data = CL01.fig.input.aob.ammpost.tbl, aes(x = CL1.GROUP.LABEL, y = LABEL.POSITION, colour = RESPONSE), label = CL01.fig.input.aob.ammpost.tbl$SIG.LABEL, size = 2.5) + scale_colour_viridis_d(name = "Association with\nPost Incubation Ammonium", labels = c("Negative", "Positive"), begin = 0.2, end = 0.8, option = "plasma", guide = FALSE) + scale_fill_viridis_d(name = "Association with\nPost Incubation Ammonium", labels = c("Negative", "Positive"), begin = 0.2, end = 0.8, option = "plasma") +   scale_y_continuous(limits = c(-1.2, 1.2), expand = expansion(mult = c(0.1, 0.1)), breaks = c(-1, -0.5, 0, 0.5, 1),labels = c("1.0", "0.5", "0", "0.5", "1.0")) + coord_flip() + scale_x_discrete(limits = rev(levels(CL01.fig.input.aob.ammpost.tbl$CL1.GROUP.LABEL))) + labs(x = NULL, y = "Proportion of ASV in 0.1 Cluster") + theme(axis.line = element_line(colour = "black", size = 0.5),panel.background = element_rect(fill = NA), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), plot.title = element_text(size = 12, colour = "black", hjust = -0.35),axis.title.x = element_text(size = 7, colour = "black"), axis.title.y = element_text(size = 7, colour = "black"), axis.text.x = element_text(size = 6, colour = "black"), axis.text.y = element_text(size = 6, colour = "black"), legend.title = element_text(size = 7, colour = "black"), legend.text = element_text(size = 6, colour = "black"), legend.key = element_blank())
ggsave("Titan_Amm_Post_AOB.pdf", width=8.5, height=10)