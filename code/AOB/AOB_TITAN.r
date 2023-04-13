library(tidyverse)
library(TITAN2)
library(ggplot2)
library(virdisLite)

phy.aoa.pruned <- readRDS("Phyloseq_AOA_Pruned.rds")
ASV.aoa.tbl <- readRDS("ASV_AOA_Table.rds")
TAX.all.aoa.tbl <- readRDS("TAX_ALL_AOA_Table.rds")



#Titan

Cl01.tbl <- phy.aob.asv.samp.df %>% inner_join(phy.aob.tax.df, ., by="ASV") %>% group_by(CL01, STAND, PLOT) %>% summarise(COUNT = sum(COUNT)) %>% ungroup() %>% select(CL01, STAND, PLOT, COUNT)
Cl03.tbl <- phy.aob.asv.samp.df %>% inner_join(phy.aob.tax.df, ., by="ASV") %>% group_by(CL03, STAND, PLOT) %>% summarise(COUNT = sum(COUNT)) %>% ungroup() %>% select(CL03, STAND, PLOT, COUNT)
Cl1.tbl <- phy.aob.asv.samp.df %>% inner_join(phy.aob.tax.df, ., by="ASV") %>% group_by(CL1, STAND, PLOT) %>% summarise(COUNT = sum(COUNT)) %>% ungroup() %>% select(CL1, STAND, PLOT, COUNT)

Cl01.perc.tbl <- Cl01.tbl  %>% group_by(CL01) %>% summarise(CL01.COUNT = sum(COUNT)) %>% ungroup() %>% mutate(TOTAL.COUNT = sum(CL01.COUNT)) %>% mutate(CL01.PERC = 100 *(CL01.COUNT/TOTAL.COUNT))%>% select(CL01, CL01.PERC)

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
