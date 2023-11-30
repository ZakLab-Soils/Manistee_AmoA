#library(tibble)
library(tidyverse)
library(DECIPHER)
library(Biostrings)

phy.aob.pruned <- readRDS("Phyloseq_AOB_Pruned.rds")
npoc <- 4

#DECIPHER to Cluster the ASV sequences
ASV.sequences <- colnames(otu_table(phy.aob.pruned))
sample.names <- rownames(otu_table(phy.aob.pruned))
dna.clustering <- Biostrings::DNAStringSet(ASV.sequences)

#Rename reads to ASV
taxa_names(phy.aob.pruned) <- paste0("ASV_", seq(ntaxa(phy.aob.pruned)))
names(dna.clustering) <- taxa_names(phy.aob.pruned)

#Cluster and make into table
set.seed(496571)
clusters <- DECIPHER::Clusterize(dna.clustering, cutoff=seq(0, 0.1, 0.01))
colnames(clusters) <- c("Cl0", "Cl01", "Cl02", "Cl03", "Cl04", "Cl05", "Cl06", "Cl07", "Cl08", "Cl09", "Cl1")

phy.aob.asv.df <- data.frame(otu_table(phy.aob.pruned)) %>% as_tibble(., rownames="SAMPLE.ID")
phy.aob.samdata.df <- data.frame(sample_data(phy.aob.pruned)) %>%  as_tibble(., rownames = "SAMPLE.ID")
phy.aob.asv.samp.df <- phy.aob.asv.df %>% pivot_longer(-SAMPLE.ID, names_to = "ASV", values_to = "COUNT") %>% inner_join(phy.aob.samdata.df, ., by="SAMPLE.ID") %>% group_by(PLOT, STAND, ASV) %>% summarise(COUNT=sum(COUNT)) %>% ungroup() %>% group_by(ASV) %>% mutate(ASV.TOTAL = sum(COUNT)) %>% ungroup() %>% filter(ASV.TOTAL >0) %>% select(-ASV.TOTAL)
phy.aob.tax.df <- data.frame(clusters) %>% as_tibble(., rownames="ASV") %>% rename_all(toupper)  %>% pivot_longer(-ASV, names_to = "RANKING", values_to="TAXON") %>% pivot_wider(id_cols = ASV, names_from="RANKING", values_from = TAXON)

Cl.all.aob.tbl <- phy.aob.asv.samp.df %>% inner_join(phy.aob.tax.df, ., by = "ASV")

saveRDS(phy.aob.asv.samp.df, file="ASV_AOB_Table.rds")
saveRDS(Cl.all.aob.tbl, file = "CL_ALL_AOB_Table.rds")
