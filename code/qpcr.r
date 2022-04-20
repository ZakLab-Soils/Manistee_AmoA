
# Load libraries
library(tidyverse)
library(mgcv)
library(itsadug)
library(ggpubr)
library(viridis)

#source("code/functions.R")

#read in soil.ph and compiled data moisture and DNA extraction weights
soil.ph <- read.csv("../Manistee_spring_2019_soil_pH.txt", sep = "\t")
soil.ph$PLOT <- paste0("Plot_", soil.ph$PLOT)

# Get soil DNA extraction masses <- think I'm just going to make this a data item. Every time I know that Will areadly published this code I feel I don't need to keep soil.DNA.extraction.tbl <- read.csv("../soil_DNA_extractino_weights.csv)
 
# Read in data - maybe as aoa.plots filter
soil.compiled.data.tbl <- readRDS(file = "data/working/soil.compiled.data.tbl.RData")

#Info for qpcr analyis
aoa.copies.tbl1 <- read_tsv(file="Text Report Data_aoa_plate1.txt")
aoa.copies.tbl1 <- aoa.copies.tbl1 %>% dplyr::filter(`Well Type` == "Unknown" & Threshold != "Reference") %>% dplyr::mutate(NUM = gsub("[A-Za-z]", "", fixed = FALSE, aoa.copies.tbl1$Well)) %>% mutate(NUM = factor(NUM, levels = c("4", "5", "6", "7", "8", "9", "10", "11", "12"))) %>% arrange(NUM) %>% mutate(Quantity = as.numeric(Quantity), Ct = as.numeric(Ct), PLOT = c(1:nrow(.)), PLOT = paste("Plot", PLOT, sep = "_"), RUN = "Run_1") %>% rename(INITIAL.COPIES = "Quantity") %>% select(PLOT, RUN, INITIAL.COPIES, Ct)

aoa.copies.tbl3 <- read_tsv(file="Text Report Data_aoa_plate3.txt")
aoa.copies.tbl3 <- aoa.copies.tbl3 %>% dplyr::filter(`Well Type` == "Unknown" & Threshold != "Reference") %>% dplyr::mutate(NUM = gsub("[A-Za-z]", "", fixed = FALSE, aoa.copies.tbl3$Well)) %>% mutate(NUM = factor(NUM, levels = c("4", "5", "6", "7", "8", "9", "10", "11", "12"))) %>% arrange(NUM) %>% mutate(Quantity = as.numeric(Quantity), Ct = as.numeric(Ct), PLOT = c(1:nrow(.)), PLOT = paste("Plot", PLOT, sep = "_"), RUN = "Run_3") %>% rename(INITIAL.COPIES = "Quantity") %>% select(PLOT, RUN, INITIAL.COPIES, Ct)

#Check Ct

aoa.copies.ctcheck <- soil.compiled.data.tbl%>% select(STAND,PLOT) %>% inner_join(., aoa.copies.tbl1, by = "PLOT")
aoa.copies.ctcheck$Ct3 <- aoa.copies.tb13$Ct
aoa.copies.ctcheck$STAND <- factor(aoa.copies.ctcheck$STAND, levels = c("Stand_20", "Stand_50", "Stand_31", "Stand_3", "Stand_58", "Stand_9", "Stand_7", "Stand_6", "Stand_41", "Stand_24", "Stand_100", "Stand_22"))
aoa.copies.ctcheck$means <- rowMeans(aoa.copies.ctcheck[,c("Ct", "Ct3")])

aoa.copies.qpcrclean<- dplyr::bind_rows(aoa.copies.tbl1, aoa.copies.tbl3) %>% inner_join(., soil.compiled.data.tbl)median.ct <- aoa.copies.qpcrclean %>% group_by(STAND,PLOT) %>% summarise(media.ct = median(Ct, na.rm = TRUE), count = sum(!is.na(Ct)))

ctrun1 <- ggplot(aoa.copies.notaveraged, aes(y=Ct, x = STAND)) +geom_boxplot() + ggtitle("AOA amoaA qPCR Ct values Run 1") + theme(axis.text.x = element_text(angle = 90))
ctrun3 <- ggplot(aoa.copies.notaveraged, aes(y=Ct3, x = STAND)) +geom_boxplot() + ggtitle("AOA amoaA qPCR Ct values Run 3") + theme(axis.text.x = element_text(angle = 90))

ggplot(aoa.copies.combine.tblB, aes(y=COPY.NUMBER, x = STAND)) +geom_boxplot() + ggtitle("AOA amoaA qPCR Average Copy numbers by Stand") + theme(axis.text.x = element_text(angle = 90))

#after ct check put in command %>% filter(PLOT != "Plot_59" & PLOT != "Plot_19") 
aoa.copies.combine.tbl <- dplyr::bind_rows(aoa.copies.tbl1, aoa.copies.tbl3) %>% dplyr::group_by(PLOT) %>% dplyr::summarise(MEAN.INITIAL.COPIES = mean(INITIAL.COPIES)) %>% dplyr::ungroup() %>% dplyr::inner_join(., soil.DNA.extraction.masses.tbl, by="PLOT") %>% dplyr::mutate(VOLUME = 1, DILUTION.FACTOR =10, EXTR.VOLUME = 400) %>% dplyr::mutate(COPY.NUMBER = (MEAN.INITIAL.COPIES * VOLUME *DILUTION.FACTOR * EXTR.VOLUME * (100/150) * (1/SOIL.DRY.MASS))) %>% dplyr:: inner_join(., soil.compiled.data.tbl, by = "PLOT") %>% dplyr::inner_join(., soil.ph, by = "PLOT")
