
library(dada2)
library(ShortRead)
library(Biostrings)

#setwd for folder where plan to sort and store output for R ("C:/Users/17348/Desktop/amoa-AOA/aoa/2023")

#set.seed(101279)

#source("functions.R") ##I’m trying to skip this so I need to see when it gets used. I think it’s only for STAND <-> PLOT assignment so I can create a file for that info (metadata file that is the compiled soil data (ph, aq2, etc) with assignments

path <- "../../reads/" #This folder contains the sequence data - can be same as wd if the reads are in that folder
list.files(path)
                    
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#Check quality of reads - this step takes a long time and may need to be broken down into sections.
plotQualityProfile(fnFs)
plotQualityProfile(fnRs)
#plotQualityProfile(fnFs[11:14])

##Primer removal and checking for primer contamination

fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
 
noN.ft <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN=0, multithread = FALSE)

cutadapt <- "/Users/17348/AppData/Local/Programs/Python/Python310/Scripts/cutadapt"
path.cut <- file.path(path, "cutadapt_processed")

if(!dir.exists(path.cut)) dir.create(path.cut)

fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD <- "ATGGTCTGGCTWAGACG"
REV <- "GCCATCCATCTGTATGTCCA"
REV.RC <- dada2:::rc(REV)
FWD.RC <- dada2:::rc(FWD)

R1.flags <- paste("-g", FWD, "-a", REV.RC)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

for(i in seq_along(fnFs)){system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "-o", fnFs.cut[i], "-p", fnRs.cut[i], fnFs.filtN[i], fnRs.filtN[i]))}

#Cutadapt Command line parameters: -g ATGGTCTGGCTWAGACG -a TGGACATACAGATGGATGGC -G GCCATCCATCTGTATGTCCA -A CGTCTWAGCCAGACCAT -n 2 -o data/cutadapt_processed/AOA7_S1_L001_R1_001.fastq.gz #-p data/cutadapt_processed/AOA7_S1_L001_R2_001.fastq.gz data/filtN/AOA7_S1_L001_R1_001.fastq.gz data/filtN/AOA7_S1_L001_R2_001.fastq.gz

#See how much the primers appear in the data for selected libraries. To check all libraries change[N] to different library numbers.

allOrients <- function(primer) {
require(Biostrings)
dna <- DNAString(primer)
orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
RevComp = reverseComplement(dna))
return(sapply(orients, toString))
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

primerHits <- function(primer, fn) {
nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
return(sum(nhits > 0))
}

rbind(FWD.FowardReads = sapply(FWD.orients, primerHits, fn=fnFs.filtN[[1]]), FWD.ReverseReads = sapply(FWD.orients, primerHits, fn=fnRs.filtN[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits, fn=fnFs.filtN[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn=fnRs.filtN[[1]]))

#                Forward Complement Reverse RevComp
#FWD.FowardReads        x          x       x       x
#FWD.ReverseReads       x          x       x       x
#REV.ForwardReads       x          x       x       x
#REV.ReverseReads       x          x       x       x

rbind(FWD.FowardReads = sapply(FWD.orients, primerHits, fn=fnFs.cut[[1]]), # Quantify error rateFWD.ReverseReads = sapply(FWD.orients, primerHits, fn=fnRs.cut[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits, fn=fnFs.cut[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn=fnRs.cut[[1]]))
#                Forward Complement Reverse RevComp
#FWD.FowardReads        0          0       0       0
#FWD.ReverseReads       0          0       0       0
#REV.ForwardReads       0          0       0       0
#REV.ReverseReads       0          0       0       0

cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))

## Filter and trim for length, quality and size; Calculate error rate.

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

filter.summary <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2,2), truncQ = 2, minLen = c(200,200), rm.phix = TRUE, compress = TRUE, multithread = FALSE)

#Double check the quality looks acceptable across reads. We felt it was not not necessary to trim the read lengths for this data set.
plotQualityProfile(filtFs)
plotQualityProfile(filtRs)

##Error rate processing, dereplication and ASV inference and sequence table formation

errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)

#Plot error rates - check if black line (estimated error) matches points (observed error)
errF.plot <- plotErrors(errF, nominalQ = TRUE)
errR.plot <- plotErrors(errR, nominalQ = TRUE)

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
dadaRs <- dada(derepRs, err = errR, multithread = FALSE)

#Sequences do not overlap for the most part, and we chose the reverse reads to use
seqtab <- makeSequenceTable(dadaRs)
dim(seqtab)
saveRDS(seqtab, file = "data/seqtab.rds")

#Remove ASVs from the data set that are not in at least 2 plots
ASVs.multisamp <- integer()
for (i in 1:365) {if(sum(seqtab[,i]==0) != 37){tmp.paste <- i; ASVs.multisamp <- append(ASVs.multisamp, tmp.paste) }}
seqtab.multisamp <- seqtab[, ASVs.multisamp]
dim(seqtab.multisamp)
      
seqtab.nochim <- removeBimeraDenovo(seqtab.multisamp, method = "pooled", multithread = FALSE, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

##Track DADA2 pipeline

getN <- function(x) sum(getUniques(x))
track <- as.data.frame(cbind(filter.summary, sapply(dadaRs, getN), rowSums(seqtab.nochim)))
colnames(track) <- c("input", "filtered", "denoisedR", 
    "nonchim")
rownames(track) <- sample.names
head(track)
##Export ASV to remove non-amoA and guided chimera removal

uniquesToFasta(getUniques(seqtab.nochim), "data/seqtab_nochim_uniques.fasta")

#Reference https://github.com/alex-bagnoud/arctic-nitrifiers#2-dada2-pipeline-in-r
#Put databases in the github folder and a read.me in there that states where the files come from (the Nature paper)

#Files for removing non-AmoA and guided chimera removal are found in the database folder. References for them are in the read.me file in that folder.

#command line for usearch to double check/remove non-amoA sequences for AOA
#usearch11.0.667.exe -usearch_global seqtab_nochim_uniques.fasta -db ..\databases\AamoA.db_nr.aln.fasta -id 0.55 -strand both -uc uclust_uniques.txt -matched uniques_match.fasta -notmatched uniques_notmatch.fasta
#usearch11.0.667.exe -uchime2_ref uniques_match.fasta -db ..\databases\AamoA_chimera.ref.db_aln.fasta -notmatched uniques_match_uchimed.fasta -mode balanced -strand plus -mindiv 1.7 -minh 0.15 -chimeras uniques_match_chimeras.fasta

## Taxonomy classification in Qiime2 - showing commands here that were used

#conda activate qiime2-2022.2

#qiime tools import --input-path Desktop/amoa_qiime/uniques_nochim_pooled_match_uchimed.fna --type 'FeatureData[Sequence]' --output-path Desktop/amoa_qiime/uniques_nochim_pooled_match_uchimed.qza

#Removed all dash marks (-) from AamoA.db_nr.aln.fasta file
#qiime tools import --type 'FeatureData[Sequence]' --input-path Desktop/amoa_qiime/AamoA.db_nr.aln.fna --output-path Desktop/amoa_qiime/AamoA.db_nr.aln.qza
#qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path Desktop/amoa_qiime/AamoA.db_nr.aln_taxonomy_qiime.txt --output-path Desktop/amoa_qiime/AamooA.db_nr.aln_taxonomy_qiime.qza
#qiime feature-classifier classify-consensus-vsearch --i-query uniques_nochim_pooled_match_uchimed.qza --i-reference-reads AamoA.db_nr.aln.qza --i-reference-taxonomy AamooA.db_nr.aln_taxonomy_qiime.qza --p-perc-identity 0.8 --p-query-cov 0.5 --p-top-hits-only --p-strand both --p-maxaccepts 1 --p-unassignable-label 'Unassigned' --o-classification taxonomy_convsea.qza
 
##Combine qiime taxonomy, ASV table to create phyloseq object

library(stringr)

asv1 <- data.frame(t(seqtab.nochim))
asv1$sequence <- rownames(asv1)
rownames(asv1) <- NULL
nrow(asv1)

#Several steps and functions from Alex # pipeline
fastaToDf <- function(fastaFile){
    dnaSeq <- readBStringSet(fastaFile)
    fasta_df <- data.frame(header = names(dnaSeq), sequence = paste(dnaSeq))
}

fasta.uchime <- fastaToDf("seqtab_nochim_uniques_match_uchimed.fasta")
asv2 <- merge(fasta.uchime, asv1)

#   sequence                header        AOA10 AOA11 AOA12
#1 AGTATTAAGTTTAGCACCTA... sq238 size=33    0     0     0
#2 AGTATTAAGTTTAGCACCTA... sq209 size=70    0     0     0
#3 AGTATTAAGTTTAGCACCTA... sq279 size=14    0     0     0

#Taxonomy file import and clean up
tax <- read.table("qiime_tax_convse.tsv", header = TRUE, sep = "\t")
tax$Consensus <- NULL
tax <- cbind(tax, str_split_fixed(tax$Taxon, ";", 11))
tax$Taxon <- NULL
tax <- as.data.frame(apply(tax, 2, function(x) gsub("^$|^$", NA, x)))
                              
col_to_remove <- c()

for (col in 1:ncol(tax)) {
    x <- sum(is.na(tax[,col]))/nrow(tax)
    if (x == 1) {
        col_to_remove <- c(col_to_remove, col)
    }
}

if (length(col_to_remove) != 0) {
    tax2 <- tax[,-col_to_remove]
} else {
    tax2 <- tax
}
                              
#Rename columns for TAX levels and generate taxa information to for all levels (expand out since naming for levels is uneven)
names(tax2)[1] <- "ASV"
names(tax2)[-1] <- paste0("TAX_LV_", 1:(ncol(tax2)-1))

for (col in 2:ncol(tax2)) {
    tax2[,col] <- as.character(tax2[,col])
}

for (col in 1:ncol(tax2)) {
    for (row in 1:nrow(tax2)) {
        if (is.na(tax2[row,col])) {
            if (!grepl("OTU", tax2[row,col-1]) & !grepl("unassigned", tax2[row,col-1])) {
                tax2[row,col] <- paste0(tax2[row,col-1], "_unassigned")
            } else {
                tax2[row,col] <- tax2[row,col-1]
            }
        }
    }
}
                              
#Formatting to get it ready to create phyloseq package for ASV, tax and 
asv_tax <-asv2
colnames(asv_tax)[2] <- "ASV"
asv_tax <- merge(asv_tax, tax2, by = "header")
#asv_tax$ASV <- str_split_fixed(asv_tax$header, ";", 2)[,1]
                           
phy.tbl <- asv_tax
row.names(phy.tbl) <- asv_tax$sequence
phy.tbl$sequence <- NULL
phy.tbl$ASV <- NULL
asv.tbl <- phy.tbl[,1:42]
asv.tbl <- t(asv.tbl)
tax.tbl <- phy.tbl[,43:length(phy.tbl)]

#Import soil environmental variables and metadata for samdf
soil.data.all <- read.csv("All_soil_compiled_data.csv")

samples.out <- rownames(asv.tbl)
aoa.plots <- paste0("Plot_", gsub("AOA", "", samples.out))
samdf <- data.frame(PLOT = aoa.plots, stringsAsFactors = FALSE)
samdf <- merge(samdf, soil.data.all, by = "PLOT")
row.names(samdf) <- samples.out

#Create phyloseq object

library(phyloseq)

phy.aoa <- phyloseq(otu_table(asv.tbl, taxa_are_rows=FALSE), tax_table(as.matrix(tax.tbl)), sample_data(samdf))

#Rename the sequence name to "ASV_" # in the matrix
taxa_names(phy.aoa) <- paste0("ASV_", seq(ntaxa(phy.aoa)))
                              
#Remove plots 49, 65 and 71 (low number of reads), plot 20 (only plot for that stand) and ASVs with < 2 occurences
phy.aoa.pruned <- prune_samples(sample_names(phy.aoa) != c("AOA49", "AOA65", "AOA71"), phy.aoa)
phy.aoa.pruned <- prune_samples(sample_names(phy.aoa.pruned) != "AOA20", phy.aoa.pruned)
#Removing plots means that some ASVs are not in multiple samples now; I will filter_taxa to remove any that are not in at least 2 plots.
library(genefilter)
flist <- filterfun_sample(kOverA(2, 1))
phy.aoa.pruned <- filter_taxa(phy.aoa.pruned, flist, TRUE)
saveRDS(phy.aoa.pruned, file = "Phyloseq_AOA_Pruned.rds")
