library(dada2)
library(ShortRead)
library(Biostrings)

#setwd for folder where plan to sort and store output for R ("C:/**/**/amoa-AOB/")

set.seed(496571)

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

cutadapt <- "~/Python/Python310/Scripts/cutadapt"
path.cut <- file.path(path, "cutadapt_processed")

if(!dir.exists(path.cut)) dir.create(path.cut)

fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD <- "GGGGTTTCTACTGGTGGT"
REV <- "CCCCTCKGSAAAGCCTTCTTC"
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

rbind(FWD.FowardReads = sapply(FWD.orients, primerHits, fn=fnFs.cut[[1]]), FWD.ReverseReads = sapply(FWD.orients, primerHits, fn=fnRs.cut[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits, fn=fnFs.cut[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn=fnRs.cut[[1]]))
#                Forward Complement Reverse RevComp
#FWD.FowardReads        0          0       0       0
#FWD.ReverseReads       0          0       0       0
#REV.ForwardReads       0          0       0       0
#REV.ReverseReads       0          0       0       0

cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))

## Filter and trim for length, quality and size; Calculate error rate.

get.sample.name <- function(fname) strsplit(basename(fname), "-|_" )[[1]][3]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))


#Still re-editing
filter.summary <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, truncLen = c(245,220), minLen = c(200,200), maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
#filter.summary <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2,2), truncQ = 2, minLen = c(200,200), rm.phix = TRUE, compress = TRUE, multithread = FALSE)

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

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)

#Investigation has shown that AOB merged reads are 452 bp long so I only use those ASVs with a +/- 1bp
seqtab.452 <- seqtab[,nchar(colnames(seqtab)) %in% 451:453]
saveRDS(seqtab.452, file="Seqtab_452_AOB.rds")

#Remove ASVs from the data set that are not in at least 2 plots
ASVs.multisamp <- integer()
for (i in 1:dim(seqtab.452)[2]) {if(sum(seqtab.452[,i]==0) != 38){tmp.paste <- i; ASVs.multisamp <- append(ASVs.multisamp, tmp.paste) }}
seqtab.multisamp <- seqtab.452[, ASVs.multisamp]
dim(seqtab.multisamp)

seqtab.nochim <- removeBimeraDenovo(seqtab.multisamp, method = "pooled", multithread = FALSE, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

##Track DADA2 pipeline

getN <- function(x) sum(getUniques(x))
track <- as.data.frame(cbind(filter.summary, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track) 

#uniquesToFasta(getUniques(seqtab.nochim), ids=paste0("ASV_", seq(length(getUniques(seqtab.nochim)))), "seqtab_nochim_uniques_renamed.fasta")
#Used Aigle paper database classification (D1-D17) - Note that these are already trimmed to 414 since some of reads do not exist longer than that.
## Taxonomy classification in Qiime2 - showing commands here that were used

#conda activate qiime2-2022.2

#qiime tools import --input-path Desktop/amoa_qiime/seqtab_nochim_uniques.fasta --type 'FeatureData[Sequence]' --output-path Desktop/amoa_qiime/seqtab_nochim_uniques.qza

#Removed all dash marks (-) from AamoA.db_nr.aln.fasta file
#qiime tools import --type 'FeatureData[Sequence]' --input-path Desktop/amoa_qiime/AigleDB_ref_reads.fasta --output-path Desktop/amoa_qiime/AigleDB_ref_reads.qza
#qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path Desktop/amoa_qiime/AigleDB_taxonomy_qiime.txt --output-path Desktop/amoa_qiime/AigleDB_taxonomy_qiime.qza
#qiime feature-classifier classify-consensus-vsearch --i-query seqtab_nochim_uniques.qza --i-reference-reads AigleDB_ref_reads.qza --i-reference-taxonomy AigleDB_taxonomy_qiime.qza --p-perc-identity 0.8 --p-query-cov 0.6 --p-top-hits-only --p-strand both --p-maxaccepts 1 --p-unassignable-label 'Unassigned' --o-classification taxonomy_convsea.qza
 
##Combine qiime taxonomy, ASV table to create phyloseq object

#asv.aob <- data.frame(t(seqtab.nochim))
#asv.aob$sequence <- rownames(asv.aob)
#rownames(asv.aob) <- NULL
#asv.aob$ASV <- paste0("ASV_", seq(length(getUniques(seqtab.nochim))))

#Taxonomy file import and clean up
#tax <- read.table("AOB_AigleDB_taxonomy.tsv", header = TRUE, sep = "\t")
#tax$Consensus <- NULL
#tax <- cbind(tax, str_split_fixed(tax$Taxon, ";", 3))
#tax$Taxon <- NULL
#tax <- as.data.frame(apply(tax, 2, function(x) gsub("^$|^$", NA, x)))
                              
#names(tax)[1] <- "ASV"
#names(tax)[-1] <- c("Lineage", "Genus", "Strain")
                              
#Formatting to get it ready to create phyloseq package for ASV, tax and 
#asv.aob.tax <-asv.aob
#asv.aob.tax <- merge(asv.aob.tax, tax, by = "ASV")
#row.names(asv.aob.tax) <- asv.aob.tax$ASV
#asv.aob.tax$sequence <- NULL
#asv.aob.tax$ASV <- NULL

#asv.tbl <- t(asv.aob.tax[,1:39])
#tax.tbl <- asv.aob.tax[,40:length(asv.aob.tax)]



#Import soil environmental variables and metadata for samdf
soil.data.all <- read.csv("All_soil_compiled_data.csv")
samples.out <- rownames(seqtab.nochim)
#samples.out <- rownames(asv.tbl)
aob.plots <- paste0("Plot_", gsub("AOB", "", samples.out))
samdf <- data.frame(PLOT = aob.plots, stringsAsFactors = FALSE)
samdf <- merge(samdf, soil.data.all, by = "PLOT")
row.names(samdf) <- samples.out

#Create phyloseq object; Keep read names until after clustering.

library(phyloseq)
phy.aob <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samdf))
#phy.aob <- phyloseq(otu_table(asv.tbl, taxa_are_rows=FALSE), tax_table(as.matrix(tax.tbl)), sample_data(samdf))                               

#Removing plots means that some ASVs are not in multiple samples now; I will filter_taxa to remove any that are not in at least 2 plots.
library(genefilter)
flist <- filterfun_sample(kOverA(2, 1))
phy.aob.pruned <- filter_taxa(phy.aob, flist, TRUE)
saveRDS(phy.aob.pruned, file = "Phyloseq_AOB_Pruned.rds")

