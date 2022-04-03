library(dada2)
library(ShortRead)
library(Biostrings)

#setwd("C:/Users/17348/Desktop/amoa-AOA/Fastq/dada2-onlyafew") - I'm mentioning this since it affects the path and path.cut variables and functions.R

source("functions.R")
path <- "data"
list.files(path)
                    
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#Check quality of reads - these can be saved
plotQualityProfile(fnFs)
plotQualityProfile(fnRs)


#AOA primers
FWD <- "ATGGTCTGGCTWAGACG"
REV <- "GCCATCCATCTGTATGTCCA"
#AOB
#FWD <- "GGGGTTTCTACTGGTGGT"
#REV <- "CCCCTCKGSAAAGCCTTCTTC"

allOrients <- function(primer) {
require(Biostrings)
dna <- DNAString(primer)
orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
RevComp = reverseComplement(dna))
return(sapply(orients, toString))
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN=0, multithread = FALSE)

cutadapt <- "/Users/17348/AppData/Local/Programs/Python/Python310/Scripts/cutadapt"
path.cut <- file.path(path, "cutadapt_processed")

if(!dir.exists(path.cut)) dir.create(path.cut)

fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

REV.RC <- dada2:::rc(REV)
FWD.RC <- dada2:::rc(FWD)

R1.flags <- paste("-g", FWD, "-a", REV.RC)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

#See how much the primers appear in the data

#primerHits function included in the "functions.R" file
#primerHits <- function(primer, fn) {
#nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#return(sum(nhits > 0))
#}

rbind(FWD.FowardReads = sapply(FWD.orients, primerHits, fn=fnFs.filtN[[1]]), FWD.ReverseReads = sapply(FWD.orients, primerHits, fn=fnRs.filtN[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits, fn=fnFs.filtN[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn=fnRs.filtN[[1]]))

#                Forward Complement Reverse RevComp
#FWD.FowardReads        x          x       x       x
#FWD.ReverseReads       x          x       x       x
#REV.ForwardReads       x          x       x       x
#REV.ReverseReads       x          x       x       x


for(i in seq_along(fnFs)){system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "-o", fnFs.cut[i], "-p", fnRs.cut[i], fnFs.filtN[i], fnRs.filtN[i]))}
#Cutadapt Command line parameters: -g ATGGTCTGGCTWAGACG -a TGGACATACAGATGGATGGC -G GCCATCCATCTGTATGTCCA -A CGTCTWAGCCAGACCAT -n 2 -o data/cutadapt_processed/AOA7_S1_L001_R1_001.fastq.gz #-p data/cutadapt_processed/AOA7_S1_L001_R2_001.fastq.gz data/filtN/AOA7_S1_L001_R1_001.fastq.gz data/filtN/AOA7_S1_L001_R2_001.fastq.gz

#Check if it worked

rbind(FWD.FowardReads = sapply(FWD.orients, primerHits, fn=fnFs.cut[[1]]), FWD.ReverseReads = sapply(FWD.orients, primerHits, fn=fnRs.cut[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits, fn=fnFs.cut[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn=fnRs.cut[[1]]))
#                Forward Complement Reverse RevComp
#FWD.FowardReads        0          0       0       0
#FWD.ReverseReads       0          0       0       0
#REV.ForwardReads       0          0       0       0
#REV.ReverseReads       0          0       0       0

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#Filter out samples based on number of N's, expected error (EE) and quality score (truncQ), length, 
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
#Upcoming edit - change minLen to show minLen= c(200,200).
filter.summary <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2,2), truncQ = 2, minLen = 100, rm.phix = TRUE, compress = TRUE, multithread = FALSE)

#Error rate processing
# Quantify error rate
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)

# Plot error rates
errF.plot <- plotErrors(errF, nominalQ = TRUE)

#Dereplicate
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Inference 
dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
dadaRs <- dada(derepRs, err = errR, multithread = FALSE)

#Merge pairs for AOA non-overlapping so use 
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate = TRUE, verbose = TRUE)
head(mergers[[1]])

#Sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
seqtab.table <- table(nchar(getSequences(seqtab)))
saveRDS(seqtab.table, file = "data/seqtab.table.rds")

#After mergers added a string of Ns so I wanted to remove them in the colnames of the seqtab
#newcolnames <- gsub("NNNNNNNNNN", "----------", colnames(seqtab))
#colnames(seqtab) <- newcolnames

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Track DADA2 pipeline
#getN <- function(x) sum(getUniques(x))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
#track <- as.data.frame(cbind(filter.summary, sapply(dadaFs, getN), rowSums(seqtab.nochim)))
#colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
#rownames(track) <- sample.names

getN <- function(x) sum(getUniques(x))
track <- cbind(filter.summary, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
    getN), rowSums(seqtab.nochim))
#If saving track file add the 'as.data.frame' function around cbind command
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
    "nonchim")
rownames(track) <- sample.names
head(track)

#Exporting sequences
uniquesToFasta(getUniques(seqtab.nochim), "data/seqtab_nochim_uniques.fasta")
#write.table(t(seqtab.nochim), "seqtab_nochim_asv_table.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

#command line for usearch to double check/remove non-amoA sequences for AOA (none for AOB yet)
#usearch11.exe -usearch_global data\seqtab_nochim_uniques.fasta -db databases\AamoA.db_all.seqs.fasta -id 0.55 -strand plus -uc data\uclust_report.txt -matched data\uniques_nochim_match.fasta -notmatched uniques_nochim_nomatch.fasta

#usearch11.0.667_win32.exe -uchime2_ref aoa\Fastq\data\data\uniquesRs_nochim_pooled_match.fasta -db aoa\Fastq\databases\AamoA_chimera.ref.db_aln.fasta -notmatched aoa\Fastq\data\data\uniquesRs_nochim_pooled_match_uchimed.fasta -mode balanced -strand plus -mindiv 1.7 -minh 0.15 -chimeras aoa\Fastq\data\data\uniquesRs_nochim_pooled_match_chimeras.fasta

#C:\Users\17348\Desktop\amoa-AOA>usearch11.0.667_win32.exe -uchime2_ref aoa\Fastq\data\data\uniques_nochim_pooled_newest_matchnoNs.fasta -db aoa\Fastq\databases\AamoA_chimera.ref.db_aln.trim.fasta -notmatched aoa\Fastq\data\data\uniques_nochim_pooled_newest_matchnoNs_uchimed02.fasta -mode balanced -mindiffs 4 -mindiv 2 -minh 0.15 -strand plus -chimeras aoa\Fastq\data\data\uniques_nochim_pooled_newest_matchnoNs_chimeras02.fasta

#taxonomy call in qiime
#two ways to try it out

#qiime2

conda activate qiime2-2022.2

#import unique read data post chimera checking (outside of dada2)
qiime tools import --input-path Desktop/amoa_qiime/uniques_nochim_pooled_match_uchimed.fna --type 'FeatureData[Sequence]' --output-path Desktop/amoa_qiime/uniques_nochim_pooled_match_uchimed.qza

#Importing taxonomy info
#Import ref seqs - had to remove (---) that were in 3 different locations in each sequence
#saved the AamoA.db_nr.aln.fasta as a .fna file after removing (---)
#qiime tools import --type 'FeatureData[Sequence]' --input-path Desktop/amoa_qiime/AamoA.db_nr.aln.fna --output-path Desktop/amoa_qiime/AamoA.db_nr.aln.qza

#Import qiime taxonomy file
#qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path Desktop/amoa_qiime/AamoA.db_nr.aln_taxonomy_qiime.txt --output-path Desktop/amoa_qiime/AamooA.db_nr.aln_taxonomy_qiime.qza

#create classifer file using full length reads (did not do primer selection for ref.reads)
#qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads Desktop/amoa_qiime/AamoA.db_nr.aln.qza --i-reference-taxonomy Desktop/amoa_qiime/AamooA.db_nr.aln_taxonomy_qiime.qza --o-classifier Desktop/amoa_qiime/classifier.qza

#assign taxonomy
#qiime feature-classifier classify-sklearn --i-classifier Desktop/amoa_qiime/classifier.qza --i-reads Desktop/amoa_qiime/uniques_nochim_pooled_match_uchimed.qza  --o-classification Desktop/amoa_qiime/AamoA_taxonomy.qza                    

#get metadat ready and output taxonomy classificatoin
#qiime metadata tabulate --m-input-file Desktop/amoa_qiime/AamoA_taxonomy.qza --o-visualization Desktop/amoa_qiime/taxonomy.qzv

#qiime tools export --input-path Desktop/amoa_qiime/AamoA_taxonomy.qza --output-path Desktop/amoa_qiime/taxonomy_output

#Tried to do a quick QC step before classify sklearn
#qiime quality-control exclude-seqs --i-query-sequences uniques_nochim_pooled_match_uchimed.qza --i-reference-sequences AamoA.db_nr.aln.qza --p-method vsearch --p-perc-identity 0.8 --o-sequence-hits QChits.qza --o-sequence-misses QCmiss.qza

#qiime feature-classifier classify-sklearn --i-classifier classifier.qza --i-reads QChits.qza --o-classification AamoA_taxonomy_afterQC.qza

#Tried a different classification because confidence was a bit "odd"; this was similar to qiime1 classification
#qiime feature-classifier classify-consensus-vsearch --i-query uniques_nochim_pooled_match_uchimed.qza --i-reference-reads AamoA.db_nr.aln.qza --i-reference-taxonomy AamooA.db_nr.aln_taxonomy_qiime.qza --p-perc-identity 0.8 --p-query-cov 0.5 --p-top-hits-only --p-strand both --p-maxaccepts 1 --p-unassignable-label 'Unassigned' --o-classification taxonomy_convsea.qza
 
#to "flatten" from alex's pipeline to make it read it correctly 
#awk -v RS='>' -v FS="\n" -v OFS="" -v ORS="" '{ if (NR>1) { printf ">%s\n", $1; $1=""; printf "%s\n", $0}}' uniques_nochim_pooled_match_uchimed.fasta > uniques_nochim_pooled_match_uchimed_flat.fasta

#Importing taxonomy file to combine with a subset of asv table that corresponds to unique asvs that were not removed via usearch and uchime
asv1 <- seqtabRs.nochim
asv1$sequence <- rownames(asv1)
rownames(asv1) <- NULL
nrow(asv1)

#function from Alex's pipeline
fastaToDf <- function(fastaFile){
    dnaSeq <- readBStringSet(fastaFile)
    fasta_df <- data.frame(header = names(dnaSeq), sequence = paste(dnaSeq))
}

fastaRs <- fastaToDf("uniquesRs_nochim_pooled_match_uchimed_flat.fasta")
asvRs2 <- merge(fastaRs, asvRs)
asvRs3 <- asvRs2[,c(2:ncol(asvRs2),1)]

#Cutting to the chase a bit - after a lot of testing, I ended up wanting things set up for the phyloseq package
#If you take the output from the asvRs2 which now has been edited to remove chimeras, non-amoa and classified (you can adjust which remains based on classification prior to this step or not)

#   sequence                header        AOA10 AOA11 AOA12
#1 AGTATTAAGTTTAGCACCTA... sq238 size=33    0     0     0
#2 AGTATTAAGTTTAGCACCTA... sq209 size=70    0     0     0
#3 AGTATTAAGTTTAGCACCTA... sq279 size=14    0     0     0

#Time for taxonomy from qiime - used the vsearch consensus here
taxRs <- read.table("Rs_tax_convse.tsv", header = FALSE, sep = "\t")
taxRs$V3 <- NULL

library(stringr)
taxRs <- cbind(taxRs, str_split_fixed(taxRs$V2, ";", 11))
taxRs$V2 <- NULL
#Put NA for empty cells
taxRs2 <- as.data.frame(apply(taxRs, 2, function(x) gsub("^$|^$", NA, x)))

#Remove any columns that have only NAs
col_to_remove <- c()

for (col in 1:ncol(taxRs2)) {
    x <- sum(is.na(taxRs2[,col]))/nrow(taxRs2)
    if (x == 1) {
        col_to_remove <- c(col_to_remove, col)
    }
}

if (length(col_to_remove) != 0) {
    taxRs3 <- taxRs2[,-col_to_remove]
} else {
    taxRs3 <- taxRs2
}
                              
#Rename columns
names(taxRs3_attempt2)[1] <- "ASV"
names(taxRs3_attempt2)[-1] <- paste0("1", 1:(ncol(taxRs3_attempt2)-1), "_tax")

#Character for taxonomy
for (col in 2:ncol(taxRs3)) {
    taxRs3[,col] <- as.character(taxRs3[,col])
}

#Rename the cells with taxonomy that have NA to whatever was in the previous cell. This is because the "levels" are different to get to OTU - it might take like out to column 11 to get to OTU or column 3.
for (col in 1:ncol(taxRs3)) {
    for (row in 1:nrow(taxRs3)) {
        if (is.na(taxRs3[row,col])) {
            if (!grepl("OTU", taxRs3[row,col-1]) & !grepl("unassigned", taxRs3[row,col-1])) {
                taxRs3[row,col] <- paste0(taxRs3[row,col-1], "_unassigned")
            } else {
                taxRs3[row,col] <- taxRs3[row,col-1]
            }
        }
    }
}

#Remove unneeded row 1 (leftover from taxonomy naming)
 taxRs3 <- taxRs3[-1,]
                              
#get it ready for phyloseq package
asvRs_taxmerger <-asvRs2
asvRs_taxmerger$ASV <- sub(" size=[0-9]* ", "", asvRs_taxmerger$header)
asvRs_taxmerger$header <- NULL
asvRs.tax <- merge(asvRs_taxmerger, taxRs3, by.x = "ASV", by.y = "ASV")

#otu_table needs to look like seqtab.nochim output, but this has only the sequences remaining after usearch, uchime and taxonomy
asvRs.otu.tbl <- asvRs.tax
row.names(asvRs.otu.tbl) <- asvRs.tax$sequence
asvRs.otu.tbl$sequence <- NULL
asvRs.otu.tbl$ASV <- NULL
asvRs.otu.tbl <- asvRs.otu.tbl[,1:42]
asvRs.otu.tbl.t <- t(asvRs.otu.tbl)

#tax_table
taxRs.tax.tbl <- asvRs.tax_norel
row.names(taxRs.tax.tbl) <- taxRs.tax.tbl$sequence
taxRs.tax.tbl <- taxRs.tax.tbl[,45:length(taxRs.tax.tbl)]
#taxRs.tax.tbl.t <- t(taxRs.tax.tbl)

#Sample data (this info included in "functions.R" that Will made
samples.out <- rownames(asvRs.otu.tbl.t)
samdf <- data.frame(SAMPLE.ID = samples.out, stringsAsFactors = FALSE)
row.names(samdf) <- samples.out

# Get sample type
samdf$TYPE <- sapply(samdf$SAMPLE.ID, sample_type, simplify = TRUE)
# Get sample plot
samdf$PLOT <- mapply(sample_plot, samdf$SAMPLE.ID, samdf$TYPE, SIMPLIFY = TRUE)
# Get sample stand
samdf$STAND <- mapply(sample_stand, samdf$PLOT, samdf$TYPE, SIMPLIFY = TRUE)
# Update plot
samdf$PLOT <- mapply(sample_plot_update, samdf$PLOT, samdf$TYPE, SIMPLIFY = TRUE)
# Update stand
samdf$STAND <- mapply(sample_stand_update, samdf$STAND, samdf$TYPE, SIMPLIFY = TRUE)

#phyloseq object                              
psRs.otu.t <- otu_table(asvRs.otu.tbl.t, taxa_are_rows = FALSE)
psRs.tax <- tax_table(as.matrix(taxRs.tax.tbl))
psRs.samdf <- sample_data(samdf2)
psRs <- phyloseq(psRs.otu.t, psRs.tax, psRs.samdf)
psRs.refseqs <- Biostrings::DNAStringSet(taxa_names(psRs))
names(psRs.refseqs) <- taxa_names(psRs)
psRs.refseqs <- Biostrings::DNAStringSet(taxa_names(psRs))
#Rename the sequence name to "ASV_" # in the matrix
psRs <- merge_phyloseq(psRs, psRs.refseqs)taxa_names(psRs) <- paste0("ASV_", seq(ntaxa(psRs)))
                              
#Remove plots 49, 65 and 71 (low number of reads)
psRs.rm496571 <- prune_samples(sample_names(psRs) != c("AOA49", "AOA65", "AOA71"), psRs)
#Remove any ASVs that have less than 2 occurences 
psRs.rm496571.less2 <- prune_taxa(taxa_sums(psRs.rm496571) > 2, psRs.rm496571)

#Transform otu abundances for analysis with deconstand
#Relative Abundance
psRs.rm496571.less2.tot <- decostand(phylose::otu_table(psRs.rm496571.less2), method = "total")
psRs.rm496571.less2.hel <- deconstand(phyloseq::otu_table(psRs.rm496571.less2), method = "hellinger")
                              
