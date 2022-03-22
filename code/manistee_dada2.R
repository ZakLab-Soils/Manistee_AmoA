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

