library(dada2)
library(ShortRead)
library(Biostrings)

#setwd("C:/Users/17348/Desktop/amoa-AOA/Fastq/dada2-onlyafew") - I'm mentioning this since it affects the path and path.cut variables and functions.R

source("functions.R")
path <- "data"
list.files(path)
                    
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#AOA primers
FWD <- "ATGGTCTGGCTWAGACG"
REV <- "GCCATCCATCTGTATGTCCA"

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

fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

REV.RC <- dada2:::rc(REV)
FWD.RC <- dada2:::rc(FWD)

R1.flags <- paste("-g", FWD, "-a", REV.RC)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

#See how much the primers appear in the data

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
