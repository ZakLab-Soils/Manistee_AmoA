# *Ammonium availability and pH control the contributions of ammonia-oxidizing bacteria and archaea to nitrification*
## Jennifer Wen<sup>1</sup>, Rima Upchurch<sup>1</sup>, Donald R. Zak<sup>1,2</sup>
1 School for Environment and Sustainability, University of Michigan, Ann Arbor, MI, USA 
2 Department of Ecology & Evolutionary Biology, University of Michigan, Ann Arbor, MI, USA 

***************

This is a repository for code and data files for processing data for the above named manuscript. Sequences of *amoA* for this experiment are available at NCBI in the Manistee Forest Soil Nitrogen Mineralization Gradient Study	BioProject (PRJNA714922) under the SRA accession numbers SRR25249670-SRR25249708 for AOA and SRR25407544 – SRR25407582 for AOB.

The R scripts used for analyzing are included in the `code` folder with subfolders for AOA and AOB. The `data` folder contains the soil environmental data (Nmin, pH, etc.), the qPCR data and treefiles used in the manuscript. The `database` folder has the files used for usearch chimera removal and qiime taxonomic assignment for AOA and the references used the phylogenetic trees for AOA and AOB *amoA*.

Commands for usearch and qiime for AOA assignments are listed as comments in the AOA1_dada2.r script. Phylogenetic trees were constructed by aligning the amoA sequences to the relevant databases (AOA, AOB) with MAFFT (https://mafft.cbrc.jp/alignment/server/index.html) using the L-INS-i method. Alignments were trimmed with msaTrim in R. Maximum likelihood treefiles were inferred using IQ-TREE with GTR+F+I+Gamma4 at 1000 bootstraps with a perturbation of 0.1 (iqtree2.exe -s mafft.alignment -bb 1000 -alrt 1000 -nt AUTO -pers 0.1 -m GTR+F+I+G4). Tree visualizations and annotations were created from iTOL (https://itol.embl.de/).

# Addition to previous release
# *Comammox Nitrospira Abundance And Contribution To Nitrification In A Forested Landscape*
## Donald R. Zak<sup>1,2</sup> & Rima Upchurch<sup>1</sup> 
1 School for Environment and Sustainability, University of Michigan, Ann Arbor, MI, USA 
2 Department of Ecology & Evolutionary Biology, University of Michigan, Ann Arbor, MI, USA 

***************

The newest release for this repository includes code and data files for processing data for the above named note submitted to _Ecosphere_. Included in the `code` folder is a subfolder for Comammox. The `data` folder now contains the qPCR data for the Comammox analyses. 
