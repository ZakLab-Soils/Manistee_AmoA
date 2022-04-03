#Some info about making the NMDS etc

#metaMDS to see what works best
nmds.decohel <- metaMDS(psRs.rm496571.less2.hel, autotransform = FALSE, distance = "bray", maxtry = 100)
nmds.decohel.k3 <- metaMDS(psRs.rm496571.less2.hel, autotransform = FALSE, distance = "bray", maxtry = 100, k=3)

#distance matrix for adonis using vegdist

#adonis (permanvoa)
adonis(bc.dist.hlgr ~ N.MIN + SOIL.PH + SPEC.TOTAL.N + VWC.MEAN + TEMP.MEAN + SOIL.CN + SPEC.TOTAL.C, data = env.input.tbl_aoaplots_rm49_65_71)

#environmental data
aoa.plots.rm49_65_71 <- sample_names(psRs.rm496571.less2)
#Data previous collected from Will
soil.compiled.data.tbl <- read.table("soil_compiled_data.txt", sep = "\t")
soil.compiled.data.tbl.aoaplots <- soil.compiled.data.tbl[soil.compiled.data.tbl$PLOT %in% aoa.plots.rm49_65_71,]
