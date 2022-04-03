#Some info about making the NMDS etc

#environmental data
aoa.plots.rm49_65_71.df <- sample_data(psRs.rm496571.less2) %>% data.frame() %>% select("PLOT") %>% mutate_if(is.factor, as.character)
aoa.plots.rm49_65_71 <- aoa.plots.rm49_65_71$PLOT
#Data previous collected from Will
soil.compiled.data.tbl <- read.table("soil_compiled_data.txt", sep = "\t")
soil.pH.tbl <- read_tsv(file = "Manistee_spring_2019_soil_pH.txt") %>%
as_tibble(.) %>%
mutate(PLOT = paste("Plot", PLOT, sep = "_"),
SOIL.PH = as.numeric(SOIL.PH))
soil.compiled.data.tbl.ph <- soil.compiled.data.tbl %>%
inner_join(., soil.pH.tbl, by = "PLOT")

#Remove any plots not in AOA analysis
soil.compiled.data.tbl.ph.aoaplots <- soil.compiled.data.tbl.ph[soil.compiled.data.tbl.ph$PLOT %in% aoa.plots.rm49_65_71,]
soil.compiled.data.tbl.ph.aoaplots <- soil.compiled.data.tbl.ph.aoaplots[order(soil.compiled.data.tbl.ph.aoaplots$PLOT),]

#metaMDS to see what works best
nmds.decohel <- metaMDS(psRs.rm496571.less2.hel, autotransform = FALSE, distance = "bray", maxtry = 100)
nmds.decohel.k3 <- metaMDS(psRs.rm496571.less2.hel, autotransform = FALSE, distance = "bray", maxtry = 100, k=3)

#distance matrix for adonis using vegdist
bc.dist.hlgr <- vegdist(psRs.rm496571.less2.hel, method = "bray")

#adonis (permanvoa)
adonis(bc.dist.hlgr ~N.MIN + SOIL.PH, data =soil.compiled.data.tbl.ph.aoaplots)
