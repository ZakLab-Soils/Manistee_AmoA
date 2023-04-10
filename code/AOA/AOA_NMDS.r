#Some info about making the NMDS etc

#Scree plot to test for stree to see k value. The point where it pivot (“elbow”) is what we use for k in adonis.
NMDS.scree <- function(x) {
plot(rep(1,10), replicate(10, metaMDS(x, autotransform = F, k= 1)$stress), xlim = c(1,10), ylim=c(0,0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
for (i in 1:10) {
 points(rep(i +1,10), replicate(10, metaMDS(x, autotransform = F, k = i +1)$stress))
}
}
NMDS.scree(phy.aoa.hel.dist)


library(tidyverse)
#environmental data
aoa.plots.rm49_65_71.df <- sample_data(psRs.rm496571.less2) %>% data.frame() %>% select("PLOT") %>% mutate_if(is.factor, as.character)
aoa.plots.rm49_65_71 <- aoa.plots.rm49_65_71df$PLOT
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
#Remove plot 20 (only member of Stand 9
scdtpar_nostand9 <- soil.compiled.data.tbl.ph.aoaplots[soil.compiled.data.tbl.ph.aoaplots$PLOT == "Plot_20",]


#transform otu_table #can do this for "total" method for relative abundance
#psRs.rm496571.less2.no20.hel <- decostand(otu_table(psRs.rm496571.less2.no20), method = "hellinger")
ps.rm496571.less2.hel.no20 <- ps.rm496571.less2.hel[row.names(ps.rm496571.less2.hel) != "AOA20",]

#distance matrix for adonis using vegdist
bc.dist.hlgr.no20 <- vegdist(psRs.rm496571.less2.no20.hel, method = "bray")

#metaMDS for k = 2 and k =3 . Scree plot (data not shown) implied 3 but maybe 2 which is pretty accurate as k=2 stress ~0.15 and k=3 stress ~0.08
nmds.decohel <- metaMDS(bc.dist.hlgr.no20, autotransform = FALSE, distance = "bray", maxtry = 100)
nmds.decohel.k3 <- metaMDS(bc.dist.hlgr.no20, autotransform = FALSE, distance = "bray", maxtry = 100, k=3)

#Assign level to Nmin to test for ANOSIM - but ended up doing mantel test anyway
ano.nlevel$Nlevel <- c("M", "M", "M", "M", "M", "M", "M", "M", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "M", "M", "M", "M", "M", "H", "H", "H", "M", "H", "H", "M", "M")
anosim(bc.dist.hlgr.no20, ano.nlevel$Nlevel)

#Mantel test for continuous variables
dist.soilph <- dist(scdtpar_nostand9$SOIL.PH, method = "euclidean")
dist.nmin <- dist(scdtpar_nostand9$N.MIN, method = "euclidean")
soilph.mantel <- mantel(bc.dist.hlgr.no20, dist.soilph, method = "spearman", permutations = 9999, na.rm = TRUE)
nmin.mantel <- mantel(bc.dist.hlgr.no20, dist.nmin, method = "spearman", permutations = 9999, na.rm = TRUE)

#scale.env <- scale(scdtpar_nostand9[,c(3:7,10)], center = TRUE, scale = TRUE)
#dist.env <- dist(scale.env, method = "euclidean")
#env.mantel <- mantel(bc.dist.hlgr.no20, dist.env, method = "spearman", permutations = 9999, na.rm = TRUE)

#No gonna lie - distance was significant but since all of our env are on a gradient, that is pretty much an obvious conclusion
#geo <- data.frame(scdtpar_nostand9$LON, scdtpar_nostand9$LAT)
#d.geo <- distm(geo, fun = distHaversine)
#dist.geo <- as.dist(d.geo)
#geo.mantel <- mantel(bc.dist.hlgr.no20, dist.geo, method = "spearman", permutations = 9999, na.rm = TRUE)

#adonis (permanvoa)
adonis2(bc.dist.hlgr.no20 ~SOIL.PH + N.MIN + SPEC.TOTAL.N + N.INORG , data =scdtpar_nostand9)

#Pull out scores to plot - used with no plot 20 here but can do it with it in
data.scores.no20 <- as.data.frame(scores(nmds.decohel.k3))
data.scores.no20$sample <- rownames(data.scores)
data.scores.no20$STAND <- as.factor(as.integer(gsub("Stand_", "", scdtpar_nost9$STAND)))

#color coded - remove a blue one if taking out plot 20
nmds.colors <- c("coral4", "coral3", "coral2", "coral1", "dodgerblue4", "dodgerblue3", "dodgerblue2", "dodgerblue1")

#Heads up this the ggplot version but I also can do a 3d version (had some issues with coloring those...) 
#Also this is for 2 v 3 but I have this in 1 v 2 and 1 v3 and just as a k=2.
p.nmdsB.no20 <- ggplot(data = data.scores.no20, aes(x = NMDS2, y = NMDS3, fill = STAND)) + geom_point(aes(), colour = "black", pch = 21, size = 3) + xlab(paste0("NMDS2")) + ylab(paste0("NMDS3")) + stat_ellipse(aes(x = NMDS2, y= NMDS3, color = STAND), level = 0.5) + scale_fill_manual(values = nmds.colors) + coord_fixed() + theme_linedraw(base_size = 18) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(face = "bold")) + theme(text=element_text(size = 20)) + ggtitle("NMDS of ASV amoA of AOA")

#You can use envfit to fit environmental variables to your ordination and display them as vectors
#env_fit_no20 <- envfit(nmds.decohel.no20$points, scdtpar_nostand9, perm =999)
#env_fit_no20_df <- as.data.frame(env_fit_no20$vectors$arrows*sqrt(env_fit_no20$vectors$r))
#env_fit_no20_df$vector <- rownames(env_fit_no20_df)

#3 dimentions (k = 3) for NMDS
#env_fit_no20_k3 <- envfit(nmds.decohel.no20$points, scdtpar_nostand9, perm =999, choices = c(1:3))
#env_fit_no20_k3_df <- as.data.frame(env_fit_no20_k3$vectors$arrows*sqrt(env_fit_no20_k3$vectors$r))
#env_fit_no20_k3_df$vector <- rownames(env_fit_no20_k3_df)

#Add envfit vectors for N.MIN and pH
#p.nmdsB.no20 + geom_segment(data= env_fit_no20_k3_df[c(1,8),], aes(x = 0, xend = env_fit_no20_k3_df$MDS2[c(1,8)], y = 0, yend=env_fit_no20_k3_df$MDS3[c(1,8)]), inherit.aes = FALSE, arrow = arrow(length = unit(0.25, "cm")), colour = "grey") + theme_linedraw(base_size = 18) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(face = "bold")) + theme(text=element_text(size = 20)) + geom_text(data=env_fit_no20_k3_df[c(1,8),], aes(env_fit_no20_k3_df$MDS2[c(1,8)], y=env_fit_no20_k3_df$MDS3[c(1,8)], label = vector), inherit.aes = FALSE, size = 3) + ggtitle("NMDS of ASV amoA of AOA")
