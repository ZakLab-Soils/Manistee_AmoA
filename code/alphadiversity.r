#Rarefy phyloseq
phy.Rs.rare <- rarefy_even_depth(phy.Rs, rngseed =496571, sample.size=min(sample_sums(phy.Rs)), replace=F)
#Prettier name for Stand
sample_data(phy.Rs.rare)$STAND.CLEAN <- gsub("Stand_*", "", fixed=FALSE, sample_data(phy.Rs.rare)$STAND)

#Alpha diversity
plot.rich.rare.data <- plot_richness(phy.Rs.rare, x="STAND.CLEAN", measures = "Shannon", color = "STAND.CLEAN")$data
plot.rich.rare.data <- plot.rich.rare.data[order(plot.rich.rare.data$STAND.LEVEL),]

plot.chao.data <- plot_richness(phy.Rs.rare, x="STAND.CLEAN", measures = "Chao1", color = "STAND.CLEAN")$data
plot.chao.data <- plot.chao.data[order(plot.chao.data$STAND.CLEAN),]

#Create Nmin level (Stand ranking low to high 58, 7, 6, 41, 24, 100, 22)
plot.rich.rare.data$STAND.LEVEL <- c(6,6,6,6,6,7,7,7,7,7,7,5,5,5,5,5,5,4,4,4,4,4,1,1,1,1,1,3,3,3,3,3,3,2,2,2,2,2)
plot.chao.data$STAND.LEVEL <- c(6,6,6,6,6,7,7,7,7,7,7,5,5,5,5,5,5,4,4,4,4,4,1,1,1,1,1,3,3,3,3,3,3,2,2,2,2,2)
nmds.colors2 <- c("coral4", "coral3", "coral2", "coral1", "dodgerblue2", "dodgerblue3", "dodgerblue4")

shannon.plot <- ggplot(data=plot.rich.rare.data, aes(x=reorder(STAND.CLEAN, STAND.LEVEL), y=value, fill = reorder(STAND.CLEAN, STAND.LEVEL))) + geom_point(aes(), colour = "black", pch = 21, size = 3) + xlab(paste0("STAND")) + ylab(paste0("Alpha Diversity")) + guides(fill=guide_legend(title="Stand"))+scale_fill_manual(values = nmds.colors2)+theme(text=element_text(size = 18)) + ggtitle("Shannon Richness of AOA")

chao.plot <- ggplot(data=plot.chao.data, aes(x=reorder(STAND.CLEAN, STAND.LEVEL), y=value, fill = reorder(STAND.CLEAN, STAND.LEVEL))) + geom_point(aes(), colour = "black", pch = 21, size = 3) + xlab(paste0("STAND")) + ylab(paste0("Alpha Diversity")) + guides(fill=guide_legend(title="Stand"))+scale_fill_manual(values = nmds.colors2)+theme(text=element_text(size = 18)) + ggtitle("Chao1 of AOA")

