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

#Boxplots
shannon.boxplot <- ggplot(data=plot.rich.rare.data, aes(x=reorder(STAND.CLEAN, STAND.LEVEL), y=value, fill = reorder(STAND.CLEAN, STAND.LEVEL))) + geom_boxplot(alpha=0.7) + xlab(paste0("STAND")) + ylab(paste0("Alpha Diversity")) + guides(fill=guide_legend(title="Stand"))+scale_fill_manual(values = nmds.colors2)+theme(text=element_text(size = 18)) + ggtitle("Shannon of AOA")

chao.boxplot <- ggplot(data=plot.chao.data, aes(x=reorder(STAND.CLEAN, STAND.LEVEL), y=value, fill = reorder(STAND.CLEAN, STAND.LEVEL))) + geom_boxplot(alpha=0.7) + xlab(paste0("STAND")) + ylab(paste0("Alpha Diversity")) + guides(fill=guide_legend(title="Stand"))+scale_fill_manual(values = nmds.colors2)+theme(text=element_text(size = 18)) + ggtitle("Chao1 of AOA")

#Rarefaction by Stand
#Made an "otu table" for each stand



rarefact.stand.named <- rarecurve(ASV.stand.mean.pw.integer.t, col = nmds.colors, step = 100, lwd=2, ylab = "ASVs", label = F, main = "Rarefaction Curve by Stand")
names(rarefact.stand.named) <- rownames(ASV.stand.mean.pw.integer.t)
protox <- mapply(FUN = function(x, y) {
    mydf <- as.data.frame(x)
    colnames(mydf) <- "value"
    mydf$SampleID <- y
    mydf$subsample <- attr(x, "Subsample")
    mydf
}, x = rarefact.stand.named, y = as.list(names(rarefact.stand.named)), SIMPLIFY = FALSE)

xy <- do.call(rbind, protox)
rownames(xy) <- NULL

ggplot(xy, aes(x=subsample, y=value, color=SampleID)) +theme_bw()+geom_line(size = 0.8)+geom_vline(xintercept = raremaxstand, color = "red", linetype = "dashed") +ggtitle("Rarefaction Curve by Stand") + theme(plot.title =element_text(hjust = 0.5)) + xlab("Subsampled Reads") + ylab("Observed ASVs") + theme(plot.title = element_text(size = 20),legend.key.size = unit(.8,"cm"),legend.text = element_text(size = 16), legend.title=element_text(size = 18), axis.text = element_text(size = 10), axis.title = element_text(size = 16)) + guides(color=guide_legend(title = "Stand"))
