##Comammox Figure 3 creation code
#Using all plots for AOA and AOB for rank order assessment
#Load data

qpcr_all <- read.csv("qpcr_aoa_aob_caob_ph.csv", header=TRUE)

##Geometric means
#To assess stands that had consistent results across all plots within that stand, we used geometric means

library(tidyverse)
library(psych)

geomeans <- qpcr_all[,c(2,5)] %>% group_by(STAND) %>% summarize(across(everything(), ~geometric.mean(.x)))

#Four sites with a geometric mean
#6, 7, 22, 24
comammox_gmsub <- subset(qpcr_all[,c(2,5:6)], qpcr_all$STAND %in% geomeans$STAND[geomeans$CAOB > 0])

#soil pH for copies comammox Figure 3
library(ggplot2)

scientific_10 <- function(x) {
  gsub("e", " x 10", scientific_format()(x))
}

ph.plot.quad <- ggplot(data=comammox_gmsub, aes(x=pH, y=CAOB)) + geom_point(col="#454545") + geom_smooth(aes(y=CAOB), method="lm", formula = y ~ x + I(x^2)) + theme_classic() + ylab ("Gene Copy / g soil") + theme(
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16)
) + scale_y_continuous(label=scientific_10)


ph.summary <- summary(lm(CAOB~pH +I(pH^2), data=comammox_gmsub))
Adj.r.square.ph <- round(ph.summary[9]$adj.r.squared,3)
pvalue.ph <- round(0.1568, 2) #Reported in the summary output

r2.pvalue.ph <- data.frame(label = c(paste("italic(R^2) == ", Adj.r.square.ph), ifelse(pvalue.ph>0.001, paste("italic(P) == ", pvalue.ph), paste("italic(P) < ", 0.001))))

Figure_3 <- ph.plot.quad + annotate(geom="text", x=4.3 ,y=c(4e+09,3.68e+09),size=4, label=r2.pvalue.ph$label, parse=TRUE)
Figure_3
ggsave("Figure3r.pdf",height=10, width=8.5, dpi=600)

