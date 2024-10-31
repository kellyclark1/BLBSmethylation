library(BiocManager)
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)
library(RColorBrewer)
library(DMRcate)
library(DMRcatedata)
library(missMethyl)
library(stringr)
library(methylclock)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(viridis)
library(rcartocolor)

up.color = "#CC6677"
down.color = "#88CCEE"
control.color = "#888888"
blbs.color = "#117733"
tail.color = "#332288"
core.color = "#DDCC77"
a.color = "#AA4499"
b.color = "#44AA99"

setwd("~/Documents/DNA_methylation_histone/IDATS")

targets <- read.metharray.sheet(".", "IEB_BLBS_x28_5-14-24.csv")
targets

targets$Sample_Group <- c("control", "control", "control", "control", 
                          "control", "control", "control", "control", 
                          "control", "control", "control", "control", 
                          "histone", "histone", "histone", "histone", 
                          "histone", "histone", "histone", "histone", 
                          "histone", "histone", "histone", "histone", 
                          "histone", "histone", "histone", "histone")

targets$age <- c(11, 11, 2, 2, 1, 1, 13, 13, 17, 17, 11, 11,
                 4, 4, 10, 10, 5, 5, 2, 2, 2, 2, 8, 8, 15, 15, 18, 18)



rgSet <- read.metharray.exp(targets=targets)
rgSet

targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID
rgSet

detP <- detectionP(rgSet)
nprobesinit <- dim(detectionP(rgSet))

head(detP)


qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group, 
         pdf="qcReport.pdf")

keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
rgSet

targets <- targets[keep,]
targets[,1:5]

detP <- detP[,keep]
dim(detP)

mSetSq <- preprocessQuantile(rgSet) 

mSetRaw <- preprocessRaw(rgSet)

detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)

mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

annEPIC <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
keep <- !(featureNames(mSetSqFlt) %in% annEPIC$Name[annEPIC$chr %in% 
                                                      c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]

mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt

mVals <- as.data.frame(getM(mSetSqFlt))
head(mVals[,1:5])

bVals <- as.data.frame(getBeta(mSetSqFlt))
head(bVals[,1:5])

mVals$control.ctrlI <- (mVals$`control.ctrl I rep1` + mVals$`control.ctrl I rep2`) / 2
mVals$control.ctrlF <- (mVals$`control.ctrl F rep1` + mVals$`control.ctrl F rep2`) / 2
mVals$control.ctrlB <- (mVals$`control.ctrl B rep1` + mVals$`control.ctrl B rep2`) / 2
mVals$control.ctrlG <- (mVals$`control.ctrl G rep1` + mVals$`control.ctrl G rep2`) / 2
mVals$control.ctrlN <- (mVals$`control.ctrl N rep1` + mVals$`control.ctrl N rep2`) / 2
mVals$control.ctrlD <- (mVals$`control.ctrl D rep1` + mVals$`control.ctrl D rep2`) / 2
mVals$histone.EB257 <- (mVals$`histone.EB257 rep1` + mVals$`histone.EB257 rep2`) / 2
mVals$histone.EB252 <- (mVals$`histone.EB252 rep1` + mVals$`histone.EB252 rep2`) / 2
mVals$histone.EB249 <- (mVals$`histone.EB249 rep1` + mVals$`histone.EB249 rep2`) / 2
mVals$histone.EB452 <- (mVals$`histone.EB452 rep1` + mVals$`histone.EB452 rep2`) / 2
mVals$histone.EB427 <- (mVals$`histone.EB427 rep1` + mVals$`histone.EB427 rep2`) / 2
mVals$histone.EB256 <- (mVals$`histone.EB256 rep1` + mVals$`histone.EB256 rep2`) / 2
mVals$histone.EB318 <- (mVals$`histone.EB318 rep1` + mVals$`histone.EB318 rep2`) / 2
mVals$histone.EB383 <- (mVals$`histone.EB383 rep1` + mVals$`histone.EB383 rep2`) / 2



bVals$control.ctrlI <- (bVals$`control.ctrl I rep1` + bVals$`control.ctrl I rep2`) / 2
bVals$control.ctrlF <- (bVals$`control.ctrl F rep1` + bVals$`control.ctrl F rep2`) / 2
bVals$control.ctrlB <- (bVals$`control.ctrl B rep1` + bVals$`control.ctrl B rep2`) / 2
bVals$control.ctrlG <- (bVals$`control.ctrl G rep1` + bVals$`control.ctrl G rep2`) / 2
bVals$control.ctrlN <- (bVals$`control.ctrl N rep1` + bVals$`control.ctrl N rep2`) / 2
bVals$control.ctrlD <- (bVals$`control.ctrl D rep1` + bVals$`control.ctrl D rep2`) / 2
bVals$histone.EB257 <- (bVals$`histone.EB257 rep1` + bVals$`histone.EB257 rep2`) / 2
bVals$histone.EB252 <- (bVals$`histone.EB252 rep1` + bVals$`histone.EB252 rep2`) / 2
bVals$histone.EB249 <- (bVals$`histone.EB249 rep1` + bVals$`histone.EB249 rep2`) / 2
bVals$histone.EB452 <- (bVals$`histone.EB452 rep1` + bVals$`histone.EB452 rep2`) / 2
bVals$histone.EB427 <- (bVals$`histone.EB427 rep1` + bVals$`histone.EB427 rep2`) / 2
bVals$histone.EB256 <- (bVals$`histone.EB256 rep1` + bVals$`histone.EB256 rep2`) / 2
bVals$histone.EB318 <- (bVals$`histone.EB318 rep1` + bVals$`histone.EB318 rep2`) / 2
bVals$histone.EB383 <- (bVals$`histone.EB383 rep1` + bVals$`histone.EB383 rep2`) / 2


mVals <- as.matrix(mVals %>% dplyr::select(29:42))
bVals <- as.matrix(bVals %>% dplyr::select(29:42))


targets <- targets[c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27),]
targets$Sample_Name <- substr(targets$Sample_Name,1,nchar(targets$Sample_Name)-5)
targets$ID <- substr(targets$ID,1,nchar(targets$ID)-5)

targets$individual <- substr(targets$Sample_Name,1,nchar(targets$Sample_Name)-5)
condition <- factor(targets$Sample_Group)

design <- model.matrix(~0+condition, data=targets)
colnames(design) <- c(levels(condition))

fit <- lmFit(mVals, design)
contMatrix <- makeContrasts(histone-control,
                            levels=design)
contMatrix

fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

annEPICSub <- annEPIC[match(rownames(mVals),annEPIC$Name),
                      c(1:4,12:19,24:ncol(annEPIC))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=annEPICSub)
head(DMPs)


mVals <- rmSNPandCH(mVals)
mVals <- rmPosReps(mVals, filter.strategy = "mean")
myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "histone - control", arraytype = "EPICv2")

str(myAnnotation)


long.bvals2 = tidyr::pivot_longer(bVals.df, cols = -CpG_name) %>% dplyr::rename(sample = "name")
long.bvals2$Group = "BLBS"
long.bvals2$Group[substr(long.bvals$sample, 1,7) == "control"] <- "control"


pdf(height = 4, width = 6, "Bvalue_histogram.pdf")
ggplot(long.bvals2, aes(value, fill = Group, color = Group)) +
  geom_density(data = long.bvals2 %>% filter(Group == "BLBS"), alpha = 0.1,
               linewidth = 1, aes(y = ..density..)) +
  geom_density(data = long.bvals2 %>% filter(Group == "control"), alpha = 0.1,
               linewidth = 1, aes(y = ..density..)) +
  theme_bw() +
  scale_fill_manual(values = c(blbs.color, control.color)) +
  scale_color_manual(values = c(blbs.color, control.color)) +
  labs(title = "Beta values", y = "Density", x = "Value") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()






ggplot(DMPs, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(alpha = 0.2) +
  geom_hline(aes(yintercept = -log10(0.05)), color = "red") +
  labs(title = "CpG probe methylation in BLBS patient-derived\nvs control dermal fibroblasts",
       x = "log(fold change)",
       y = "-log10(FDR)") +
  theme_bw() +
  xlim(-4.3, 4.3) +
  theme(plot.title = element_text(hjust = 0.5))

top_DMPs <- DMPs %>% arrange(P.Value) %>% head(1000)

bVals.df <- as.data.frame(bVals)

bVals.filt <- as.matrix(bVals.df %>% dplyr::filter(rownames(bVals.df) %in% top_DMPs$Name))
colnames(bVals.filt) <- c("control I", "control F", "control B",
                          "control G", "control N",
                          "control D", "H3-3A p.T45I",
                            "H3-3A p.R17G",  "H3-3B p.G34V", 
                            "H3-3A p.T6S",  "H3-3A p.S31F",  
                            "H3-3A p.G90R",   "H3-3B p.V117V", 
                              "H3-3B p.P121R") 


m2 <- mVals
colnames(m2) <- c("control I", "control F", "control B", "control G", 
                  "control N", "control D",
                  "H3-3A p.T45I",   "H3-3A p.R17G",
                  "H3-3B p.G34V",   "H3-3A p.T6S",  
                  "H3-3A p.S31P",  "H3-3A p.G90R",  "H3-3B p.V117V",
                  "H3-3B p.P121R")
par(1,1)
p12 <- plotMDS(m2, top=1000, gene.selection="common",
               col=c(rep(control.color, 6), rep(blbs.color, 8)), dim=c(1,2))



p12 <- as.data.frame(p12) %>% dplyr::select(x, y)

targets$Sample_Group[targets$Sample_Group == "control"] <- "Control"
targets$Sample_Group[targets$Sample_Group == "histone"] <- "BLBS"

ggplot(p12, aes(x, y, color = targets$Sample_Group)) + geom_point(aes(shape = targets$Sample_Group), size = 3) + theme_bw() +
  geom_text_repel(aes(label = rownames(p12)), nudge_y = 0.03, show.legend = F) +
  scale_color_manual(values = c(blbs.color, control.color)) + scale_shape_manual(values = c(15,19)) +
  labs(x = "PC1 (27%)",
       y = "PC2 (22%)",
       title = "PCA of BLBS patient-derived and control\ndermal fibroblasts based on filtered M-values",
       color = "Group", shape = "Group") +
  #geom_polygon(stat = "ellipse", aes(color = targets$Sample_Group), fill = NA) + scale_color_manual(values=c(blbs.color, control.color)) +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(p12, aes(x, y, color = targets$age)) + geom_point(aes(shape = targets$Gender), size = 3) + theme_bw() +
  geom_text_repel(aes(label = rownames(p12)), nudge_y = 0.03, show.legend = F) +
  scale_color_gradient(low = "grey", high = "black") +
  scale_shape_manual(values = c(15,19)) +
  labs(x = "PC1 (27%)",
       y = "PC2 (22%)",
       title = "PCA of BLBS patient-derived and control\ndermal fibroblasts based on filtered M-values",
       color = "Age", shape = "Sex") +
  theme(plot.title = element_text(hjust = 0.5))



bVals.df <- as.data.frame(bVals)
bVals.df$CpG_name <- rownames(bVals.df)
long.bvals2 = tidyr::pivot_longer(bVals.df, cols = -CpG_name) %>% dplyr::rename(sample = "name")
long.bvals2$Group = "Control"
long.bvals2$Group[substr(long.bvals2$sample, 1,7) == "histone"] <- "BLBS"



ggplot(long.bvals2, aes(value, fill = Group, color = Group)) +
  geom_density(data = long.bvals2 %>% dplyr::filter(Group == "Control"), alpha = 0.1,
               linewidth = 1, aes(y = ..density..)) +
  geom_density(data = long.bvals2 %>% dplyr::filter(Group == "BLBS"), alpha = 0.1,
               linewidth = 1, aes(y = ..density..)) +
  theme_bw() +
  scale_fill_manual(values = c(blbs.color, control.color)) +
  scale_color_manual(values = c(blbs.color, control.color)) +
  labs(title = "Beta values in BLBS patient-derived and\ncontrol dermal fibroblasts", y = "Density", x = "Value") +
  theme(plot.title = element_text(hjust = 0.5))





