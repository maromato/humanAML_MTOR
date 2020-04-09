
library(pheatmap)
library(RColorBrewer)
library(edgeR)
library(maptools)
library(ggplot2)
library(gplots)
library(made4)
library(calibrate)
library(fgsea)
library(survival)
library(survminer)
library(dplyr)
library(data.table)
library(genefilter)
library(reshape2)
library(dendextend)

rm(list = ls())
getwd()
setwd("/Users/USERNAME/YOUR PATH/YOUR FOLDER")

############## set up the data set############
# read trascriptome data 
x <- read.table("./GSE12417_mod_4.txt", row.names=1,header=T, sep="\t") #old data file with new annotation

dim(x)
head(x)
y <- NULL
y <- x[,3:ncol(x)]
y <- as.matrix(y)
rownames(y) <- x[,2]
head(y)
dim(y)

# read meatadata
s <- read.table("GSE12417-GPL96_meta.txt",header=T, sep="\t")   
head(s)

# read gene set data
b = read.table("RAPAMYCIN_SENSITIVE_GENES_mod.txt",sep="\t",header=F)
c = read.table("HALLMARK_MTORC1_SIGNALING_mod.txt",sep="\t",header=F)
d = read.table("DEG_venus_high_low.txt",sep="\t",header=F)
e = read.table("CUNNINGHAM_RAPAMYCIN_DN.txt",sep="\t",header=F)
f = read.table("DEG_NT_TX_1.txt",sep="\t",header=F)
g = read.table("DEG_venus_low_high_all.txt",sep="\t",header=F)
h = read.table("Myc_core.txt",sep="\t",header=F)
i = read.table("YUAN_MTOR_SIG.txt",sep="\t",header=F)


################ Heatmap and survival curve for mouse DEG ############################################
z3 <- y[rownames(y) %in% toupper(d[,1]),]


heatplot(as.matrix(z3),labRow = FALSE, cexCol=0.4)

# extract col dendrogram information
sTree<-heatplot(as.matrix(z3), returnSampleTree=TRUE)
plot(sTree)

# extraxt clustering information
col_clust <- as.hclust(sTree)
mTORhigh <- cutree(col_clust,k=3,order_clusters_as_data = FALSE)

dend <-as.dendrogram(col_clust)
clust.cutree <- dendextend:::cutree(dend, k=3, order_clusters_as_data = FALSE)
idx <- order(as.numeric(names(clust.cutree)))
clust.cutree <- clust.cutree[idx]
dend1 <- color_branches(dend, k = 3, groupLabels = TRUE )
plot(dend1)

# pick up mTOR high populaiton
mTORhigh <- as.matrix(mTORhigh)
mTORhigh

MTH <- NULL
MTH<- mTORhigh == 1
MTH
head(MTH)
MTH <- cbind(rownames(MTH),MTH)
colnames(MTH) <- c("ID","MTOR_DEG")
head(MTH)

# merging mTOR high and low group information to meta data
t <- NULL
t <- merge(s,MTH, by="ID")
dim(t)
head(t)
attach(t)

#survival curve
fit <- survfit(Surv(OS, Status) ~ MTOR_DEG, data = t,type="kaplan-meier")
ggsurvplot(fit, data= t, legend.labs = c("DGE low", "DGE high"),
           pval = TRUE)

pchisq(survdiff(Surv(OS, Status)~MTOR_DEG, data=t)$chisq, 1, lower.tail=FALSE)


#p value calculation
p_all <- pchisq(survdiff(Surv(OS, Status)~MTOR_DEG, data=t)$chisq, 1, lower.tail=FALSE)# overall

st = 300
for (i in 1:length(t$Status)) {
  if (t$OS[i] >= st) {
    t$OS300[i] <- st
    t$Status300[i] <- 0
  }
  else {
    t$OS300[i] <- t$OS[i]
    t$Status300[i] <- t$Status[i]
  }
}
p_300 <- pchisq(survdiff(Surv(OS300, Status300)~MTOR_DEG, data=t)$chisq, 1, lower.tail=FALSE)# pval at day300

st = 600
for (i in 1:length(t$Status)) {
  if (t$OS[i] >= st) {
    t$OS600[i] <- st
    t$Status600[i] <- 0
  }
  else {
    t$OS600[i] <- t$OS[i]
    t$Status600[i] <- t$Status[i]
  }
}

p_600 <- pchisq(survdiff(Surv(OS600, Status600)~MTOR_DEG, data=t)$chisq, 1, lower.tail=FALSE) # pval at day600

st = 900
for (i in 1:length(t$Status)) {
  if (t$OS[i] >= st) {
    t$OS900[i] <- st
    t$Status900[i] <- 0
  }
  else {
    t$OS900[i] <- t$OS[i]
    t$Status900[i] <- t$Status[i]
  }
}

p_900 <- pchisq(survdiff(Surv(OS900, Status900)~MTOR_DEG, data=t)$chisq, 1, lower.tail=FALSE) #pvalue at day 900

p_total <- rbind(p_all, p_300, p_600, p_900)
colnames(p_total) <- "MTOR_DEG"
p_MTOR_DEG <- p_total
p_total


write.table(t, "G12417_metadata_MTOR_DEG.txt",sep="\t")
write.table(p_total, "G12417_pvalue_MTOR_DEG.txt",sep="\t")

##################  DEG of DEG #######################
x2 <- x[,3:ncol(x)]
head(x)
head(t)


t[t$MTOR_DEG == TRUE,1]


length(t[t$MTOR_DEG == TRUE,1])
MTORhigh <- NULL
MTORimlow <- NULL
MTORhigh<- x2[,colnames(x2) %in% t[t$MTOR_DEG == TRUE,1]]
MTORimlow <- x2[,! colnames(x2) %in% colnames(MTORhigh) ]
dim(MTORimlow)
dim(MTORhigh)

hMTOR <- cbind(MTORhigh, MTORimlow)
ncol(MTORhigh)
torgroup <- NULL
torgroup <- c(rep(1,ncol(MTORhigh)),rep(2,ncol(MTORimlow)))
length(torgroup)
torgroup
head(hMTOR)
hMTOR <- as.matrix(hMTOR)

out <- rowttests(hMTOR,factor(torgroup))
q.value <- p.adjust(out$p.value,method="BH")
out2 <- cbind(out, q.value)
head(out2)

out2 <- cbind(out2,x[,2])

head(out2)
out2 <- out2[order(out2$p.value),]

colnames(out2)[ncol(out2)] <- "Gene"
head(out2)
write.table(out2,"humanAML_MTOR_DEG.txt",sep="\t")
head(out2)

out2 <- out2[order(out2$p.value),]
head(out2,10)
dim(out2)
out3 <- out2[out2$p.value < 0.05,]
dim(out3)
head(out3)
write.table(out3,"humanAML_MTOR_DEG_p_05.txt",sep="\t")

tail(out2)
out4 <- out2[(out2$q.value < 0.05),]
out4 <- out4[out4$dm >0,]
dim(out4)
intersect(toupper(e[,1]),out4$Gene)

write.table(intersect(toupper(e[,1]),out4$Gene),"human_mouse_MTOR_DEG_high_DEG.txt",sep="\t")

########## GSEA for mTORC high and low of MTOR_DEG ###########

target <- intersect(toupper(e[,1]),out4$Gene)

MTOR_H  <- y[rownames(y) %in% target,colnames(MTORhigh)]  
MTOR_L  <- y[rownames(y) %in% target,colnames(MTORimlow)]  
dim(MTOR_H)
dim(MTOR_L)
rowMeans(MTOR_H)
FC <- NULL
FC <- log(rowMeans(MTOR_H)/rowMeans(MTOR_L))
names(FC) <- rownames(y[rownames(y) %in% target,])
head(FC)
length(FC)
FC[order(names(FC), -abs(FC) )] ### sort first
FC[ !duplicated(names(FC)) ]  ### Keep highest
length(FC)

#ep <-gmtPathways("Hallmarks.gmt")
ep <- gmtPathways("Oki-tools-ver5.gmt")
#ep <-gmtPathways("Oki-tools-ver6.gmt")

er = FC
fgseaRes <- fgsea(pathways = ep, 
                  stats = er,
                  minSize= 20,
                  maxSize=500,
                  nperm=10000)

head(fgseaRes[order(pval), ],30)

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

dev.off()

#topPathways <- fgseaRes[head(order(pval), n=20), pathway]

plotGseaTable(ep[topPathways], er, fgseaRes, 
              gseaParam = 0.5)



write.table(unlist(fgseaRes[order(pval), ]),"GSEA_mouse_human_MTOR_DEG.txt",sep="\t")
