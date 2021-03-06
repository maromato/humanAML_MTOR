
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

R.version
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
#d = read.table("DEG_venus_high_low.txt",sep="\t",header=F)
e = read.table("CUNNINGHAM_RAPAMYCIN_DN.txt",sep="\t",header=F)


################ Heatmap with CUNNINGHAM_RAPAMYCIN_DN.txt #######################
z4 <- y[rownames(y) %in% e[,1],]
head(z4)

heatplot(as.matrix(z4),labRow = FALSE, cexCol=0.4)

# extract col dendrogram information
sTree<-heatplot(as.matrix(z4), returnSampleTree=TRUE)
plot(sTree,labels=NULL)

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

colnames(MTH) <- c("ID","MTOR_CUNNING")
head(MTH)

# merging mTOR high and low group information to meta data
t <- NULL
t <- merge(s,MTH, by="ID")
dim(t)
head(t)
tail(t)

#survival curve
fit <- survfit(Surv(OS, Status) ~ MTOR_CUNNING, data = t,type="kaplan-meier")
ggsurvplot(fit, data= t, legend.labs = c("mTORC low", "mTORC high"),
           pval = TRUE)

p_all <- pchisq(survdiff(Surv(OS, Status)~MTOR_CUNNING, data=t)$chisq, 1, lower.tail=FALSE)

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
p_300 <- pchisq(survdiff(Surv(OS300, Status300)~MTOR_CUNNING, data=t)$chisq, 1, lower.tail=FALSE)

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

p_600 <- pchisq(survdiff(Surv(OS600, Status600)~MTOR_CUNNING, data=t)$chisq, 1, lower.tail=FALSE)

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

p_900 <- pchisq(survdiff(Surv(OS900, Status900)~MTOR_CUNNING, data=t)$chisq, 1, lower.tail=FALSE)

p_total <- rbind(p_all, p_300, p_600, p_900)
colnames(p_total) <- "MTOR_CUNNING"
p_MTOR_CUNNING <- p_total
p_total

write.table(t, "G12417_metadata_MTOR_CUNNING.txt",sep="\t")
write.table(p_total, "G12417_pvalue_MTOR_CUNNING.txt",sep="\t")



################# DEG of MTOR_CUNNING ######################
x2 <- x[,3:ncol(x)]
head(x)
head(t)


t[t$MTOR_CUNNING == TRUE,1]


length(t[t$MTOR_CUNNING == TRUE,1])
MTORhigh <- NULL
MTORimlow <- NULL
MTORhigh<- x2[,colnames(x2) %in% t[t$MTOR_CUNNING == TRUE,1]]
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
write.table(out2,"humanAML_MTOR_CUNNING.txt",sep="\t")
head(out2)

out2 <- out2[order(out2$p.value),]
head(out2,10)
dim(out2)
out3 <- out2[out2$p.value < 0.05,]
dim(out3)
head(out3)
write.table(out3,"humanAML_MTOR_CUNNING_p_05.txt",sep="\t")

tail(out2)
out4 <- out2[(out2$q.value < 0.05),]
out4 <- out4[out4$dm >0,]
dim(out4)
intersect(toupper(e[,1]),out4$Gene)

write.table(intersect(toupper(e[,1]),out4$Gene),"human_mouse_MTOR_CUNNING_high_DEG.txt",sep="\t")


