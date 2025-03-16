
## 1. calculate gene FPKM ####
library(GenomicFeatures)
library(stringr)

# load read count matrix
cts = read.delim("PolysomeSeq.raw_counts.txt")
colnames(cts) <- str_replace(colnames(cts), "u5.3u5.4", "u5")


#读入统计counts的gtf
txdb <- makeTxDbFromGFF("Araport11.gtf",format="gtf")
columns(txdb)
#统计基因的exons
exons.list.per.gene <- exonsBy(txdb,by="gene")
#calculate effective length for each gene
exonic.gene.sizes <- as.data.frame(sum(width(reduce(exons.list.per.gene))))
exonic.gene.sizes$Gene.ID <- rownames(exonic.gene.sizes)
colnames(exonic.gene.sizes)[1] <- "len"
exonic.gene.sizes <- exonic.gene.sizes[,c(2,1)]

#merge read count and effective length for each gene
m <- merge(cts, exonic.gene.sizes, by="Gene.ID")

countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

countToEffCounts <- function(counts, len, effLen)
{
  counts * (len / effLen)
}

#calculate FPKM
m$FPKM.WT.input_1 <- with(m, countToFpkm(WT_input_1,len))
m$FPKM.WT.input_2 <- with(m, countToFpkm(WT_input_2,len))
m$FPKM.WT.input_3 <- with(m, countToFpkm(WT_input_3,len))
m$FPKM.WT.polysome_1 <- with(m, countToFpkm(WT_ribo_1,len))
m$FPKM.WT.polysome_2 <- with(m, countToFpkm(WT_ribo_2,len))
m$FPKM.WT.polysome_3 <- with(m, countToFpkm(WT_ribo_3,len))

m$FPKM.u5.input_1 <- with(m, countToFpkm(u5_input_1,len))
m$FPKM.u5.input_2 <- with(m, countToFpkm(u5_input_2,len))
m$FPKM.u5.input_3 <- with(m, countToFpkm(u5_input_3,len))
m$FPKM.u5.polysome_1 <- with(m, countToFpkm(u5_ribo_1,len))
m$FPKM.u5.polysome_2 <- with(m, countToFpkm(u5_ribo_2,len))
m$FPKM.u5.polysome_3 <- with(m, countToFpkm(u5_ribo_3,len))

m.fpkm <- m[,c("Gene.ID",colnames(m)[str_detect(colnames(m), pattern="FPKM.")])]

#write.table(m.fpkm, "PolysomeSeq.FPKM_matrix.txt", sep = "\t", quote = F, row.names = F)

## 2. calculate TE: FPKM(polysome)/FPKM(input) ####

# Translation efficiency (TE) was derived by normalizing polysome-associated RNA FPKM to cytoplasmic RNA FPKM.
library(tidyr)
library(dplyr)
library(openxlsx)

m.fpkm <- read.delim("PolysomeSeq.FPKM_matrix.txt")

# A pseudocount of 1e-5 was added to FPKM values to avoid computational failure of translation efficiency (TE) due to zero FPKM.
m.fpkm$TE.WT_1 <- (m.fpkm$FPKM.WT.polysome_1 + 1e-5)/(m.fpkm$FPKM.WT.input_1 + 1e-5)
m.fpkm$TE.WT_2 <- (m.fpkm$FPKM.WT.polysome_2 + 1e-5)/(m.fpkm$FPKM.WT.input_2 + 1e-5)
m.fpkm$TE.WT_3 <- (m.fpkm$FPKM.WT.polysome_3 + 1e-5)/(m.fpkm$FPKM.WT.input_3 + 1e-5)

m.fpkm$TE.u5_1 <- (m.fpkm$FPKM.u5.polysome_1 + 1e-5)/(m.fpkm$FPKM.u5.input_1 + 1e-5)
m.fpkm$TE.u5_2 <- (m.fpkm$FPKM.u5.polysome_2 + 1e-5)/(m.fpkm$FPKM.u5.input_2 + 1e-5)
m.fpkm$TE.u5_3 <- (m.fpkm$FPKM.u5.polysome_3 + 1e-5)/(m.fpkm$FPKM.u5.input_3 + 1e-5)

m.fpkm$Avg.TE.u5 <- (m.fpkm$TE.u5_1 + m.fpkm$TE.u5_2 + m.fpkm$TE.u5_3)/3
m.fpkm$Avg.TE.WT <- (m.fpkm$TE.WT_1 + m.fpkm$TE.WT_2 + m.fpkm$TE.WT_3)/3

m.fpkm$log2FC_TE.u5_WT <- log2(m.fpkm$Avg.TE.u5/m.fpkm$Avg.TE.WT)


FPKM.s.t.TE <- m.fpkm[,c("Gene.ID", str_subset(colnames(m.fpkm), "^TE"))]
colnames(FPKM.s.t.TE) <- str_remove(colnames(FPKM.s.t.TE), "TE.")

FPKM.s.t.TE.long <- pivot_longer(FPKM.s.t.TE, cols = 2:ncol(FPKM.s.t.TE), names_to = "sample",
                                 values_to = "TE")
FPKM.s.t.TE.long$material <- str_split_fixed(FPKM.s.t.TE.long$sample, "_",2)[,1]


library(rstatix)
FPKM.s.t.TE1 <- FPKM.s.t.TE[apply(FPKM.s.t.TE[,2:4], 1, var) !=0, ]
FPKM.s.t.TE1 <- FPKM.s.t.TE1[apply(FPKM.s.t.TE1[,5:7], 1, var) !=0, ]

FPKM.s.t.TE.u5 <- subset(FPKM.s.t.TE.long, Gene.ID %in% FPKM.s.t.TE1$Gene.ID 
                         & material %in% c("u5","WT"))
stat.u5.TE.oneside.less <- FPKM.s.t.TE.u5 %>% group_by(Gene.ID) %>%
  t_test(TE ~ material, alternative = "less")
stat.u5.TE.oneside.greater <- FPKM.s.t.TE.u5 %>% group_by(Gene.ID) %>%
  t_test(TE ~ material, alternative = "greater")

FPKM.s.t.stat1 <- left_join(m.fpkm, stat.u5.TE.oneside.less[,c("Gene.ID","p")], by="Gene.ID")
colnames(FPKM.s.t.stat1)[ncol(FPKM.s.t.stat1)] <- "p.oneside.less"
FPKM.s.t.stat1 <- left_join(FPKM.s.t.stat1, stat.u5.TE.oneside.greater[,c("Gene.ID","p")], by="Gene.ID")
colnames(FPKM.s.t.stat1)[ncol(FPKM.s.t.stat1)] <- "p.oneside.greater"
FPKM.s.t.stat1$pvalue[FPKM.s.t.stat1$Avg.TE.u5 < FPKM.s.t.stat1$Avg.TE.WT] <- FPKM.s.t.stat1$p.oneside.less[FPKM.s.t.stat1$Avg.TE.u5 < FPKM.s.t.stat1$Avg.TE.WT] 
FPKM.s.t.stat1$pvalue[FPKM.s.t.stat1$Avg.TE.u5 > FPKM.s.t.stat1$Avg.TE.WT] <- FPKM.s.t.stat1$p.oneside.greater[FPKM.s.t.stat1$Avg.TE.u5 > FPKM.s.t.stat1$Avg.TE.WT] 

ga <- read.xlsx("Araport11_annotation.xlsx")
#ga <- ga[,c(1,7:10)]
colnames(ga)[1] <- "Gene.ID"
ga$computational_description <- NULL
FPKM.s.TE <- left_join(FPKM.s.t.stat1, ga, by="Gene.ID")

#write.xlsx(FPKM.s.TE, "PolysomeSeq_TE_info.xlsx")
