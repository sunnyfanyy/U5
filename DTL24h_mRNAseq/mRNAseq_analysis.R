library(GenomicFeatures)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggtext)

## 1. calculate gene FPKM ####
library(GenomicFeatures)

# load read count matrix
cts = read.delim("D_DTL12h_DTL24h_DTL48h.raw_counts.txt")

colnames(cts)
cts.m <- cts
cts.m$D.WT <- cts.m$Dark_WT_1 + cts.m$Dark_WT_2 + cts.m$Dark_WT_3
cts.m$DTL12h.WT <- cts.m$DTL12h_WT_1 + cts.m$DTL12h_WT_2 + cts.m$DTL12h_WT_3
cts.m$DTL24h.WT <- cts.m$DTL24h_WT_1 + cts.m$DTL24h_WT_2 + cts.m$DTL24h_WT_3
cts.m$DTL48h.WT <- cts.m$DTL48h_WT_1 + cts.m$DTL48h_WT_2 + cts.m$DTL48h_WT_3
cts.m$DTL24h.u5 <- cts.m$DTL24h_u5.3u5.4_1 + cts.m$DTL24h_u5.3u5.4_2 + cts.m$DTL24h_u5.3u5.4_3

cts.m <- cts.m[,c("D.WT","DTL12h.WT","DTL24h.WT","DTL48h.WT","DTL24h.u5")]

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
m <- merge(cts.m, exonic.gene.sizes, by="Gene.ID")

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
m$FPKM.D.WT <- with(m, countToFpkm(D.WT,len))
m$FPKM.DTL12h.WT <- with(m, countToFpkm(DTL12h.WT,len))
m$FPKM.DTL24h.WT <- with(m, countToFpkm(DTL24h.WT,len))
m$FPKM.DTL48h.WT <- with(m, countToFpkm(DTL48h.WT,len))
m$FPKM.DTL24h.u5 <- with(m, countToFpkm(DTL24h.u5,len))

m.fpkm <- m[,c("Gene.ID",colnames(m)[str_detect(colnames(m), pattern="FPKM.")])]

#write.csv(m.fpkm, "D_DTL12h_DTL24h_DTL48h.calculated_gene_FPKM.csv", row.names = F)

## 2. DESeq2 analysis ####
library(DESeq2)
library(openxlsx)
# load read count matrix
cts = read.delim("D_DTL12h_DTL24h_DTL48h.raw_counts.txt", row.names = "Gene.ID")

# delete non-expressed genes
cts <- cts[apply(cts[,1:ncol(cts)]>0, 1, any),]   

#构建每个sample所对应的条件
light <- rep(c("D","DTL12h","DTL24h","DTL48h","DTL24h"), each=3)
material <- c(rep("WT", times=12), rep("u5", times=3))

coldata <- data.frame(row.names = colnames(cts),
                      light=light, material=material)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ material)
ddsMF <- dds
ddsMF$group = factor(paste0(ddsMF$light,".",ddsMF$material))
design(ddsMF) = ~0+group
ddsMF <- DESeq(ddsMF)

#Comparing the differences between WT and u5-3 u5-4 double mutants under DTL24h
res.DTL24h.u5_WT <- results(ddsMF, contrast = c("group","DTL24h.u5","DTL24h.WT")) #log2FC为DTL24h.cn1cn2 - DTL24h.WT

re.u5 <- res.DTL24h.u5_WT[,c("log2FoldChange","pvalue","padj")]
colnames(re.u5) <- c("DTL24h.u5_WT.FC","DTL24h.u5_WT.p","DTL24h.u5_WT.q")

#Comparing the differences between WT in Dark and DTL12h
res.WT.DTL12h_D <- results(ddsMF, contrast = c("group","DTL12h.WT","D.WT")) #log2FC为DTL12h.WT - D.WT
#Comparing the differences between WT in Dark and DTL24h
res.WT.DTL24h_D <- results(ddsMF, contrast = c("group","DTL24h.WT","D.WT")) #log2FC为DTL24h.WT - D.WT
#Comparing the differences between WT in Dark and DTL48h
res.WT.DTL48h_D <- results(ddsMF, contrast = c("group","DTL48h.WT","D.WT")) 

#Comparing the differences between WT in DTL12h and DTL24h
res.WT.DTL24h_DTL12h <- results(ddsMF, contrast = c("group","DTL24h.WT","DTL12h.WT")) #log2FC为DTL24h.WT - DTL12h.WT
#Comparing the differences between WT in DTL12h and DTL48h
res.WT.DTL48h_DTL12h <- results(ddsMF, contrast = c("group","DTL48h.WT","DTL12h.WT")) #log2FC为DTL48h.WT - DTL12h.WT
#Comparing the differences between WT in DTL24h and DTL48h
res.WT.DTL48h_DTL24h <- results(ddsMF, contrast = c("group","DTL48h.WT","DTL24h.WT")) #log2FC为DTL48h.WT - DTL24h.WT


re.DTL <- cbind(res.WT.DTL12h_D[,c("log2FoldChange","padj")],
              res.WT.DTL24h_D[,c("log2FoldChange","padj")],
              res.WT.DTL48h_D[,c("log2FoldChange","padj")],
              res.WT.DTL24h_DTL12h[,c("log2FoldChange","padj")],
              res.WT.DTL48h_DTL12h[,c("log2FoldChange","padj")],
              res.WT.DTL48h_DTL24h[,c("log2FoldChange","padj")])
colnames(re.DTL) <- c("WT.DTL12h_D.FC","WT.DTL12h_D.q",
                    "WT.DTL24h_D.FC","WT.DTL24h_D.q",
                    "WT.DTL48h_D.FC","WT.DTL48h_D.q",
                    "WT.DTL24h_DTL12h.FC","WT.DTL24h_DTL12h.q",
                    "WT.DTL48h_DTL12h.FC","WT.DTL48h_DTL12h.q",
                    "WT.DTL48h_DTL24h.FC","WT.DTL48h_DTL24h.q")


# read FPKM information
fpkm.DTL <- read.csv("D_DTL12h_DTL24h_DTL48h.calculated_gene_FPKM.csv")

# merge FPKM and DESeq2 results
m.DTL <- merge(fpkm.DTL[,c(1:5)], re.DTL, by="Gene.ID")
m.u5 <- merge(fpkm.DTL[,c("Gene.ID", "FPKM.DTL24h.WT", "FPKM.DTL24h.u5")],
              re.u5, by="Gene.ID")

# add gene annotation
gt <- read.csv("../Subcellular_RNAseq/merged_gene_annotation.csv")
m.DTL <- merge(m.DTL, gt, by="Gene.ID", all.x=T, sort=F)
m.u5 <- merge(m.u5, gt, by="Gene.ID", all.x=T, sort=F)

# save FPKM and DESeq2 results for WT under different light condition
#write.xlsx(m.DTL, "D_DTL12h_DTL24h_DTL48h.merged_info.xlsx")

# save FPKM and DESeq2 results for WT and u5 under DTL24h
#write.xlsx(m.u5, "DTL24h_WT_u5.merged_info.xlsx")

# get normalized reads
vsd <- vst(dds, blind=FALSE)
vsd <- assay(vsd) %>% as.data.frame()
vsd.keep <- c(7:9, 13:15)
#write.csv(vsd.keep, "DTL24h.VST_DESeq2_normalized_counts_sample.csv", quote = F)


## 3. Differentially expressed genes regulated by U5 ####
library(openxlsx)

m.u5 <- read.xlsx("DTL24h_WT_u5.merged_info.xlsx")

# DEGs definition: |FoldChange| > 1.5 & p < 0.05
u5.24h.DEGs <- subset(m.u5, abs(DTL24h.u5_WT.FC) > log2(1.5) & DTL24h.u5_WT.p < 0.05)
u5.24h.DEGs$regulated_in_u5[u5.24h.DEGs$DTL24h.u5_WT.FC > 0] <- "Up-regulated"
u5.24h.DEGs$regulated_in_u5[u5.24h.DEGs$DTL24h.u5_WT.FC < 0] <- "Down-regulated"
table(u5.24h.DEGs$regulated_in_u5)
# Down-regulated   Up-regulated 
# 480            460 
#write.xlsx(u5.24h.DEGs, "DTL24h_u5_vs._WT_DEGs_list.xlsx")

## 4. differential alternative splicing events (ASE) ####
# ASE definition: abs(IncLevelDifference) > 0.1 & FDR < 0.05

f <- list.files(path = "rMATS")
fa <- list.files(path = "rMATS", pattern = ".JC.txt",recursive = T, full.names = T)
# read gene annotation
ga <- read.csv("../Subcellular_RNAseq/merged_gene_annotation.csv")
for (i in f) {
  print(i)
  e <- dir(path = paste0("rMATS/",i), pattern = ".JC.txt", full.names = T)
  # print(e)
  e1 <- read.delim(e[1])
  e1$Event <- "A3SS"
  e2 <- read.delim(e[2])
  e2$Event <- "A5SS"
  e3 <- read.delim(e[3])
  e3$Event <- "MXE"
  e4 <- read.delim(e[4])
  e4$Event <- "RI"
  e5 <- read.delim(e[5])
  e5$Event <- "SE"
  
  e1 <- filter(e1, abs(IncLevelDifference) > 0.1 & FDR < 0.05)
  e2 <- filter(e2, abs(IncLevelDifference) > 0.1 & FDR < 0.05)
  e3 <- filter(e3, abs(IncLevelDifference) > 0.1 & FDR < 0.05)
  e4 <- filter(e4, abs(IncLevelDifference) > 0.1 & FDR < 0.05)
  e5 <- filter(e5, abs(IncLevelDifference) > 0.1 & FDR < 0.05)
  
  e1 <- merge(e1, ga, by.x="GeneID", by.y="ID", sort=F)
  e2 <- merge(e2, ga, by.x="GeneID", by.y="ID", sort=F)
  e3 <- merge(e3, ga, by.x="GeneID", by.y="ID", sort=F)
  e4 <- merge(e4, ga, by.x="GeneID", by.y="ID", sort=F)
  e5 <- merge(e5, ga, by.x="GeneID", by.y="ID", sort=F)
  
  wb <- createWorkbook() # 创建workbook
  addWorksheet(wb,"A3SS") # 给workbook增加worksheet
  addWorksheet(wb,"A5SS") # 再增加一个worksheet
  addWorksheet(wb,"MXE") # 再增加一个worksheet
  addWorksheet(wb,"RI") # 再增加一个worksheet
  addWorksheet(wb,"SE") # 再增加一个worksheet
  # 保存数据到相应的worksheet
  writeData(wb,"A3SS",e1) 
  writeData(wb,"A5SS",e2)
  writeData(wb,"MXE",e3)
  writeData(wb,"RI",e4)
  writeData(wb,"SE",e5)
  # 保存workbook到xlsx文件
  saveWorkbook(wb, paste(i,"_ASE.xlsx",sep = ""), overwrite = T)
}

## 5. ASE statistics for *smd3a* and *smd3b* mutant: Figure S8D,E ####
library(openxlsx)
library(ggplot2)
library(ggpubr)

smd3.A3SS <- read.xlsx("smd3a_smd3b_ASE_Type.xlsx", sheet = "A3SS") 
smd3.A5SS <- read.xlsx("smd3a_smd3b_ASE_Type.xlsx", sheet = "A5SS") 
smd3.RI <- read.xlsx("smd3a_smd3b_ASE_Type.xlsx", sheet = "RI") 
smd3.SE <- read.xlsx("smd3a_smd3b_ASE_Type.xlsx", sheet = "SE") 


#统计个数
smd3a.A3SS.num <- subset(smd3.A3SS, Event =="A3SS" & Type %in% c("SmD3a", "Co-regulated"))
smd3a.A5SS.num <- subset(smd3.A5SS, Event =="A5SS" & Type %in% c("SmD3a", "Co-regulated"))
smd3a.RI.num <- subset(smd3.RI, Event =="RI" & Type %in% c("SmD3a", "Co-regulated"))
smd3a.SE.num <- subset(smd3.SE, Event =="SE" & Type %in% c("SmD3a", "Co-regulated"))

smd3b.A3SS.num <- subset(smd3.A3SS, Event =="A3SS" & Type %in% c("SmD3b", "Co-regulated"))
smd3b.A5SS.num <- subset(smd3.A5SS, Event =="A5SS" & Type %in% c("SmD3b", "Co-regulated"))
smd3b.RI.num <- subset(smd3.RI, Event =="RI" & Type %in% c("SmD3b", "Co-regulated"))
smd3b.SE.num <- subset(smd3.SE, Event =="SE" & Type %in% c("SmD3b", "Co-regulated"))

smd3_ASE.mf <- data.frame(material=rep(c("smd3a","smd3b"),each=4),
                          event=rep(c("SE","RI","A5SS","A3SS"),times=2),
                          number=c(nrow(smd3a.SE.num),nrow(smd3a.RI.num),nrow(smd3a.A5SS.num),nrow(smd3a.A3SS.num),
                                   nrow(smd3b.SE.num),nrow(smd3b.RI.num),nrow(smd3b.A5SS.num),nrow(smd3b.A3SS.num)
                          )
)

smd3_ASE.mf$event <- factor(smd3_ASE.mf$event, levels = c("SE", "RI", "A5SS", "A3SS"))
ggplot(smd3_ASE.mf, aes(x=material, y=number, fill=event))+ xlab("") + ylab("ASE") +
  geom_col(color="black", width = 0.5) +  ggtitle("DTL-24") +
  scale_y_continuous(expand = c(0,0), limits = c(0,600), breaks = seq(0,600,200)) +
  scale_fill_manual(values = c("#df352e", "#186dae", "#fcb86b", "#f8e369")) +
  theme_classic() + scale_x_discrete(labels=c("smd3a/WT", "smd3b/WT")) +
  theme(axis.text = element_text(color = "black"),legend.title = element_blank(),
        legend.position= c(0.25, 0.75), plot.title = element_text(hjust = 0.5))
ggsave("FigureS8D.pdf", width = 6, height = 8, units = "cm")

# scatter plot of delta_PSI for smd3a and smd3b 

smd3_deltaPSI <- rbind(smd3.A3SS[,c("Event", "smd3a.deltaPSI", "smd3b.deltaPSI", "Type")],
                       smd3.A5SS[,c("Event", "smd3a.deltaPSI", "smd3b.deltaPSI", "Type")],
                       smd3.RI[,c("Event", "smd3a.deltaPSI", "smd3b.deltaPSI", "Type")],
                       smd3.SE[,c("Event", "smd3a.deltaPSI", "smd3b.deltaPSI", "Type")])

table(smd3_deltaPSI$Type)
# Co-regulated        SmD3a        SmD3b 
# 165           60          302 

  
# convert to percentage 
smd3_deltaPSI$smd3a.p <- smd3_deltaPSI$smd3a.deltaPSI * 100
smd3_deltaPSI$smd3b.p <- smd3_deltaPSI$smd3b.deltaPSI * 100

smd3a_only <- smd3_deltaPSI[smd3_deltaPSI$Type == 'SmD3a',]
smd3b_only <- smd3_deltaPSI[smd3_deltaPSI$Type == 'SmD3b',]
co_regulate <- smd3_deltaPSI[smd3_deltaPSI$Type == 'Co-regulated',]
smd3_deltaPSI$Type <- factor(smd3_deltaPSI$Type,
                             levels = c("SmD3a", "SmD3b", "Co-regulated"))

library(showtext)
showtext_auto(enable=TRUE)

ggplot(smd3_deltaPSI, aes(smd3b.p, smd3a.p)) +
  geom_point(aes(color = Type), size=2.5)+
  scale_size(range = c(0, 4)) +
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = 'black', fill = 'transparent')) +
  #stat_cor(data=smd3_deltaPSI, method = "spearman", label.x = -50, label.y = 50, color = 'black')+ 
  stat_cor(data=smd3a_only, method = "spearman", label.x = 20, label.y = -40, color = '#186dae')+
  stat_cor(data=smd3b_only, method = "spearman", label.x = 20, label.y = -45, color = '#df352e')+
  stat_cor(data=co_regulate, method = "spearman", label.x = 20, label.y = -50, color = '#FEB451')+
  ylim(-50, 50)+ xlim(-50, 50)+ 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x="smd3b/WT(∆PSI)",y="smd3a/WT(∆PSI)",title = "DTL-24")+
  scale_size(range = c(0, 4)) +
  scale_color_manual(limits = c( 'SmD3a', 'SmD3b',  'Co-regulated'),
                     values = c('#1072bd','red', '#FEB451')) 
ggsave("FigureS8E.pdf", width = 18, height = 15, units = "cm")
