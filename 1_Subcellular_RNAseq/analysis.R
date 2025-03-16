library(GenomicFeatures)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggtext)

## 1. calculate gene FPKM and TPM ####

# load read count matrix
c <- read.delim("subcellular.raw_counts.txt")
c.m <- c
c.m$Col.D.N <- c.m$Dark_N_1 + c.m$Dark_N_2 + c.m$Dark_N_3
c.m$Col.DL.N <- c.m$DTL6h_N_1 + c.m$DTL6h_N_2 + c.m$DTL6h_N_3
c.m$Col.WL.N <- c.m$WL_N_1 + c.m$WL_N_2 + c.m$WL_N_3

c.m$Col.D.F <- c.m$Dark_F_1 + c.m$Dark_F_2 + c.m$Dark_F_3
c.m$Col.DL.F <- c.m$DTL6h_F_1 + c.m$DTL6h_F_2 + c.m$DTL6h_F_3
c.m$Col.WL.F <- c.m$WL_F_1 + c.m$WL_F_2 + c.m$WL_F_3

c.m$Col.D.SM <- c.m$Dark_SM_1 + c.m$Dark_SM_2 + c.m$Dark_SM_3
c.m$Col.DL.SM <- c.m$DTL6h_SM_1 + c.m$DTL6h_SM_2 + c.m$DTL6h_SM_3
c.m$Col.WL.SM <- c.m$WL_SM_1 + c.m$WL_SM_2 + c.m$WL_SM_3

c.m$Col.D.P <- c.m$Dark_P_1 + c.m$Dark_P_2 + c.m$Dark_P_3
c.m$Col.DL.P <- c.m$DTL6h_P_1 + c.m$DTL6h_P_2 + c.m$DTL6h_P_3
c.m$Col.WL.P <- c.m$WL_P_1 + c.m$WL_P_2 + c.m$WL_P_3

c.m <- c.m[,c("Gene.ID","Col.D.N","Col.DL.N","Col.WL.N",
                  "Col.D.F","Col.DL.F","Col.WL.F",
                  "Col.D.P","Col.DL.P","Col.WL.P",
                  "Col.D.SM","Col.DL.SM","Col.WL.SM")]
dim(c.m)
write.csv(c.m, "merged_bam_gene_count.csv", row.names = F)

#读入统计counts的gtf
txdb <- makeTxDbFromGFF("merged_correct_exon.gtf",format="gtf")
columns(txdb)
#统计基因的exons
exons.list.per.gene <- exonsBy(txdb,by="gene")
#calculate effective length for each gene
exonic.gene.sizes <- as.data.frame(sum(width(reduce(exons.list.per.gene))))
exonic.gene.sizes$Gene.ID <- rownames(exonic.gene.sizes)
colnames(exonic.gene.sizes)[1] <- "len"
exonic.gene.sizes <- exonic.gene.sizes[,c(2,1)]

dim(exonic.gene.sizes)

c <- read.csv("merged_bam_gene_count.csv") 

#merge read count and effective length for each gene
m <- merge(c,exonic.gene.sizes,by="Gene.ID")

countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

countToEffCounts <- function(counts, len, effLen)
{
  counts * (len / effLen)
}

#计算各sample的FPKM
m$FPKM.Col.D.N <- with(m, countToFpkm(Col.D.N,len))
m$FPKM.Col.DL.N <- with(m, countToFpkm(Col.DL.N,len))
m$FPKM.Col.WL.N <- with(m, countToFpkm(Col.WL.N,len))

m$FPKM.Col.D.F <- with(m, countToFpkm(Col.D.F,len))
m$FPKM.Col.DL.F <- with(m, countToFpkm(Col.DL.F,len))
m$FPKM.Col.WL.F <- with(m, countToFpkm(Col.WL.F,len))

m$FPKM.Col.D.P <- with(m, countToFpkm(Col.D.P,len))
m$FPKM.Col.DL.P <- with(m, countToFpkm(Col.DL.P,len))
m$FPKM.Col.WL.P <- with(m, countToFpkm(Col.WL.P,len))

m$FPKM.Col.D.SM <- with(m, countToFpkm(Col.D.SM,len))
m$FPKM.Col.DL.SM <- with(m, countToFpkm(Col.DL.SM,len))
m$FPKM.Col.WL.SM <- with(m, countToFpkm(Col.WL.SM,len))

#计算各sample的TPM
m$TPM.Col.D.N <- with(m, countToTpm(Col.D.N,len))
m$TPM.Col.DL.N <- with(m, countToTpm(Col.DL.N,len))
m$TPM.Col.WL.N <- with(m, countToTpm(Col.WL.N,len))

m$TPM.Col.D.F <- with(m, countToTpm(Col.D.F,len))
m$TPM.Col.DL.F <- with(m, countToTpm(Col.DL.F,len))
m$TPM.Col.WL.F <- with(m, countToTpm(Col.WL.F,len))

m$TPM.Col.D.P <- with(m, countToTpm(Col.D.P,len))
m$TPM.Col.DL.P <- with(m, countToTpm(Col.DL.P,len))
m$TPM.Col.WL.P <- with(m, countToTpm(Col.WL.P,len))

m$TPM.Col.D.SM <- with(m, countToTpm(Col.D.SM,len))
m$TPM.Col.DL.SM <- with(m, countToTpm(Col.DL.SM,len))
m$TPM.Col.WL.SM <- with(m, countToTpm(Col.WL.SM,len))

m.fpkm <- m[,c("Gene.ID",colnames(m)[str_detect(colnames(m), pattern="FPKM.")])]
m.tpm <- m[,c("Gene.ID",colnames(m)[str_detect(colnames(m), pattern="TPM.")])]

write.csv(m.fpkm, "calculated_gene_FPKM.csv", row.names = F, quote = F)
write.csv(m.tpm, "calculated_gene_TPM.csv", row.names = F, quote = F)

## 2. TPM stack plot: Figure S1A ####

tpm <- read.csv("calculated_gene_TPM.csv")
colnames(tpm) <- str_remove(colnames(tpm), pattern = "TPM.")

# read gene type information
gt <- read.csv("merged_gene_annotation.csv")
gx <- merge(gt, tpm, by="Gene.ID")
colnames(gx)
table(gx$detailed_type)

gx.plot <- melt(gx, id=c(1:10), measure=11:ncol(gx), variable.name = "material")
head(gx.plot)

#reorder and rename levels
levels(gx.plot$detailed_type)
gx.plot$detailed_type <- factor(gx.plot$detailed_type,
                                levels = c("intergenic","intragenic",
                                           "snoRNA","snRNA","tRNA","mRNA"),
                                labels = c("inter-ncRNA","intra-ncRNA",
                                           "snoRNA","snRNA","tRNA","mRNA"))
## facet
gx.plot$subcellular[gx.plot$material %in% c("Col.D.N","Col.DL.N","Col.WL.N")] <- "Nuclei"  
gx.plot$subcellular[gx.plot$material %in% c("Col.D.F","Col.DL.F","Col.WL.F")] <- "Cytoplasm"  
gx.plot$subcellular[gx.plot$material %in% c("Col.D.SM","Col.DL.SM","Col.WL.SM")] <- "Monosome"  
gx.plot$subcellular[gx.plot$material %in% c("Col.D.P","Col.DL.P","Col.WL.P")] <- "Polysome"  

gx.plot$subcellular %>% as.factor() %>% levels()
gx.plot$subcellular <- factor(gx.plot$subcellular,
                              levels = c("Nuclei","Cytoplasm",
                                         "Monosome","Polysome"))

ggplot(gx.plot,aes(x=material, y=value, fill=detailed_type))+
  geom_bar(stat = "identity")+
  facet_grid(. ~ subcellular, scales = "free_x")+ #并排排列
  labs(x="", y="Transcripts per million (TPM)")+  #更改坐标轴标题
  ylim(0,1e6)+   #设置y坐标轴起始位置
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_text(size=rel(1)),  #更改坐标轴标题字体大小
        # axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1,size = rel(0.9)),
        axis.text.y = element_text(size = rel(0.9)),
        legend.text = element_text(size = rel(0.9)),
        legend.title = element_text(size = rel(0.9)),
        panel.spacing = unit(0.2, "lines"))+ #调整两个分面间的距离
  scale_fill_manual(values = c("#65499E","#BFABCF","#F58229","#37A849","#98C9DD","#207CB5"))+ #自定义颜色
  scale_x_discrete(labels=rep(c("D","DL","WL"),times=4))+ #修改横坐标轴标签
  labs(fill="Type") #修改图例标题
ggsave("FigureS1A.pdf",width = 16,height = 10, units = "cm")

## 3. scatter plot for snRNA, snoRNA and tRNA: Figure1B ####
fpkm <- read.csv("calculated_gene_FPKM.csv")
colnames(fpkm) <- str_remove(colnames(fpkm), pattern = "FPKM.")

# only keep genes with FPKM > 1 
fpkm1 <- fpkm[apply(fpkm[,2:ncol(fpkm)] > 1, 1, any),] 

fpkm1.gt <- merge(gt, fpkm1, by.y="Gene.ID")
table(fpkm1.gt$detailed_type)
fpkm1.sm <- subset(fpkm1.gt, detailed_type %in% c("snoRNA","snRNA","tRNA"))

#select Nucleus and Cytoplasm
cn.fpkm <- fpkm1.sm[,c("Gene.ID","detailed_type",colnames(fpkm1.sm)[str_detect(colnames(fpkm1.sm), pattern="[.N.F]$")])]
table(cn.fpkm$detailed_type)
cn.fpkm$Type <- factor(cn.fpkm$detailed_type, levels = c("snRNA","snoRNA","tRNA"))

cn.fpkm.D <- cn.fpkm[,c("Gene.ID","Type","Col.D.N","Col.D.F")]
cn.fpkm.D$condition <- "D"
colnames(cn.fpkm.D)[3:4] <- c("Nu","Cy")

cn.fpkm.DL <- cn.fpkm[,c("Gene.ID","Type","Col.DL.N","Col.DL.F")]
cn.fpkm.DL$condition <- "DtoL"
colnames(cn.fpkm.DL)[3:4] <- c("Nu","Cy")

cn.fpkm.WL <- cn.fpkm[,c("Gene.ID","Type","Col.WL.N","Col.WL.F")]
cn.fpkm.WL$condition <- "WL"
colnames(cn.fpkm.WL)[3:4] <- c("Nu","Cy")

cn.fpkm.plot <- rbind(cn.fpkm.D, cn.fpkm.DL, cn.fpkm.WL)
cn.fpkm1.plot <- cn.fpkm.plot[!str_detect(cn.fpkm.plot$Gene.ID, "ATCG|ATMG"),]

#构建注释文本的位置
ann_text <- data.frame(Nu=2^12, Cy=2^-2, Type="snRNA",
                       condition = factor("WL", levels = c("D","DtoL","WL")))
#title position
title.pos <- data.frame(Nu=rep(-6,3), Cy=rep(16,3), Type="snRNA",
                        condition = c("D","DtoL","WL"))

ggplot(cn.fpkm1.plot, aes(x=log2(Nu), y=log2(Cy),colour=Type)) + geom_point(size=0.3) +
  geom_abline(size = 0.25) + facet_grid(condition ~ .) +
  scale_x_continuous(limits = c(-8,18), breaks = c(seq(-5,15,5))) +
  scale_y_continuous(limits = c(-8,18), breaks = c(seq(-5,15,5))) +
  xlab(expression(paste('log'[2]*'Nu'["FPKM"]))) +
  ylab(expression(paste('log'[2]*'Cy'["FPKM"]))) +
  theme_test()+ scale_color_manual(values = c("#EE312D","#FFCB40","#0D74BB"))+
  geom_richtext(data=ann_text,
                label = "<span style = 'font-size:6pt; color:#EE312D;'>snRNA</span><br>
                <span style = 'font-size:6pt; color:#FFD92F;'>snoRNA</span><br>
                <span style = 'font-size:6pt; color:#0D74BB;'>tRNA</span>",
                fill = NA, label.colour = NA) +
  geom_text(data = title.pos, aes(x=Nu,y=Cy, label=condition), color="black", size=2) +
  theme(text = element_text(color = "black", size=7),
        line = element_line(size = 0.25),
        legend.position = "none",
        #panel.spacing = unit(0.1, "lines"), #调整两个分面间的距离
        strip.text = element_blank(),
        #strip.placement = "inside",
        strip.background = element_blank()
  )
ggsave("Figure1B.pdf", width = 4, height = 10, units = "cm")

## 4. bar plot for high expressed U5 on ribosomes: Figure S1D ####
u5.fpkm <- subset(fpkm1.gt, annotation=="U5")
u5.fpkm.ribo <- u5.fpkm[,c("Gene.ID", "Col.D.P", "Col.DL.P", "Col.WL.P", "Col.D.SM", "Col.DL.SM", "Col.WL.SM")]

u5.fpkm.ribo.plot <- melt(u5.fpkm.ribo, id.vars = "Gene.ID", variable.name = "condition",value.name = "FPKM")

u5.fpkm.ribo.plot$material[str_detect(u5.fpkm.ribo.plot$condition, pattern = "Col.D[.]")] <- "DK"
u5.fpkm.ribo.plot$material[str_detect(u5.fpkm.ribo.plot$condition, pattern = "Col.DL[.]")] <- "DTL-6"
u5.fpkm.ribo.plot$material[str_detect(u5.fpkm.ribo.plot$condition, pattern = "Col.WL[.]")] <- "WL"

u5.fpkm.ribo.plot$material <- factor(u5.fpkm.ribo.plot$material,
                                   levels = c("DK", "DTL-6", "WL"))

u5.fpkm.ribo.plot$subcellular[u5.fpkm.ribo.plot$condition %in% c("Col.D.SM","Col.DL.SM","Col.WL.SM")] <- "Monosome"  
u5.fpkm.ribo.plot$subcellular[u5.fpkm.ribo.plot$condition %in% c("Col.D.P","Col.DL.P","Col.WL.P")] <- "Polysome"  

u5.fpkm.ribo.plot$subcellular <- factor(u5.fpkm.ribo.plot$subcellular,
                                           levels = c("Polysome", "Monosome"))
u5.fpkm.ribo.plot$Gene.ID <- str_remove(u5.fpkm.ribo.plot$Gene.ID, "^at")
# keep top five expressed U5
u5.fpkm.ribo.plot.keep <- subset(u5.fpkm.ribo.plot, Gene.ID %in% c("U5.1b", "U5-3", "U5-5", "U5-4", "U5-6")) %>% droplevels()
ggplot(u5.fpkm.ribo.plot.keep, aes(x=reorder(Gene.ID, -FPKM), y=FPKM, fill=subcellular))+
  geom_col(position="stack") +
  xlab("")+ylab("FPKM") + facet_grid(. ~ material) +
  scale_y_continuous(expand = c(0,0),limits = c(0,5800)) +
  scale_fill_manual(values = c("#EE312D","#207CB5")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90),
        panel.grid = element_blank(), #remove grid lines
        legend.position = c(0.85,0.75),
        legend.title = element_blank(),
        legend.background = element_blank(), #remove legend background
        legend.key = element_blank()) 
ggsave("Figure1D.pdf", width =10, height = 8,units = "cm")

## 5. DESeq2: analyze transcripts with changes in subcellular localization ####
library(DESeq2)

# load read count matrix
cts <- read.delim("subcellular.raw_counts.txt", row.names = "Gene.ID")

# delete non-expressed genes
cts <- cts[apply(cts[,1:ncol(cts)]>0, 1, any),]  

#构建每个sample所对应的条件
subcellular <- rep(c("N","F","P","SM"), each=9)
material <- rep(c("WT_D","WT_DL","WT_WL"),each=3,times=4)

coldata <- data.frame(row.names = colnames(cts),
                      subcellular=subcellular, material=material)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ material)
ddsMF <- dds
ddsMF$group = factor(paste0(ddsMF$subcellular,".",ddsMF$material))
design(ddsMF) = ~0+group
ddsMF <- DESeq(ddsMF)

#compare nuclei and cytoplasm in WT-D
res.WTD_N_F <- results(ddsMF, contrast = c("group","F.WT_D","N.WT_D")) #log2FC为F-N
#compare nuclei and cytoplasm in WT-DL
res.WTDL_N_F <- results(ddsMF, contrast = c("group","F.WT_DL","N.WT_DL")) #log2FC为F-N
#compare nuclei and cytoplasm in WT-WL
res.WTWL_N_F <- results(ddsMF, contrast = c("group","F.WT_WL","N.WT_WL")) #log2FC为F-N

re.N_F <- cbind(res.WTD_N_F[,c("log2FoldChange","pvalue","padj")],
            res.WTDL_N_F[,c("log2FoldChange","pvalue","padj")],
            res.WTWL_N_F[,c("log2FoldChange","pvalue","padj")])
colnames(re.N_F) <- c("WTD_N_F.FC","WTD_N_F.p","WTD_N_F.q",
                  "WTDL_N_F.FC","WTDL_N_F.p","WTDL_N_F.q",
                  "WTWL_N_F.FC","WTWL_N_F.p","WTWL_N_F.q")
write.csv(re.N_F, file = "N_F.DESeq2.csv", quote = F)


#compare monosome and cytoplasm in WT-D
res.WTD_SM_P <- results(ddsMF, contrast = c("group","P.WT_D","SM.WT_D")) #log2FC为F-N
#compare monosome and cytoplasm in WT-DL
res.WTDL_SM_P <- results(ddsMF, contrast = c("group","P.WT_DL","SM.WT_DL")) #log2FC为F-N
#compare monosome and cytoplasm in WT-WL
res.WTWL_SM_P <- results(ddsMF, contrast = c("group","P.WT_WL","SM.WT_WL")) #log2FC为F-N

re.SM_P <- cbind(res.WTD_SM_P[,c("log2FoldChange","pvalue","padj")],
                res.WTDL_SM_P[,c("log2FoldChange","pvalue","padj")],
                res.WTWL_SM_P[,c("log2FoldChange","pvalue","padj")])
colnames(re.SM_P) <- c("WTD_SM_P.FC","WTD_SM_P.p","WTD_SM_P.q",
                      "WTDL_SM_P.FC","WTDL_SM_P.p","WTDL_SM_P.q",
                      "WTWL_SM_P.FC","WTWL_SM_P.p","WTWL_SM_P.q")
write.csv(re.SM_P, file = "SM_P.DESeq2.csv", quote = F)

# get log2(normalized counts)
x <- mcols(ddsMF,use.names = T) %>% as.data.frame()
nc <- x[,str_detect(colnames(x),"^group")]
colnames(nc) <- str_replace(string = colnames(nc), pattern = "group", replacement = "")
nc <- nc[,c("N.WT_D","N.WT_DL","N.WT_WL",
            "F.WT_D","F.WT_DL","F.WT_WL",
            "P.WT_D","P.WT_DL","P.WT_WL",
            "SM.WT_D","SM.WT_DL","SM.WT_WL")]
write.csv(nc, "all_log2normalized_count.csv")

## 6. Heatmap for snRNA: Figure1C ####
library(dplyr)
library(ComplexHeatmap)
library(circlize)

fpkm <- read.csv("calculated_gene_FPKM.csv")
colnames(fpkm) <- str_remove(colnames(fpkm), pattern = "FPKM.")

fpkm.sn <- fpkm[str_detect(fpkm$Gene.ID, "^atU"), ] 

# select snRNA with FPKM median > 1
fpkm.sn1 <- fpkm.sn %>% 
  mutate(median_per_row = apply(fpkm.sn[,2:ncol(fpkm.sn)], 1, median)) %>% 
  filter(median_per_row > 1) 

# load log2 normalized counts
norm.c <- read.csv("all_log2normalized_count.csv") #39882
colnames(norm.c) <- colnames(fpkm)
normc.sn <- norm.c[norm.c$Gene.ID %in% fpkm.sn1$Gene.ID,] #caculate FPKM:23816

de.N_F <- read.csv("N_F.DESeq2.csv")
de.N_F.normc.sn <- merge(normc.sn, de.N_F, by="Gene.ID")
dim(de.N_F.normc.sn)

# read gene type information
gt <- read.csv("merged_gene_annotation.csv")
de.N_F.normc.gt.sn <- merge(gt, de.N_F.normc.sn, by="Gene.ID")

# order as U1, U2, U4, U5, U6
clus.sn <- arrange(de.N_F.normc.gt.sn, annotation, WTD_N_F.FC)
rownames(clus.sn) <- clus.sn$symbol
rownames(clus.sn) <- str_remove(rownames(clus.sn), "^at")
keep.order <- rownames(clus.sn) 
write.table(keep.order, "Cy_Nu_FC_heatmap.order", row.names = F,
            col.names = F, quote = F)
clus.sn.N_F <- clus.sn[,c("WTD_N_F.FC", "WTDL_N_F.FC", "WTWL_N_F.FC")]
colnames(clus.sn.N_F) <- c("DK","DTL-6","WL")

fc.sn.N_F.rot <- t(clus.sn.N_F)
p1 <- Heatmap(fc.sn.N_F.rot,
  gap = unit(2,"mm"),
  cluster_rows = F,
  cluster_columns = F,
  col = colorRamp2(c(-3,0,3), c("#0D74BB","white", "#EE312D")), 
  row_title_rot = 0, 
  row_title_gp = gpar(fontsize=11),
  rect_gp = gpar(col = "white",
                 lty = 1,
                 lwd = 0.5),
  na_col = "darkgrey",name = "Cy/Nu",
  column_names_rot = 90, column_names_side = "top",
  border = F,
  show_row_names = T,row_names_side = "left", row_names_gp = gpar(fontsize=10),
  )
gb.p1 = grid.grabExpr(draw(p1))
cowplot::plot_grid(gb.p1)

# Monosome and Polysome
de.ribo <- read.csv("SM_P.DESeq2.csv")
colnames(de.ribo)[1] <- "Gene.ID"
de.ribo.sn <- de.ribo[str_detect(de.ribo$Gene.ID, "^atU"), ] 
de.ribo.sn$Gene.ID <- str_remove(de.ribo.sn$Gene.ID, "^at")
#order as Nu/Cy heatmap
order <- read.table("Cy_Nu_FC_heatmap.order")
colnames(order) <- "Gene.ID"
sn.ribo <- left_join(order, de.ribo.sn[,c("Gene.ID","WTD_SM_P.FC", "WTDL_SM_P.FC", "WTWL_SM_P.FC")], 
                     by="Gene.ID")
rownames(sn.ribo) <- sn.ribo$Gene.ID
sn.ribo.m <- sn.ribo[,c(2:4)]
colnames(sn.ribo.m) <- c("DK","DTL-6","WL")

#rotation
sn.ribo.m.rot <- t(sn.ribo.m)
p2 <- sn.ribo.m.rot  %>%  Heatmap(
  gap = unit(2,"mm"),
  cluster_rows = F,
  cluster_columns = F,
  col = colorRamp2(c(-3,0,4), c("#EE312D","lightgrey", "#0D74BB")),
  #split = clus.sn$clus.1,
  row_title_rot = 0, #调整clus.1方向
  row_title_gp = gpar(fontsize=11),
  #row_order = order,
  rect_gp = gpar(col = "white",
                 lty = 1,
                 lwd = 0.5),
  na_col = "darkgrey",name = "P/M",
  column_names_rot = 90, column_names_side = "top",
  border = F,
  show_column_names = F,
  show_row_names = T,row_names_side = "left", row_names_gp = gpar(fontsize=10),
  )
gb.p2 = grid.grabExpr(draw(p2))

p.N_F <- ggplotify::as.ggplot(gb.p1)
p.SM_P <- ggplotify::as.ggplot(gb.p2)
p.N_F + p.SM_P

pdf("Figure1C.pdf", width = 10, height = 4)
cowplot::plot_grid(gb.p1, gb.p2, align = "hv", ncol = 1,
                   rel_widths = c(1, 1), rel_heights = c(1.5,1))
dev.off()


