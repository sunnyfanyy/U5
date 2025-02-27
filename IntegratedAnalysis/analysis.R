library(openxlsx)
library(stringr)
library(dplyr)
library(ComplexHeatmap)

## venn diagram: DEGs and changed TE&Proteome in u5 mutants, Figure3D ####
library(VennDiagram)

# DTL24h RNAseq
rna <- read.xlsx("../DTL24h_mRNAseq/DTL24h_WT_u5.merged_info.xlsx")
u5.DEGs <- subset(rna, abs(DTL24h.u5_WT.FC) > log2(1.5) & DTL24h.u5_WT.p < 0.05)
u5.DEGs.up <- subset(u5.DEGs, DTL24h.u5_WT.FC > 0) # 460
u5.DEGs.down <- subset(u5.DEGs, DTL24h.u5_WT.FC < 0) # 480

# Translation efficiency
FPKM.s.TE <- read.xlsx("../Polysome_seq/PolysomeSeq_TE_info.xlsx")
# significantly changed TE
u5.TE <- subset(FPKM.s.TE, pvalue < 0.05)
u5.TE.down <- subset(u5.TE, log2FC_TE.u5_WT < 0)
u5.TE.up <- subset(u5.TE, log2FC_TE.u5_WT > 0)

# TMT-based proteome results
prot <- read.xlsx("../Proteome/DTL24h_TMT_proteome.xlsx")
colnames(prot)

prot.u5.down <- subset(prot, `t.test_24h.p` < 0.05 & `DTL24h.u5/WT` < 0.88)
prot.u5.up <- subset(prot, `t.test_24h.p` < 0.05 & `DTL24h.u5/WT` > 1.2)

# TE and TMT：both up-regulated in u5
u5.TE_TMT.up <- u5.TE.up[u5.TE.up$Gene.ID %in% prot.u5.up$Gene.ID, ] 
# TE and TMT：both down-regulated in u5
u5.TE_TMT.down <- u5.TE.down[u5.TE.down$Gene.ID %in% prot.u5.down$Gene.ID, ]
dim(u5.TE_TMT.up) # 41
dim(u5.TE_TMT.down) #242


u5.TE.down$`TE_u5/WT` <- u5.TE.down$Avg.TE.u5/u5.TE.down$Avg.TE.WT

TE_TMT.dn.both <- merge(u5.TE.down[,c("Gene.ID","Avg.TE.WT","Avg.TE.u5","TE_u5/WT","pvalue")], 
                        prot.u5.down[,c("Gene.ID", "Avg.24h.WT","Avg.24h.u5", "DTL24h.u5/WT", "t.test_24h.p",
                                        "symbol", "full_name", "Note", "curator_summary", "Alias")], by="Gene.ID")
colnames(TE_TMT.dn.both)[5:9] <- c("TE.pvalue", "Avg.TMT.WT", "Avg.TMT.u5", "TMT_u5/WT", "TMT.pvalue") 
# [save results] TE and TMT：both down-regulated in u5
#write.xlsx(TE_TMT.dn.both, "TE_TMT_both_down_in_u5.xlsx")
#write.table(TE_TMT.dn.both$Gene.ID, "TE_TMT_down_g242.genelist.txt",
#            quote = F, row.names = F, col.names = F)


library(VennDiagram)
library(cowplot)
library(ggplotify)
venn.plot.up <- VennDiagram::venn.diagram(list(ΔmRNA.up=u5.DEGs.up$Gene.ID,
                                                 ΔmRNA.dn=u5.DEGs.down$Gene.ID,
                                                 ΔTE.up=u5.TE_TMT.up$Gene.ID), 
                                            col=c("#009E73","#E69F00", "#56B4E9"),
                                            cat.cex=0.8, # 标签的大小
                                            cat.dist=c(0.1,0.1,0.05), #调节标签位置,离circle的距离
                                            ext.length=0.8, # 延长线的长度
                                            lwd=1.5, # circle边缘线的粗细
                                            #cat.pos=c(-20,20,0), #调节标签位置,在圆形的几点钟方向
                                            width=800,height = 700, resolution = 300,
                                            #main.cex = 1.1,main.just = c(0.5,1), #调整标题大小和位置
                                            main.fontfamily = "sans", cat.fontfamily = "sans", # Helvetica字体
                                            fontfamily = "sans",
                                            margin=c(0.11),
                                            filename=NULL)
as.ggplot(plot_grid(grobTree(venn.plot.up)))

venn.plot.down <- VennDiagram::venn.diagram(list(ΔmRNA.up=u5.DEGs.up$Gene.ID,
                                            ΔmRNA.dn=u5.DEGs.down$Gene.ID,
                                            ΔTE.dn=u5.TE_TMT.down$Gene.ID), 
                                       col=c("#009E73","#E69F00", "#56B4E9"),
                                       cat.cex=0.8, # 标签的大小
                                       cat.dist=c(0.1,0.1,0.05), #调节标签位置,离circle的距离
                                       ext.length=0.8, # 延长线的长度
                                       lwd=1.5, # circle边缘线的粗细
                                       #cat.pos=c(-20,20,0), #调节标签位置,在圆形的几点钟方向
                                       width=800,height = 700, resolution = 300,
                                       #main.cex = 1.1,main.just = c(0.5,1), #调整标题大小和位置
                                       main.fontfamily = "sans", cat.fontfamily = "sans", # Helvetica字体
                                       fontfamily = "sans",
                                       margin=c(0.11),
                                       filename=NULL)
as.ggplot(plot_grid(grobTree(venn.plot.down)))

## barplot for genes essential for chloroplast biogenesis and photosynthesis, Figure3H ####
library(dplyr)
library(tidyr)
library(ggplot2)

gl <- c("RPS1","GHS1", "PnsB2", "PnsB4", "PnsL4", "LHCB3", "LHCA4","FLU","FAD7","SKL1","UBQ5",
        "CHLM", "RBCS1B", "ALB3")

# RNAseq FPKM log2FC
fpkm.DTL24 <- read.xlsx("../DTL24h_mRNAseq/DTL24h_WT_u5.merged_info.xlsx")
fpkm.DTL24$fpkm.ratio <- (fpkm.DTL24$FPKM.DTL24h.u5+1e-5)/(fpkm.DTL24$FPKM.DTL24h.WT+1e-5)
fpkm.DTL24$fpkm.log2ratio <- log(fpkm.DTL24$fpkm.ratio, 2)

# TMT proteome log2FC
prot <- read.xlsx("../Proteome/DTL24h_TMT_proteome.xlsx")
prot$TMT.log2ratio <- log(prot$`DTL24h.u5/WT`, 2)

fpkm.prot.all <- inner_join(fpkm.DTL24, prot)


gl.all <- fpkm.prot.all[fpkm.prot.all$symbol %in% gl, ]

gl.keep <- gl.all[,c("symbol", "fpkm.log2ratio", "TMT.log2ratio")]

gl.keep.plot <- pivot_longer(gl.keep, cols = 2:3, names_to = "log2ratio")

gl.keep.plot$log2ratio <- factor(gl.keep.plot$log2ratio, 
                                 labels = c("mRNA", "protein") )
gl.keep.plot <- subset(gl.keep.plot, ! symbol %in% c("FLU", "FAD7", "SKL1", "CHLM"))
gl.keep.plot$symbol <- factor(gl.keep.plot$symbol, 
                              levels = c("UBQ5", "RPS1","GHS1", 
                                         "PnsB4","PnsL4",  "ALB3", "PnsB2",
                                         "LHCB3",  "RBCS1B", "LHCA4" ) )

ggplot(gl.keep.plot, aes(x=symbol, y=value, fill=log2ratio)) +
  geom_col() + 
  coord_flip() + #旋转坐标轴
  scale_y_reverse(position = "right", limits=c(1.3,-1.3)) +
  labs(x="", y=expression(paste(log[2]," fold change"))) +
  scale_fill_manual(values = c("#DF352F","#196EAF") ) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        legend.title = element_blank(),
        legend.position = c(0.82,0.1),
        legend.background = element_blank(), #移除整体的边框
        legend.key = element_blank(), #移除每个图例项目周围的边框
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
  )
ggsave("Figure3H.pdf", width = 10, height = 10, units = "cm")

## Heatmap: time-series expression profiles of transcription and translation, Figure3I&J ####
library(openxlsx)
library(stringr)
library(ComplexHeatmap)
library(circlize)

# TMT and TE both downregulted in u5
g242 <- read.table("TE_TMT_down_g242.genelist.txt")
colnames(g242) <- "Gene.ID"
norm.Avg1.anno <- read.xlsx("../Proteome/eLife_all_normalized_proteome.xlsx")

norm.Avg1.anno$Avg_24h_12h.rela <- norm.Avg1.anno$Avg_24h.rela/norm.Avg1.anno$Avg_12h.rela
prot.24h_12h.up <- filter(norm.Avg1.anno, Avg_24h_12h.rela > 1.25) 
dim(prot.24h_12h.up) #1189

g242.norm <- norm.Avg1.anno[norm.Avg1.anno$Gene.ID %in% g242$Gene.ID, ] #188

#keep records with peptides >=2 
g242.norm <- subset(g242.norm, Peptide.count >=2) #188

g242.rela <- g242.norm[,c("Gene.ID",str_subset(colnames(g242.norm), "rela"))]

colnames(g242.rela)[-1] <- c("0","4","8","12","24","48","72","96")
g242.rela <- na.omit(g242.rela)
# only keep 0h, 8h, 12h, 24h, 48h
g242.rela <- g242.rela[,c("Gene.ID", "0", "8", "12", "24", "48")]
rownames(g242.rela) <- g242.rela$Gene.ID

g242.rela$group[g242.rela$Gene.ID %in% prot.24h_12h.up$Gene.ID] <- "Group1"
g242.rela$group[is.na(g242.rela$group)] <- "Group2"
table(g242.rela$group)

# FPKM and DESeq2 results
m <- read.xlsx("D_WL.DTL12h_24h_48h.merged_info.xlsx")
WT.12h_24h.nosig <- filter(m, abs(WT.DTL24h_DTL12h.FC) < 1 ) # 27844

group1.prot <- subset(g242.rela, group=="Group1") 
g108 <- group1.prot[group1.prot$Gene.ID %in% WT.12h_24h.nosig$Gene.ID, ] 
g108.prot <- group1.prot[group1.prot$Gene.ID %in% g108$Gene.ID, ]

set.seed(20231126)
p.g108.rela <- g108.prot[,2:6]  %>%  Heatmap(
  cluster_rows = T,
  #col = colorRamp2(c(0, 3), c("lightblue1", "blue")),
  col = colorRamp2(c(0,1,2,2.5,3,4), 
                   c("#EEEEEE","#E7EDDC","#F6E4AB","#FBA35C","#F06640","#D42D26")),
  #  col = colorRamp2(c(0,1,2,4,6,8), 
  #                   c("#EEEEEE","#E7EDDC","#F6E4AB","#FBA35C","#F06640","#D42D26")),
  cluster_columns = F,
  #split = g242.rela$group,
  show_row_dend = F,
  #row_order = row_order(p.g108.fpkm.rela), 
  #rect_gp = gpar(col = "white",
  #               lty = 1,
  #               lwd = 0.5),
  na_col = "darkgrey", name = "Ratio",
  #heatmap_legend_param = list(title = expression(log[2]*(FPKM))), # []下标，*连接符
  #column_labels = rep(c("D","12h","24h","48h"),3), #将列名替换为指定字符
  column_names_rot = 0, column_names_gp = gpar(fontsize=11),
  column_names_side = "bottom",
  border = F,
  show_row_names = F,row_names_gp = gpar(fontsize=4),
  #width = unit(7, "cm"),
  #height = unit(9, "cm")
)
p.g108.rela <- draw(p.g108.rela)

# RNAseq: time-series expression fold change
# rRNA-depleted RNAseq
fpkm.6h <- read.csv("DTL6h_calculated_gene_FPKM.csv")
colnames(fpkm.6h)
colnames(fpkm.6h) <- str_remove(colnames(fpkm.6h), "FPKM.WT_")
fpkm.6h$DTL6h.rela <- (fpkm.6h$DTL6h+1e-5)/(fpkm.6h$D+1e-5)

# mRNAseq
fpkm.DTL12_24_48 <- read.xlsx("D_WL.DTL12h_24h_48h.merged_info.xlsx")
fpkm.DTL12_24_48 <- fpkm.DTL12_24_48[,c("Gene.ID", "FPKM.D.WT", "FPKM.DTL12h.WT",
                                        "FPKM.DTL24h.WT", "FPKM.DTL48h.WT")]
colnames(fpkm.DTL12_24_48) <- str_remove(colnames(fpkm.DTL12_24_48), "FPKM.")
colnames(fpkm.DTL12_24_48) <- str_remove(colnames(fpkm.DTL12_24_48), ".WT")

fpkm.DTL12_24_48$D.rela <- 1
fpkm.DTL12_24_48$DTL12h.rela <- fpkm.DTL12_24_48$DTL12h/(fpkm.DTL12_24_48$D+1e-5)
fpkm.DTL12_24_48$DTL24h.rela <- fpkm.DTL12_24_48$DTL24h/(fpkm.DTL12_24_48$D+1e-5)
fpkm.DTL12_24_48$DTL48h.rela <- fpkm.DTL12_24_48$DTL48h/(fpkm.DTL12_24_48$D+1e-5)
fpkm.DTL12_24_48.rela <- fpkm.DTL12_24_48[,c("Gene.ID",str_subset(colnames(fpkm.DTL12_24_48), "rela"))]

fpkm.all.rela <- inner_join(fpkm.6h[,c("Gene.ID","DTL6h.rela")], fpkm.DTL12_24_48.rela, by="Gene.ID")
colnames(fpkm.all.rela) <- str_remove(colnames(fpkm.all.rela), ".rela")
fpkm.all.rela <- fpkm.all.rela[,c(1, 3, 2, 4:6)]
g108.fpkm <- fpkm.all.rela[fpkm.all.rela$Gene.ID %in% g108$Gene.ID, ]

colnames(g108.fpkm)[2:6] <- c("0", "6", "12", "24", "48")
set.seed(20231126)
p.g108.fpkm.rela <- g108.fpkm[,2:6]  %>%  Heatmap(
  cluster_rows = F, #col = colorRamp2(c(0, 3), c("lightblue1", "blue")),
  col = colorRamp2(c(0,1,2,2.5,3,4), 
                   c("#EEEEEE","#E7EDDC","#F6E4AB","#FBA35C","#F06640","#D42D26")),
  cluster_columns = F,
  row_order = row_order(p.g108.rela),
  #split = fpkm.g242.rela$group,
  na_col = "darkgrey", name = "Ratio",
  #heatmap_legend_param = list(title = expression(log[2]*(FPKM))), # []下标，*连接符
  #column_labels = rep(c("D","12h","24h","48h"),3), #将列名替换为指定字符
  column_names_rot = 0, column_names_gp = gpar(fontsize=11),
  column_names_side = "bottom", column_names_centered = T,
  border = F,
  show_row_names = F,row_names_gp = gpar(fontsize=4),
  #width = unit(7, "cm"),
  #height = unit(9, "cm")
)
p.g108.fpkm.rela <- draw(p.g108.fpkm.rela)

gb.p.g108.fpkm.rela = grid.grabExpr(draw(p.g108.fpkm.rela))
gb.p.g108.rela = grid.grabExpr(draw(p.g108.rela))
cowplot::plot_grid(gb.p.g108.fpkm.rela, gb.p.g108.rela, ncol=2)
ggsave("Figure3I&J.pdf", width=15, height = 10, units = "cm")


## Heatmap for U5 targets: DTL24h RNAseq and TMT_proteome, Figure S3A&E ####
# U5 targets: RPS1, GHS1, LHCB3, PnsB2, PnsB4, PnsL4
gl.c <- c("RPS1","GHS1","LHCB3","PnsB2","PnsB4","PnsL4")


u5.t <- prot[prot$symbol %in% gl.c,]

u5.t1 <- u5.t[,c("Gene.ID","symbol","DTL24h.u5/WT","t.test_24h.p",
                 "Norm.DTL24h.WT_1","Norm.DTL24h.WT_2","Norm.DTL24h.WT_3",
                 "Norm.DTL24h.u5_1","Norm.DTL24h.u5_2","Norm.DTL24h.u5_3")]
colnames(u5.t1) <- str_replace(colnames(u5.t1),"Norm.","")
colnames(u5.t1)[5:10] <- c("WT_R1","WT_R2","WT_R3","*u5-3 u5-4*_R1","*u5-3 u5-4*_R2","*u5-3 u5-4*_R3")

# Create a matrix
u5.t1.m <- u5.t1[,5:10]  %>% 
  as.matrix()
# assign rownames
rownames(u5.t1.m) <- u5.t1$symbol

# z-score
u5.t1.m <- u5.t1.m %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()


ht1 = Heatmap(u5.t1.m, name="z-score",show_row_names = T, 
              cluster_rows = T, cluster_columns = F,
              column_labels = gt_render(colnames(u5.t1)[5:10]), 
              column_title = "Proteome") 

# DTL24h RNAseq
# sample normalized read counts
nc <- read.csv("../DTL24h_mRNAseq/DTL24h.VST_DESeq2_normalized_counts_sample.csv")
u5.t2 <- nc[nc$Gene.ID %in% u5.t1$Gene.ID,]
u5.t2 <- inner_join(u5.t2, u5.t1[,c("Gene.ID","symbol")],by="Gene.ID")
colnames(u5.t2)[2:7] <- c("WT_R1","WT_R2","WT_R3","*u5-3 u5-4*_R1","*u5-3 u5-4*_R2","*u5-3 u5-4*_R3")

# Create a matrix
u5.t2.m <- u5.t2[,2:7]  %>% 
  as.matrix()
# assign rownames
rownames(u5.t2.m) <- u5.t2$symbol

# z-score
u5.t2.m <- u5.t2.m %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

ht2 =   Heatmap(u5.t2.m, name="z-score",show_row_names = T, 
                cluster_rows = T, cluster_columns = F,
                column_labels = gt_render(colnames(u5.t2)[2:7]),
                column_title = "RNAseq")
ht2 + ht1


#pdf("FigureS3A&E.pdf",width = 15/2.54, height = 12/2.54)
draw(ht2 + ht1, auto_adjust = FALSE)
#dev.off()

## Heatmap for U5 targets: time-series RNAseq and TMT_proteome, Figure S4A&B ####
gl <- data.frame(Gene=c("AT1G64770_PnsB2","AT5G54270_LHCB3", 
                        "AT5G30510_RPS1","AT2G29550_TUB7"))
gl$Gene.ID <- str_split_fixed(gl$Gene, "_", 2)[,1]
gl$symbol <- str_split_fixed(gl$Gene, "_", 2)[,2]

fpkm.u5.gl.rela <- right_join(fpkm.all.rela, gl[,c("Gene.ID", "symbol")], by="Gene.ID")
rownames(fpkm.u5.gl.rela) <- fpkm.u5.gl.rela$symbol

colnames(fpkm.u5.gl.rela)[2:6] <- c("0", "6", "12", "24", "48")

set.seed(123)
p <-  fpkm.u5.gl.rela[2:6] %>%  Heatmap(
  cluster_rows = T, column_title = "RNA-seq",
  col = colorRamp2(c(0,1,7), c("blue", "#EEEEEE", "red")),
  cluster_columns = F,
  rect_gp = gpar(col = "white",
                 lty = 1,
                 lwd = 0.5),
  na_col = "darkgrey", name = "Ratio",
  column_names_rot = 0, column_names_gp = gpar(fontsize=10),
  column_names_side = "bottom", column_names_centered = T,
  border = T,
  show_row_names = T,row_names_gp = gpar(fontsize=10),
  #width = unit(7, "cm"),
  #height = unit(9, "cm")
)
pdf("FigureS4A.pdf", height = 5/2.54, width = 12/2.54)
p
dev.off()


U5.t.norm <- norm.Avg1.anno[norm.Avg1.anno$Gene.ID %in% 
                         gl$Gene.ID,]
U5.t.norm <- left_join(U5.t.norm, gl[,c("Gene.ID","symbol")], by="Gene.ID")

U5.t.plot <- U5.t.norm[,c(str_subset(colnames(U5.t.norm),"rela"))]
rownames(U5.t.plot) <- U5.t.norm$symbol
colnames(U5.t.plot) <- c("0","4","8","12","24","48","72","96")

set.seed(123)
p <- U5.t.plot[,c(1,3:6)]  %>%  Heatmap(
  cluster_rows = T, column_title = "Proteome",
  col = colorRamp2(c(0,1,4), c("blue", "#EEEEEE", "red")),
  cluster_columns = F,
  rect_gp = gpar(col = "white",
                 lty = 1,
                 lwd = 0.5),
  na_col = "darkgrey", name = "Ratio",
  column_names_rot = 0, column_names_gp = gpar(fontsize=11),
  column_names_side = "bottom",
  border = T,
  show_row_names = T,row_names_gp = gpar(fontsize=11)
)
pdf("FigureS4B.pdf", height = 5/2.54, width = 12/2.54)
p
dev.off()
