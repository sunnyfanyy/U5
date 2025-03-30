library(openxlsx)
library(stringr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

##  Histogram: RNAseq and Proteome in u5 mutants, Figure3A ####
# DTL24h RNAseq
rna <- read.xlsx("../2_DTL24h_mRNAseq/DTL24h_WT_u5.merged_info.xlsx")

# TMT-based proteome results
prot <- read.xlsx("../4_Proteome/DTL24h_TMT_proteome.xlsx")
rna_prot <- inner_join(prot[,c("Gene.ID","Avg.24h.WT","Avg.24h.u5","DTL24h.u5/WT","t.test_24h.p")], 
                       rna[,c(1:5)], by="Gene.ID") # 8636
rna_prot.uniq <- rna_prot[!duplicated(rna_prot$Gene.ID), ] # 8620

rna_prot.keep <- subset(rna_prot.uniq, (FPKM.DTL24h.WT > 0 | FPKM.DTL24h.u5 > 0) & 
                          `t.test_24h.p` < 0.05 & !is.na(rna_prot.uniq$`DTL24h.u5/WT`)) # 3129

# RNA average and variance
mean(rna_prot.keep$DTL24h.u5_WT.FC) # 0.0161
sd(rna_prot.keep$DTL24h.u5_WT.FC) # 0.313
# Proteome average and variance
mean(log2(rna_prot.keep$`DTL24h.u5/WT`), na.rm = T) # 0.0941
sd(log2(rna_prot.keep$`DTL24h.u5/WT`)) # 0.636

rna_prot.keep$log2ProFC <- log2(rna_prot.keep$`DTL24h.u5/WT`)

ggplot()+geom_histogram(data = rna_prot.keep,aes(x=DTL24h.u5_WT.FC),stat = 'bin',bins = 100,fill='#4f81bd',color='white',alpha=0.7)+
  geom_smooth(data = rna_prot.keep,aes(x=DTL24h.u5_WT.FC),stat = 'bin',color='#4f81bd',bins = 100,size=0.6)+
  geom_histogram(data = rna_prot.keep,aes(x=log2ProFC),stat = 'bin',bins = 100,fill='#c0504d',color='white',alpha=0.5)+
  geom_smooth(data = rna_prot.keep,aes(x=log2ProFC),stat = 'bin',color='#c0504d',bins = 100,size=0.6)+
  scale_x_continuous(expand = c(0, 0),limits = c(-1.5,1.5)) +scale_y_continuous(expand = c(0, 0),limits = c(0,250))+
  xlab(expression(log[2]*"FC"))+ylab("Number of genes") +
  theme_classic() + theme(axis.text = element_text(color = 'black'))

#ggsave("Figure3A.pdf",width = 8, height = 8, units = "cm")

##  Histogram: RNAseq and TE in u5 mutants, Figure3B ####
# Translation efficiency
FPKM.s.TE <- read.xlsx("../3_Polysome_seq/PolysomeSeq_TE_info.xlsx")
colnames(FPKM.s.TE)

rna_TE <- inner_join(FPKM.s.TE[,1:23], 
                     rna[,c(1:5)], by="Gene.ID") # 29538
subset(rna_TE, Gene.ID %in% c("atU1-4", "atU1-6"))

#保留TE pvalue < 0.05 并且 FPKM任一sample > 0 的基因
rna_TE.sig <- subset(rna_TE, pvalue < 0.05 & 
                       (FPKM.DTL24h.WT > 0 | FPKM.DTL24h.u5 > 0)) # 2725
#去除Polysome-seq中3个及以上sample（共12个sample）的FPKM值为0的基因。
rna_TE.keep <- rna_TE.sig[which(rowSums(rna_TE.sig[,2:13]==0) < 3),]  # 2673

ggplot()+geom_histogram(data = rna_TE.keep,aes(x=DTL24h.u5_WT.FC),stat = 'bin',bins = 100,fill='#4f81bd',color='white',alpha=0.7)+
  geom_smooth(data = rna_TE.keep,aes(x=DTL24h.u5_WT.FC),stat = 'bin',color='#4f81bd',bins = 100,size=0.6)+
  geom_histogram(data = rna_TE.keep,aes(x=log2FC_TE.u5_WT),stat = 'bin',bins = 100,fill='#c0504d',color='white',alpha=0.5)+
  geom_smooth(data = rna_TE.keep,aes(x=log2FC_TE.u5_WT),stat = 'bin',color='#c0504d',bins = 100,size=0.6)+
  scale_x_continuous(expand = c(0, 0),limits = c(-2.5,2.5)) +scale_y_continuous(expand = c(0, 0),limits = c(0,350))+
  xlab(expression(log[2]*"FC"))+ylab("Number of genes") +
  theme_classic() + theme(axis.text = element_text(color = 'black'))

#ggsave("Figure3B.pdf",width = 8, height = 8, units = "cm")

#TE: up-regulated and down-regulated
rna_TE.keep.up <- subset(rna_TE.keep, log2FC_TE.u5_WT > 0) # 747
rna_TE.keep.dn <- subset(rna_TE.keep, log2FC_TE.u5_WT < 0) #1926


# RNA average and variance
mean(rna_TE.keep$DTL24h.u5_WT.FC) # 0.026
sd(rna_TE.keep$DTL24h.u5_WT.FC) # 0.318
# TE average and variance
mean(rna_TE.keep$log2FC_TE.u5_WT, na.rm = T) # -0.2229
sd(rna_TE.keep$log2FC_TE.u5_WT, na.rm = T) # 0.6625

## dotplot: upregulated and downregulated TE genes in u5 mutants, Figure3C ####
# 将B图中的2673个基因分为RS上下调和TE上下调的组合。
# RNAseq criteria: abs(log2FC) > log2(1.5) & p < 0.05
# TE criteria：abs(log2FC) >= 0.5 & p < 0.05
rna_TE.keep1 <- rna_TE.keep[,c(1, 22:23, 26:27)] # 2673

rna_TE.keep1.TE.sig <- subset(rna_TE.keep1, abs(log2FC_TE.u5_WT) >= 0.5 & pvalue < 0.05) # 956

rna_TE.keep1.TE.up <- subset(rna_TE.keep1.TE.sig, log2FC_TE.u5_WT >= 0.5 ) # 284
rna_TE.keep1.TE.dn <- subset(rna_TE.keep1.TE.sig, log2FC_TE.u5_WT <= -0.5 ) #672

rna_TE.keep1.rna.up <- subset(rna_TE.keep1.TE.sig, DTL24h.u5_WT.FC > log2(1.5) & DTL24h.u5_WT.p < 0.05) # 18
rna_TE.keep1.rna.dn <- subset(rna_TE.keep1.TE.sig, DTL24h.u5_WT.FC < -log2(1.5) & DTL24h.u5_WT.p < 0.05) # 24
rna_TE.keep1.rna.ns <- subset(rna_TE.keep1.TE.sig, abs(DTL24h.u5_WT.FC) < log2(1.5) | DTL24h.u5_WT.p > 0.05) # 914

rna_TE.keep1.TE.sig$type[rna_TE.keep1.TE.sig$Gene.ID %in% rna_TE.keep1.TE.dn$Gene.ID &
                           rna_TE.keep1.TE.sig$Gene.ID %in% rna_TE.keep1.rna.dn$Gene.ID] <- "TE_down.RS_down"
rna_TE.keep1.TE.sig$type[rna_TE.keep1.TE.sig$Gene.ID %in% rna_TE.keep1.TE.dn$Gene.ID &
                           rna_TE.keep1.TE.sig$Gene.ID %in% rna_TE.keep1.rna.up$Gene.ID] <- "TE_down.RS_up"
rna_TE.keep1.TE.sig$type[rna_TE.keep1.TE.sig$Gene.ID %in% rna_TE.keep1.TE.dn$Gene.ID &
                           rna_TE.keep1.TE.sig$Gene.ID %in% rna_TE.keep1.rna.ns$Gene.ID] <- "TE_down.RS_nc"

rna_TE.keep1.TE.sig$type[rna_TE.keep1.TE.sig$Gene.ID %in% rna_TE.keep1.TE.up$Gene.ID &
                           rna_TE.keep1.TE.sig$Gene.ID %in% rna_TE.keep1.rna.dn$Gene.ID] <- "TE_up.RS_down"
rna_TE.keep1.TE.sig$type[rna_TE.keep1.TE.sig$Gene.ID %in% rna_TE.keep1.TE.up$Gene.ID &
                           rna_TE.keep1.TE.sig$Gene.ID %in% rna_TE.keep1.rna.up$Gene.ID] <- "TE_up.RS_up"
rna_TE.keep1.TE.sig$type[rna_TE.keep1.TE.sig$Gene.ID %in% rna_TE.keep1.TE.up$Gene.ID &
                           rna_TE.keep1.TE.sig$Gene.ID %in% rna_TE.keep1.rna.ns$Gene.ID] <- "TE_up.RS_nc"

table(rna_TE.keep1.TE.sig$type)

ggplot(rna_TE.keep1.TE.sig, aes(log2FC_TE.u5_WT, DTL24h.u5_WT.FC)) +
  geom_point(aes(color = type))+ ylim(c(-4, 4)) +
  scale_size(range = c(0, 4)) +
  scale_color_manual(limits = c('TE_down.RS_down','TE_down.RS_up','TE_down.RS_nc', 'TE_up.RS_up', 'TE_up.RS_down', 'TE_up.RS_nc'),
                     values = c("#b4deec", "#7ab2d4", "#196eaf", "#f9e469" ,"#fdb96b", "#df352f")) +
  theme(panel.grid = element_blank(), axis.text = element_text(color="black"), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.title = element_blank(), legend.background = element_blank(), legend.key = element_blank(),
        legend.position =c(0.5,0.88) 
  ) + guides(color=guide_legend(ncol = 2, byrow=F)) +
  xlab(expression(Log[2]*"FC(TE)")) + ylab(expression(Log[2]*"FC(RS)"))

#ggsave("Figure3C.pdf", width = 10, height = 10, units = "cm")

## venn diagram: DEGs and changed TE&Proteome in u5 mutants, Figure3D ####
library(VennDiagram)

# DTL24h RNAseq
rna <- read.xlsx("../2_DTL24h_mRNAseq/DTL24h_WT_u5.merged_info.xlsx")
u5.DEGs <- subset(rna, abs(DTL24h.u5_WT.FC) > log2(1.5) & DTL24h.u5_WT.p < 0.05)
u5.DEGs.up <- subset(u5.DEGs, DTL24h.u5_WT.FC > 0) # 460
u5.DEGs.down <- subset(u5.DEGs, DTL24h.u5_WT.FC < 0) # 480

# Translation efficiency
FPKM.s.TE <- read.xlsx("../3_Polysome_seq/PolysomeSeq_TE_info.xlsx")
# significantly changed TE
u5.TE <- subset(FPKM.s.TE, pvalue < 0.05)
u5.TE.down <- subset(u5.TE, log2FC_TE.u5_WT < 0)
u5.TE.up <- subset(u5.TE, log2FC_TE.u5_WT > 0)

# TMT-based proteome results
prot <- read.xlsx("../4_Proteome/DTL24h_TMT_proteome.xlsx")
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
fpkm.DTL24 <- read.xlsx("../2_DTL24h_mRNAseq/DTL24h_WT_u5.merged_info.xlsx")
fpkm.DTL24$fpkm.ratio <- (fpkm.DTL24$FPKM.DTL24h.u5+1e-5)/(fpkm.DTL24$FPKM.DTL24h.WT+1e-5)
fpkm.DTL24$fpkm.log2ratio <- log(fpkm.DTL24$fpkm.ratio, 2)

# TMT proteome log2FC
prot <- read.xlsx("../3_Proteome/DTL24h_TMT_proteome.xlsx")
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
#ggsave("Figure3H.pdf", width = 10, height = 10, units = "cm")

## Heatmap: time-series expression profiles of transcription and translation, Figure3N&O ####
library(openxlsx)
library(stringr)
library(ComplexHeatmap)
library(circlize)

# TMT and TE both downregulted in u5
g242 <- read.table("TE_TMT_down_g242.genelist.txt")
colnames(g242) <- "Gene.ID"

# read proteome results from Pipitone et al. (2021)
norm.Avg1.anno <- read.xlsx("../4_Proteome/eLife_all_normalized_proteome.xlsx")

g242.norm <- norm.Avg1.anno[norm.Avg1.anno$Gene.ID %in% g242$Gene.ID, ] #188

#keep records with peptides >=2 
g242.norm <- subset(g242.norm, Peptide.count >=2) #188

g242.rela <- g242.norm[,c("Gene.ID",str_subset(colnames(g242.norm), "rela"))]

colnames(g242.rela)[-1] <- c("0","4","8","12","24","48","72","96")
g242.rela <- na.omit(g242.rela)
# only keep 0h, 8h, 12h, 24h, 48h
g242.rela <- g242.rela[,c("Gene.ID", "0", "8", "12", "24", "48")]
rownames(g242.rela) <- g242.rela$Gene.ID

## Proteome Heatmap
set.seed(123)
p.g242.prot <- g242.rela[,2:6]  %>%  Heatmap(
  cluster_rows = T, 
  #col = colorRamp2(c(0, 3), c("lightblue1", "blue")),
  col = colorRamp2(c(0,1,2,2.5,3,4), 
                   c("#EEEEEE","#E7EDDC","#F6E4AB","#FBA35C","#F06640","#D42D26")),
  cluster_columns = F,  column_title = "Proteome",
  #split = g242.rela$group,
  show_row_dend = F,
  #row_order = row_order(p.g108.fpkm.rela), 
  na_col = "darkgrey", name = "Ratio",
  column_names_rot = 0, column_names_gp = gpar(fontsize=11),
  column_names_side = "bottom",
  border = F,
  show_row_names = F,row_names_gp = gpar(fontsize=4),
  #width = unit(7, "cm"),
  #height = unit(9, "cm")
)
p.g242.prot1 <- draw(p.g242.prot)


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

# keep genes in proteome results
fpkm.g242.rela <- fpkm.all.rela[fpkm.all.rela$Gene.ID %in% g242.norm$Gene.ID,]
colnames(fpkm.g242.rela)[2:6] <- c("0", "6", "12", "24", "48")

## RNAseq Heatmap
set.seed(123)
p.g242.fpkm.rela <- fpkm.g242.rela[,2:6]  %>%  Heatmap(
  cluster_rows = T, show_row_dend = F,
  col = colorRamp2(c(0,1,2,2.5,3,4), 
                   c("#EEEEEE","#E7EDDC","#F6E4AB","#FBA35C","#F06640","#D42D26")),
  cluster_columns = F, column_title = "RNAseq",
  #row_order = row_order(p.g242.prot1),
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
p.g242.fpkm.rela1 <- draw(p.g242.fpkm.rela)

gb.p.g242.fpkm.rela = grid.grabExpr(draw(p.g242.fpkm.rela))
gb.p.g242.prot.rela = grid.grabExpr(draw(p.g242.prot))
cowplot::plot_grid(gb.p.g242.fpkm.rela, gb.p.g242.prot.rela, ncol=2)
#ggsave("Figure3N&O.pdf", width=15, height = 10, units = "cm")


## Heatmap for U5 targets: time-series RNAseq and TMT_proteome, Figure S2C&D ####
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
#pdf("FigureS2C.pdf", height = 5/2.54, width = 12/2.54)
p
#dev.off()


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
#pdf("FigureS2D.pdf", height = 5/2.54, width = 12/2.54)
p
#dev.off()
