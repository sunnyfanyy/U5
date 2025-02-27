#### 1. elife normalized abundance ####
library(openxlsx)
library(tidyr)
library(stringr)

proteome <- read.xlsx("elife-62709-fig5-data1-v1 (1).xlsx",
                      sheet = "normalized_abundance", startRow=2)
View(proteome)
colnames(proteome)[c(6:9, 10:14, 
                     15:19, 20:24, 25:29, 30:33, 34:38,
                     39:43)] <- c("Rep1_0h","Rep2_0h","Rep3_0h","Avg_0h",
                                  "Rep1_4h","Rep2_4h","Rep3_4h","Rep4_4h","Avg_4h",
                                  "Rep1_8h","Rep2_8h","Rep3_8h","Rep4_8h","Avg_8h",
                                  "Rep1_12h","Rep2_12h","Rep3_12h","Rep4_12h","Avg_12h",
                                  "Rep1_24h","Rep2_24h","Rep3_24h","Rep4_24h","Avg_24h",
                                  "Rep1_48h","Rep2_48h","Rep3_48h","Avg_48h",
                                  "Rep1_72h","Rep2_72h","Rep3_72h","Rep4_72h","Avg_72h",
                                  "Rep1_96h","Rep2_96h","Rep3_96h","Rep4_96h","Avg_96h")
a.sp <- strsplit(proteome$Accession, "[;]")

#将Majority.protein.IDs分别和拆分后的列表的每个元素合并
a.spc <- mapply(cbind, proteome$Accession, a.sp)

#最后将列表a.spc转换为数据框，按行合并，即可得到目标排列的数据
a <- do.call(rbind.data.frame, a.spc)
colnames(a) <- c("Accession", "Locus.ID")

a <- separate(a, Locus.ID, into = c("Gene.ID","dup"), sep = "[.]")
a$dup <- NULL
#a <- a[!duplicated(proteome$Gene.ID),]

b <- left_join(a, proteome, by="Accession")
b <- b[!duplicated(b),] # 5678

#读取Araport注释
anno <- read.xlsx("../Polysome_seq/Araport11_annotation.xlsx")
anno$computational_description <- NULL
anno$locus_type <- NULL #由于下面是质谱的结果，所以都是protein coding

norm <- left_join(b, anno, by=c("Gene.ID"="ID"))

norm <- arrange(norm, Gene.ID, location_consensus)
norm <- norm[!duplicated(norm$Gene.ID),]
norm.Avg1 <- norm[,c(colnames(norm)[1:6], str_subset(colnames(norm), "Avg"))]

norm.Avg1$Avg_0h.rela <-  norm.Avg1$Avg_0h/norm.Avg1$Avg_0h
norm.Avg1$Avg_4h.rela <-  norm.Avg1$Avg_4h/norm.Avg1$Avg_0h
norm.Avg1$Avg_8h.rela <-  norm.Avg1$Avg_8h/norm.Avg1$Avg_0h
norm.Avg1$Avg_12h.rela <-  norm.Avg1$Avg_12h/norm.Avg1$Avg_0h
norm.Avg1$Avg_24h.rela <-  norm.Avg1$Avg_24h/norm.Avg1$Avg_0h
norm.Avg1$Avg_48h.rela <-  norm.Avg1$Avg_48h/norm.Avg1$Avg_0h
norm.Avg1$Avg_72h.rela <-  norm.Avg1$Avg_72h/norm.Avg1$Avg_0h
norm.Avg1$Avg_96h.rela <-  norm.Avg1$Avg_96h/norm.Avg1$Avg_0h

norm.Avg1.anno <- left_join(norm.Avg1, norm[,c(2,45:49)],by="Gene.ID")
#write.xlsx(norm.Avg1.anno, "eLife_all_normalized_proteome.xlsx")

