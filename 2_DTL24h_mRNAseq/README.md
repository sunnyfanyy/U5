# RNAseq of WT and mutants grown under DTL-24h

The raw sequencing data of WT and *u5-3 u5-4* double mutant have been deposited at the Gene Expression Omnibus (GEO) with the accession number [GSE266658](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266658).

The raw RNA-seq data of Sm protein mutants (*smd3a* mutants and *smd3b* mutants) to GEO under accession number [GSE283246](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283246). 

## Quality control and preprocessing of raw sequence data

```sh
DIR=fastp_report
fastp -i Cleandata/$i/${i}_R1.fq.gz \
	-o fastq/${i}_R1.fq.gz \
	-I Cleandata/$i/${i}_R2.fq.gz \
	-O fastq/${i}_R2.fq.gz \
	-l 20 --thread 15 -h ${DIR}/${i}_fastp.html
```

## Mapping to *Arabidipsis thaliana* TAIR10 genome

```sh
# bulid TAIR10 genome index with hisat2
hisat2-build -p 6 TAIR10.fa --exon tair10.exon --ss tair10.splice_sites TAIR_genome_tran

INDEX=database/tair10/TAIR_genome_tran
hisat2 -p 18 --dta -x ${INDEX} \
	--min-intronlen 50 --max-intronlen 6000 \
	-1 ${id}_R1.fq.gz -2 ${id}_R2.fq.gz -S ${name}.sam
# keep the unique and high mapping quality reads
    grep -E '^@|NH:i:1' ${name}.sam | samtools view -bhS -@ 18 -q 20 | samtools sort -@ 18 -o ${DIR}/highQ_bam/${name}.highQ.bam
    samtools index ${DIR}/highQ_bam/${name}.highQ.bam
```

## Read summarization for genes with featureCounts

```sh
GTF=Araport11.gtf
featureCounts -T 15 -p -t exon -g gene_id -O -M -a ${GTF} \
  -o ${OUT}/sample_gene_count_20210324.txt \
  MAPPING_RESULT/highQ_bam/*.bam
  
```

The `sample_gene_count_20210324.txt` file can be used to calculate gene expression levels and can also serve as an input file for DESeq2 software to identify differentially expressed genes in *u5-3 u5-4* double mutants.

## Differential Alternative Splicing Analysis with rMATS

WT: Col-0 ecotype; cn1cn2: *u5-3 u5-4* double mutant; D3a: *smd3a* mutant; D3b: *smd3b* mutant.

```sh
GTF=Araport11.gtf
##### DTL24h ###
echo WT vs. cn1cn2
python ${dir}/rmats.py --b1 bam_DTL24h_WT.txt --b2 bam_DTL24h_cn1cn2.txt \
  --gtf ${GTF} --od rMATS/DTL24h_WT_cn1cn2 -t paired --readLength 150 --libType fr-unstranded --nthread 15 --tstat 15 --anchorLength 6
  
echo 24h-WT vs. 24h-D3a
python ${dir}/rmats.py --b1 24h-WT.txt --b2 24h-D3a.txt \
  --gtf ${GTF} --od rMATS/24h-WT_D3a --tmp rMATS/24h-WT_D3a  -t paired --readLength 150 --libType fr-unstranded --nthread 15 --tstat 15 --anchorLength 6

echo 24h-WT vs. 24h-D3b
python ${dir}/rmats.py --b1 24h-WT.txt --b2 24h-D3b.txt \
  --gtf ${GTF} --od rMATS/24h-WT_D3b --tmp rMATS/24h-WT_D3b  -t paired --readLength 150 --libType fr-unstranded --nthread 15 --tstat 15 --anchorLength 6

```

*All downstream analyses were implemented in the R script `mRNAseq_analysis.R`.*
