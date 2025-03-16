# Polysome-seq of WT and *u5-3 u5-4* double mutant grown under DTL-24h

The raw sequencing data have been deposited at the Gene Expression Omnibus (GEO) with the accession number [GSE266868](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266868).

## Quality control and preprocessing of raw sequence data

```sh
DIR=fastp_report
fastp -i Rawdata/$i/${i}_R1.fq.gz \
	-o fastq/${i}_R1.fq.gz \
	-I Rawdata/$i/${i}_R2.fq.gz \
	-O fastq/${i}_R2.fq.gz \
	--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
	--adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	-l 20 --thread 16 -h ${DIR}/${i}_fastp.html
```

## rRNA reads filteration and mapping to TAIR10 genome 

To exclude rRNA contamination, Bowtie2 was employed to remove reads mapped to rRNA sequences. 

```sh
DIR=rRNA_depleted
#bowtie2 v2.3.5
INDEX=database/tair10/bowtie2/rRNA_bowtie2
bowtie2 -p 16 -x ${INDEX} \
 --un-conc-gz ${DIR}/${i} -1 fastq/${i}_R1.fq.gz -2 fastq/${i}_R2.fq.gz -S ${i}.sam
```

The remaining reads were subsequently mapped to the Arabidopsis thaliana (TAIR10) reference genome.

```sh
DIR=MAPPING_RESULT_rRNA_depleted
INDEX=database/tair10/TAIR_genome_tran
hisat2 -p 18 --dta -x ${INDEX} \
	--min-intronlen 50 --max-intronlen 6000 \
	-1 rRNA_depleted/${id}.1.fq.gz -2 rRNA_depleted/${id}_R2.fq.gz -S ${id}.sam
# keep the unique and high mapping quality reads
    grep -E '^@|NH:i:1' ${id}.sam | samtools view -bhS -@ 18 -q 20 | samtools sort -@ 18 -o ${DIR}/highQ_bam/${id}.highQ.bam
    samtools index ${DIR}/highQ_bam/${id}.highQ.bam
```

## Read summarization for genes with featureCounts

```sh
GTF=Araport11.gtf
featureCounts -T 15 -p -s 2 -t exon -g gene_id -O -M -a ${GTF} \
  -o ${OUT}/PolysomeSeq.raw_counts.txt \
 MAPPING_RESULT_rRNA_depleted/highQ_bam/*.bam
```

## Translation efficiency (TE) calculation

Translation efficiency (TE) was derived by normalizing polysome-associated RNA FPKM to cytoplasmic RNA FPKM.

*All downstream analyses were implemented in the R script `PolysomeSeq_analysis.R`.*

