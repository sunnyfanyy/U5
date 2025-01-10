# 1. Subcellular RNAseq

The raw sequencing data have been deposited at the Gene Expression Omnibus (GEO) with the accession number [GSE266659](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266659).

## 1.1 QC of raw reads with fastqc

Quality control of raw reads was carried out using fastqc.  

```sh
# threads option (-t below) may need to be adjusted for your machine.
 fastqc Raw/${id}.1.clean.fq.gz Raw/${id}.2.clean.fq.gz -t 20 -o fastqc
```

## 1.2 remove adaptor with cutadapt

```sh
R1=Raw/${id}.1.clean.fq.gz
R2=Raw/${id}.2.clean.fq.gz
cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -m 20 \
 -o Raw/cutadapt/${id}.1.fq.gz -p Raw/cutadapt/${id}.2.fq.gz ${R1} ${R2}
```

## 1.3 remove rRNA reads with bowtie2

```sh
# build TAIR10 index with bowtie2
bowtie2-build Arabidopsis_rRNA.fasta database/tair10/bowtie2/rRNA_bowtie2

INDEX=tair10/bowtie2/rRNA_bowtie2
R1=Raw/cutadapt/${id}.1.fq.gz
R2=Raw/cutadapt/${id}.2.fq.gz

bowtie2 -p 10 -x ${INDEX} \
 --un-conc-gz Raw/cutadapt_bowtie2_rRNA_depleted/${id} -1 ${R1} -2 ${R2} -S ${id}.sam
```

## 1.4 mapping to *Arabidipsis thaliana* TAIR10 genome

```sh
# bulid TAIR10 genome index with hisat2
hisat2-build -p 10 TAIR10.fa database/tair10/hisat2/TAIR_genome
INDEX=tair10/hisat2/TAIR_genome
DIR=hisat2_MAPPING

hisat2 -p 20 -x ${INDEX} --dta \
  --min-intronlen 50 --max-intronlen 6000 \
  -1 ${id}.1.fq.gz -2 ${id}.2.fq.gz -S ${name}.sam

# keep perfect match reads ##
grep -E '^@|XM:i:0' ${name}.sam | ${samtools} view -bhS -@ 20 | ${samtools} sort -@ 20 -o ${DIR}/${name}.bam
samtools index ${DIR}/${name}.bam

```

## 1.5 transcript assembly with stringtie 

Merge biological replicates BAM files, for example:

```sh
# merge biological replicates BAM files and create index for it
samtools merge -@ 20 merged_bam/Col-D-F.bam Col-D-1.F.bam Col-D-2.F.bam Col-D-3.F.bam
samtools index merged_bam/Col-D-F.bam
```

Assemble transcripts:

```sh
# transcript assembly
GTF="database/Araport/Araport11_all_20190918.gtf"
i=merged_bam/${id}.bam
DIR=stringtie

stringtie $i -p 20 -G ${GTF} --rf \
 -m 100 -j 5 -c 5  -o ${DIR}/${name}.gtf
 # --rf assume stranded library fr-firststrand
 # -m: minimum assembled transcript length, default:200
 # -j: minimum junction coverage (default: 1)
 # -c: minimum reads per bp coverage to consider for transcript assembly  (default: 2.5)
```

Merge GTF files into `stringtie_merged.gtf`。

```sh
ls stringtie/*gtf > mergelist.txt

GTF="database/Araport/Araport11_all_20190918.gtf"
stringtie --merge -p 20 -G ${GTF} -o ${DIR}/stringtie_merged.gtf mergelist.txt
```

Fetch new assembled transcripts and identify novel ncRNA.

 ```sh
 grep 'transcript_id "MSTRG' stringtie_merged.gtf  > new_transcript.gtf
 gtfToGenePred new_transcript.gtf tmp
 genePredToBed tmp new_transcript.bed
 ```

Novel ncRNA identification with CPC2 and COME software.

```sh
# generate FASTA file for CPC2 analysis
bedtools getfasta -fi database/tair10/TAIR10.fa \
-bed new_transcript.bed -split -s -name > new_transcript.forCPC.fa

# CPC2:
python software/CPC2-beta/bin/CPC2.py -i ./new_transcript.forCPC.fa -o ./new_transcript.CPC.result
grep 'noncoding' new_transcript.CPC.result | cut -d ':' -f 1  > new_transcript.cpc.noncoding.list

# COME
bash software/COME-master/bin/COME_main.sh new_transcript.gtf new_transcript.COME software/COME-master/bin plant plant.model
mv result.txt new_transcript_COME.result
grep 'noncoding' new_transcript_COME.result | awk '{print $1}' > new_transcript.come.noncoding.list

```

Obtain the intersection of CPC and COME results. This part is done in R.

```R
cpc <- read.table("new_transcript.cpc.noncoding.list")
come <- read.table("new_transcript.come.noncoding.list")
m <- merge(cpc,come,by="V1", all=F)
write.table(m,"new_transcript.common.noncoding.list",row.names=F,col.names=F,quote=F)
```

We generated the GTF file `merged_correct.gtf` , which including correct snRNAs loci, novel ncRNAs.

## 1.6 read summarization for genes with featureCounts

The BAM files including: Col-D-1.N.bam, Col-D-2.N.bam, Col-D-3.N.bam, Col-D-1.F.bam, Col-D-2.F.bam, Col-D-3.F.bam, Col-D-1.SM.bam, Col-D-2.SM.bam, Col-D-3.SM.bam, Col-D-1.P.bam, Col-D-2.P.bam, Col-D-3.P.bam.

```sh
featureCounts -T 20 -p -s 2 -t exon -g gene_id -O -M -a merged_correct.gtf \
 -o stringtie/merged_count.txt *.bam
```

The `merged_count.txt` file can be used to calculate gene expression levels and can also serve as an input file for DESeq2 software to analyze transcripts with changes in subcellular localization.

## 1.7 Figures in manuscript

Figure1B, Figure1C, Figure1D and FigureS1A can be reproduced by R script `analysis.R`.
