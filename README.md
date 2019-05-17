# GeneQC: Gene expression level Quality Control #

## Abstract ##

One of the main benefits of using modern RNA-Sequencing (RNA-Seq) technology is the more accurate gene expression estimations compared with previous generations of expression data, such as the microarray. However, numerous issues can result in the possibility that an RNA-Seq read can be mapped to multiple locations on the reference genome with the same alignment scores, which occurs in plant, animal, and metagenome samples. Such a read is so-called a multiple-mapping read (MMR). The impact of these MMRs is reflected in gene expression estimation and all downstream analyses, including differential gene expression, functional enrichment, etc. Current analysis pipelines lack the tools to effectively test the reliability of gene expression estimations, thus are incapable of ensuring the validity of all downstream analyses. Our investigation into 95 RNA-Seq datasets from seven plant and animal species (totaling 1,951 GB) indicates an average of roughly 22% of all reads are MMRs. Here we present a machine learning-based tool called GeneQC (Gene expression Quality Control), which can accurately estimate the reliability of each gene’s expression level derived from an RNA-Seq dataset. The underlying algorithm is designed based on extracted genomic and transcriptomic features, which are then combined using elastic-net regularization
and mixture model fitting to provide a clearer picture of mapping uncertainty for each gene. GeneQC allows researchers to determine reliable expression estimations and conduct further analysis on the gene expression that is of sufficient quality. This tool also
enables researchers to investigate continued re-alignment methods to determine more accurate gene expression estimates for those with low reliability. Application of GeneQC reveals high level of mapping uncertainty in plant samples and limited, severe mapping uncertainty in animal samples. GeneQC is freely available at http://bmbl.sdstate.edu/GeneQC/home.html.


**Citing us:** McDermaid, A., Chen, X., Zhang, Y., Wang, C., Gu, S., Xie, J., & Ma, Q. (2018). A New Machine Learning-Based Framework for Mapping Uncertainty Analysis in RNA-Seq Read Alignment and Gene Expression Estimation. *Frontiers in Genetics*, 9. doi:10.3389/fgene.2018.00313

## Environment ##

The GeneQC package requires Python 3 to execute, includes the blast+ and SAMtools libraries, GeneQC takes Reference Genome, Annotation file, Read alignment (mapped by HISAT2, compressed to BAM file by SAMtools) as input, you may take the following steps to obtain the preliminary files. The GeneQC will generate feature extraction and modeling results(D-score).

The package contains:

1. Minimum package of Blast for makeblastdb and Blastn

2. Samtools package version 1.2.1

3. GeneQC python code

4. A special script "extract_transcript_seq_gff.py" which used to transform the genome sequences to transcript sequences and the corresponding annotation file.

**Hardware Requirement:**

Hardware requirements are based on the size of Reference Genome and Annotation file. High performance cluster preferred.

**Software Requirement:**

GeneQC package requires Python 3 to execute. We recommended user use anaconda 5 to execute the GeneQC. Most of high performance clusters have already installed the anaconda, so user can load anaconda directly using following code.

1. Check available modules in the cluster:
```{r,engine='bash',eval=FALSE}
module avail
```
2. find anaconda (4 or above) in available modules and load it. (for example, if anaconda5/5.0.0-3.6 is in the available modules list):
```{r,engine='bash',eval=FALSE, module}
module load anaconda5/5.0.0-3.6
```

# Usage #

## Installation
The source code of GeneQC is freely available at: https://github.com/maqin2001/GeneQC

To install GeneQC, first download the zip file manually from github, or use the code below in Unix:
```{r,engine='bash',eval=FALSE, download}
cd your_folder_path
wget https://github.com/maqin2001/geneqc/archive/master.zip
```
Unzip the file:
```{r,engine='bash',eval=FALSE, unzip}
unzip master.zip
```

## Input data preperation

Seven sample species of data (Reference Genome, Annotation file, Read alignment (mapped by HISAT2, compressed to BAM file by SAMtools) can be downloaded from our website (http://bmbl.sdstate.edu/GeneQC/result.html).

For example, here we want to download RNA-Seq data of Homo sapines specie

```{r,engine='bash',eval=FALSE}
cd your_folder_path_of_"GeneQC_Python"
wget http://bmbl.sdstate.edu/GeneQC/GeneQC_files/data/Homo_sapiens/Humo_raw_data.tar.gz
gunzip -c Humo_raw-data.tar.gz | tar xopf -
```
## Instructions (Plant Genome)

For running GeneQC for plant data: Three inputs (data) are required: reference genome, annotation file, sam or bam file should be uploaded to folder of "GeneQC_Python" under the GeneQC-master folder in the cluster.

Move to the path of folder of "GeneQC_Python"
```{r,engine='bash',eval=FALSE}
cd your_folder_path_of_"GeneQC_Python"
module load anaconda(chose_your_version)
```

Run GeneQC:
```{r,engine='bash',eval=FALSE}
python GeneQC2.py [1] [reference genome] [standard gff annotation file] [sam or bam file]
```

Example: A.thaliana
```{r,engine='bash',eval=FALSE}
python GeneQC2.py 1 Athaliana_167_TAIR9.fa Athaliana_167_TAIR10.gene.gff3 ERR1297323.bam
```

The outputs will be generated in this folder as well. ERR1297323_out.txt will be feature extraction results. ERR1297323_out.csv will be D-scoure results.

## Instructions (Amimal Transcript)

For running GeneQC for animal data: Three inputs (data) are required: reference genome, annotation file, fastq file should be uploaded to folder of "GeneQC_Python" under the GeneQC-master folder in the clust

Move to the path of folder of "GeneQC_Python"
```{r,engine='bash',eval=FALSE}
cd your_folder_path_of_"GeneQC_Python"
module load anaconda(chose_your_version)
```

Step1: Create new defined transcripts and new defined transcripts annotation:
```{r,engine='bash',eval=FALSE}
python extract_transcript_seq_gff.py [reference genome] [standard gff annotation file] [new defined transcripts sequence file] [new defined transcripts gff annotation file]
```

Example: Humo sapiens
```{r,engine='bash',eval=FALSE}
python extract_transcript_seq_gff.py GCF_000001405.37_GRCh38.p11_genomic.fna GCF_000001405.37_GRCh38.p11_genomic.gff human_transcripts_seq.fa human_transcripts_seq.gff
```

### For Bulk RNA-seq data:
Step2: Do RNA-seq mapping work with the new mapping results, use following commands (the example used aligner HISAT2):
```{r,engine='bash',eval=FALSE}
module load hisat2
hisat2-build -f human_transcripts_seq.fa ./hisatindex/Humo
hisat2 -x ./hisatindex/Humo -k 10 -p 40 -1 SRR6029567_1.fastq -2 SRR6029567_2.fastq -S SRR6029567.sam
```

Step3: Run GeneQC:
```{r,engine='bash',eval=FALSE}
python GeneQC2.py 2 human_transcripts_seq.fa human_transcripts_seq.gff SRR6029567.sam
```

### For Single cell RNA-seq data:
Step2: Do RNA-seq mapping work with the new mapping results, use following commands (the example used aligner HISAT2):
```{r,engine='bash',eval=FALSE}
module load hisat2
hisat2-build -f human_transcripts_seq.fa ./hisatindex/Humo
hisat2 -x ./hisatindex/Humo -k 10 SRR491087.fastq -S SRR491087.sam
```

Step3: Run GeneQC:
```{r,engine='bash',eval=FALSE}
python GeneQC2.py 2 human_transcripts_seq.fa human_transcripts_seq.gff SRR491087.sam
```

The outputs will be generated in this folder as well. SRR6029567_out.txt will be feature extraction results. SRR6029567_out.csv will be D-scoure results.
