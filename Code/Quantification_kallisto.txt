# command line RNA-Seq

# on raw data
module load fastqc

fastqc *fastq.gz


module load multiqc

multiqc -n rnaseq



# trim files
module load trimgalore

trim_galore --paired --phred33 --fastqc A1_R1.fastq.gz A1_R2.fastq.gz # etc



# kallisto index
module load kallisto

kallisto index -i AipIndex.idx aip.genome_models.no_isoforms.no_duplication.mRNA.fa



# kallisto quantification
module load kallisto

kallisto quant -i AipIndex.idx -o ./quant/A1 -b 100 -t 8 ./trimmed/A1_R1.fastq.gz ./trimmed/A2_R1.fastq.gz # etc
