# Recurrent stroke R15 genomics pipeline

Michael Brewer
2018-11-14

*_running on @XXX.25_*

```
mkdir 2018_R15_NGS_data
mkdir hg38_refs
mkdir appreci8_stuff
```

Downloaded from Vanderbuilt Box.com account via Firefox

## QC

```
for FileName in *R1_001.fastq.gz
  do name=`ls $FileName | cut -d"_" -f1,2`
  fastp -i $name"_R1_001.fastq.gz" -I $name"_R2_001.fastq.gz" -o $name"_R1_001_FASTP.fastq.gz" -O $name"_R2_001_FASTP.fastq.gz" -z 4 -q 20 -c
done
```

### Resynchronizing paired-end files

**NOT NECESSARY; REMOVED**

```
#for FileName in QCed/*R1*.gz
#  do name=`ls $FileName | cut -d"_" -f1,2`
#  fastqCombinePairedEnd.py $name"_R1_001_FASTP.fastq.gz" $name"_R2_001_FASTP.fastq.gz"
#done
```

```
mkdir QCed
mv *FASTP.fastq.gz QCed/
```

## Read mapping

### Prepare databases for mapping

```
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
tar -xvf Homo_sapiens_UCSC_hg38.tar.gz
```

### Mapping

**bwa-mem_mapper.sh**

```
#!/bin/bash

name=`ls $1 | cut -d"_" -f1,2`

header=$(zcat $1 | head -n 1)
id=$(echo $name | cut -d"_" -f1 | cut -d"/" -f2)
#id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
sm=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+$")
echo "Read Group @RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA"

bwa mem \
-M \
-t 8 \
-R $(echo "@RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA") \
~/hg38_refs/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa \
$1 $2 | samtools sort -O bam -o $name"_mapped-bwa.bam" -
```

```
#for FileName in QCed/*R1*.gz; do name=`ls $FileName | cut -d"_" -f1,2`
#  bwa mem -M -R {FIGURE THIS OUT} -t 8 #../hg38_refs/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa $name"_R1_001_FASTP.fastq.gz" #$name"_R2_001_FASTP.fastq.gz" | samtools sort -o $name"_bwa-aln.bam" -
#done
```

```
for FileName in QCed/*R1*.gz; do name=`ls $FileName | cut -d"_" -f1,2`
  bwa-mem_mapper.sh $name"_R1_001_FASTP.fastq.gz" $name"_R2_001_FASTP.fastq.gz"
done
```

### Filter for unique mappings and removing duplicates

```
for FileName in QCed/*.bam; do
  name=`ls $FileName | cut -d"." -f1`
  sambamba view -t 8 -h -f bam -F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)" $FileName -o $name"_unique.bam"
  sambamba markdup -r -t 8 $name"_unique.bam" $name"_unique_DUPremoved.bam"
done
```

### Mark/Remove Duplicates

**NOT NEEDED; REMOVED**

```
#for FileName in *unique.bam; do
#  java -jar MarkDuplicates.jar REMOVE_DUPLICATES=true METRICS_FILE=$FileName"_dup.txt" INPUT=$FileName OUTPUT=$FileName"_DUPremoved.bam"
#done
```

## Variant calling

### GATK4

*Haplotypecaller*

Getting files in place

```
cp ~/hg38_refs/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.* ~/2018_R15_NGS_data/
```

Launch GATK4 in Docker

```
docker run -v ~/2018_R15_NGS_data:/gatk/my_data -it broadinstitute/gatk:latest
```

Run GATK4 HaplotypeCaller to call variants

**In GATK4 Docker**

```
for FileName in my_data/QCed/*DUPremoved.bam; do
  java -jar gatk-package-4.0.11.0-local.jar HaplotypeCaller -ERC GVCF -R my_data/variant_calling/genome.fa -I $FileName -O $FileName"_GATKvariants.vcf"
done

mkdir my_data/gvcfs
cp my_data/QCed/*.vcf my_data/QCed/*.idx my_data/gvcfs
```

Select only SNPS (no indels or MNPs)

**In GATK4 Docker**

```
#for FileName in my_data/gvcfs/*.vcf; do
#  java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -jar gatk-package-4.0.11.0-local.jar SelectVariants -R my_data/genome.fa -V $FileName --output $FileName"_onlySNPs.vcf" --select-type-to-include SNP
#done
```

Combine and Genomtype individual gVCFs in to final VCF

**In GATK4 Docker**

```
cp my_data/CombineGVCFs_run* .
bash CombineGVCFs_run_Chr1.sh
bash CombineGVCFs_run_Chr2.sh

#java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -jar gatk-package-4.0.11.0-local.jar GenotypeGVCFs -R my_data/genome.fa -V gendb:my_data/gvcfs/gworkspace -G StandardAnnotation --use-new-qual-calculator -O my_data/1707-KK-Combined1-95_Genotyped.vcf

java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -jar gatk-package-4.0.11.0-local.jar GenotypeGVCFs -R my_data/genome.fa -V gendb://my_data/gvcfs/gworkspace_Chr1 -G StandardAnnotation --use-new-qual-calculator -O my_data/1707-KK-Combined1-95_Gentotyped_Chr1.vcf

java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -jar gatk-package-4.0.11.0-local.jar GenotypeGVCFs -R my_data/genome.fa -V gendb://my_data/gvcfs/gworkspace_Chr2 -G StandardAnnotation --use-new-qual-calculator -O my_data/1707-KK-Combined1-95_Gentotyped_Chr2.vcf
```
