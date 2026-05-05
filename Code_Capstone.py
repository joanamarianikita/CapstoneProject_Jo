Mdkir capstone

cd capstone
cd case_3
conda activate 1.readq	
fastqc case3_R1.fastq.gz
fastqc case3_R1.fastq.gz

trimmomatic PE \
case3_R1.fastq.gz case3_R2.fastq.gz \
case3_R1_paired.fastq.gz case3_R1_unpaired.fastq.gz \
case3_R2_paired.fastq.gz case3_R2_unpaired.fastq.gz \
ILLUMINACLIP:TruSeq3-SE.fa:2:30:10:2:True LEADING:30 TRAILING:30 MINLEN:36 HEADCROP:15 TAILCROP:5

fastqc case3_R1_paired.fastq.gz case3_R1_unpaired.fastq.gz
fastqc case3_R2_paired.fastq.gz case3_R2_unpaired.fastq.gz

bwa mem -t 1 \
-R '@RG\tID:a1\tSM:nor\tPL:ILLUMINA\tLB:lib1\tPU:f_1' \
~/capstone/reference/Homo_sapiens_assembly38.fasta \
~/capstone/case_3/case3_R1_paired.fastq.gz \
~/capstone/case_3/case3_R2_paired.fastq.gz \
> ~/capstone/outfile/sample_case3.sam

samtools view -Shb \
-o ~/capstone/outfile/sample_case3.bam \
~/capstone/outfile/sample_case3.sam

-o ~/capstone/outfile/sample3-sorted.bam \
~/capstone/outfile/germline.bam

samtools index \
~/capstone/outfile/sample_case3-sorted.bam

picard MarkDuplicates \
INPUT=$HOME/capstone/outfile/sample_case3-sorted.bam \
METRICS_FILE=$HOME/capstone/outfile/sample_case3-met.txt \
OUTPUT=$HOME/capstone/outfile/sample_case3-dedup.bam \
CREATE_INDEX=true

gatk BaseRecalibrator \
-R $HOME/capstone/reference/Homo_sapiens_assembly38.fasta \
-I $HOME/capstone/outfile/sample_case3-dedup.bam \
--known-sites $HOME/capstone/reference/Homo_sapiens_assembly38.known_indels.vcf.gz \
-O $HOME/capstone/outfile/sample_case3.grp

gatk ApplyBQSR -R ~/capstone/reference/Homo_sapiens_assembly38.fasta \
-I  ~/capstone/outfile/sample_case3-dedup.bam \
-bqsr ~/capstone/outfile/sample_case3.grp \
-O ~/capstone/outfile/sample_case3-bqsr.bam

gatk CollectAlignmentSummaryMetrics \
> -R ~/capstone/reference/Homo_sapiens_assembly38.fasta \
>--INPUT ~/capstone/outfile/sample_case3-bqsr.bam \
>--OUTPUT ~/capstone/outfile/sample_case3-bqsr_Alignment_Summary2.txt

gatk CollectInsertSizeMetrics \
--INPUT ~/capstone/outfile/sample_case3-bqsr.bam \
--OUTPUT ~/capstone/outfile/sample_case3-bqsr_Size_Metrics.txt \
--Histogram_FILE ~/capstone/outfile/sample_case3-bqsr_Size_Histo.pdf

cd capstone
cd reference
conda activate 3.varcall

gatk HaplotypeCaller \
-L chr13 \
-L chr17 \
-R ~/capstone/reference/Homo_sapiens_assembly38.fasta \
-I ~/capstone/outfile/sample_case3-sorted.bam \
-O ~/capstone/outfile/sample_case3.vcf

gatk VariantFiltration \
-R ~/capstone/reference/Homo_sapiens_assembly38.fasta \
-V ~/capstone/outfile/sample_case3.vcf \
-O ~/capstone/outfile/sample_case3_filtered.vcf \
--filter-name "QD2" --filter-expression "QD < 2.0" \
--filter-name "FS60" --filter-expression "FS > 60.0" \
--filter-name "MQ40" --filter-expression "MQ < 40.0"

bcftools view \
-f PASS \
-O z \
-o ~/capstone/outfile/sample_case3_pass.vcf.gz \
~/capstone/outfile/sample_case3_filtered.vcf

bgzip ~/capstone/outfile/sample3_filtered.vcf
tabix -p vcf ~/capstone/outfile/sample_case3_filtered.vcf.gz

bcftools annotate \
-a ~/capstone/reference/clinvar_chr.vcf.gz \
-c ID,INFO \
-h <(bcftools view -h ~/capstone/reference/clinvar_chr.vcf.gz | grep "^##INFO=") \
~/capstone/outfile/sample3_filtered.vcf.gz \
-Oz -o ~/capstone/outfile/sample_case3_filtered_annotated.vcf.gz

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ALLELEID\t%INFO/RS\t%INFO/CLNHGVS\t%INFO/GENEINFO\t%INFO/MC\t%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/ONC\t%INFO/SCI\t%INFO/SCIDN\n' sample_case3_filtered_annotated.vcf.gz > sample_case3_filtered_annotated_extracted.tsv

bcftools view \
-i 'INFO/CLNSIG = "Pathogenic" || INFO/CLNSIG = "Likely_pathogenic"' \
~/capstone/outfile/sample_case3_filtered_annotated.vcf.gz \
-Oz -o ~/capstone/outfile/BRCA_case3_pathogenic.vcf.gz




