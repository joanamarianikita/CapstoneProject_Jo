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
