#!/bin/bash
<<"COMMENT"
The code should do the following :
Run alignments :
	- bwa mem on the input fastq files for read1 and read2 in that order
	- sort the file and add read groups
then go through preprocessing from GATK :
	- realign indels (include npm1 indels?)
	- cut ends
	- fix mate information
	- DON'T BECAUSE IT'S A SMALL TARGETED EXPERIMENT <100Mb recalibrate (include known variations in leukemia?)
=======
License
=======
This code is released under the GNU General Public License 3.0. A copy
of this license is in the LICENSE.txt file.
copyright Irina Krier 2015
This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
COMMENT

usage ()
{
  echo 'Usage : run_alignments_onefile.sh <read1.fastq> <read2.fastq> <Libname>'
  exit
}

if [ "$#" -ne 3 ]
then
  usage
fi

echo "Running alignments on" $1 "and" $2 "to create "$3".bam"

#Just in case this affects anything let's remove the adaptors :
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT "$1" "$2" -o $3"_R1.fastq.gz" -p $3"_R2.fastq.gz" -m 20 &>$3".out"

# Read groups are necessary for downstream processing
# Working with single files, not all are useful, but for batch-processing it will be
bwa mem -t 4 -M -v 2 -R "@RG\tID:$3\tSM:$3\tPL:ILLUMINA\tLB:$3" /data/genomes/Broadhs37/hs37d5.fa.gz $3"_R1.fastq.gz" $3"_R2.fastq.gz" 2>$3.out|samtools view -bSu - |samtools sort - -f $3".bam"
samtools index $3".bam"

bedfile=/data/TruSight_analysis/trusight-myeloid-amplicon-track_forBrowser.bed

# Realigning around indels
# First we have to create the targets :
java -jar -Xmx4g -Djava.io.tmpdir=tmp /data/software/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /data/genomes/Broadhs37/hs37d5.fa -I $3".bam" -known /data/test_gatk/bundle/Mills_and_1000G_gold_standard.indels.b37.vcf -known /data/test_gatk/bundle/1000G_phase1.indels.b37.vcf -o forIndelRealigner_$3".intervals" -dt NONE &>>$3.out

# Then we have to actually realign :
java -jar -Xmx4g -Djava.io.tmpdir=tmp /data/software/GenomeAnalysisTK.jar -T IndelRealigner -R /data/genomes/Broadhs37/hs37d5.fa -I $3".bam" -targetIntervals forIndelRealigner_$3".intervals" -known /data/test_gatk/bundle/Mills_and_1000G_gold_standard.indels.b37.vcf -known /data/test_gatk/bundle/1000G_phase1.indels.b37.vcf -o $3"_realigned.bam" --maxReadsForRealignment 500000 -maxInMemory 500000 &>>$3.out

# We cut the ends
/data/TruSight_analysis/scripts/remove_probes.py $bedfile $3"_realigned.bam" $3"_realigned_cutends.bam"

# We fix the mate intormation and resort the files
java -jar -Xmx4g /data/software/picard/dist/picard.jar FixMateInformation I=$3"_realigned_cutends.bam" O=$3"_realigned_cutends.srt.bam" SORT_ORDER=coordinate ADD_MATE_CIGAR=true TMP_DIR=`pwd`/tmp &>>$3.out
samtools index $3"_realigned_cutends.srt.bam"

#We mark the duplicates and create another file :
java -jar /data/software/picard/dist/picard.jar MarkDuplicates I=$3"_realigned_cutends.srt.bam" O=$3"_realigned_cutends.srt.dups.bam" M=$3"_realigned_cutends.srt.bam.metrics" &>>$3.out
samtools index $3"_realigned_cutends.srt.dups.bam"

# We DON'T create the recalibration data
#java -Xmx4g -jar /data/software/GenomeAnalysisTK.jar    -T BaseRecalibrator    -R /data/genomes/Broadhs37/hs37d5.fa   -I $3"_realigned_cutends.srt.bam" -knownSites /data/test_gatk/bundle/1000G_phase1.indels.b37.vcf -knownSites /data/test_gatk/bundle/Mills_and_1000G_gold_standard.indels.b37.vcf  -knownSites /data/test_gatk/bundle/dbsnp_138.b37.vcf -o $3"_recalibration_report.grp" &>>$3.out

# We DON'T recalibrate the scores
#java -Xmx4g -jar /data/software/GenomeAnalysisTK.jar -T PrintReads -R /data/genomes/Broadhs37/hs37d5.fa -I $3"_realigned_cutends.srt.bam" -BQSR $3"_recalibration_report.grp" -o $3"_final.bam" &>>$3.out

