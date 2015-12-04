#!/bin/bash
<<"COMMENT"
The code should do the following :
Given a properly pre-processed file :
- Run GATK HaplotypeCaller on the bam file in discovery mode
- DON'T Prepare the Recalibration files with VariantRecalibrator and run recalibration on SNPs using databases of known constitutive SNPS
- DON'T Idem for INDELS
Future improvements planned :
- Include cancer-specific mutations
- Perform initial calling with stringent parameters then run VQSR using this set of high confidence mutations
  to avoid underestimating the mutational load in cancer samples
- Possibly iterative INDEL calling
Disclaimer : GATK isn't typically good to call somatic variants because it relies on an assumption of heterozygozity or homozygosity, i.e. AF or 0, 0.5 or 1
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
  echo 'Usage : ./gatk_calling.sh <bam_file.bam> <library_number>'
  exit
}

if [ "$#" -ne 2 ]
then
  usage
fi

echo "Running gatk pipekine on" $1 " "to create GATK_$2".vcf.gz"

# Choose the proper amplicons file
bedfile=/data/TruSight_analysis/trusight-myeloid-amplicon-track.bed

java -jar -Djava.io.tmpdir=tmp /data/software/GenomeAnalysisTK.jar -T HaplotypeCaller -R /data/genomes/Broadhs37/hs37d5.fa -I $1 -L $bedfile  --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o GATK_$2".vcf" -rf BadCigar -dt NONE

#java -jar -Djava.io.tmpdir=tmp /data/software/GenomeAnalysisTK.jar -T VariantRecalibrator    -R /data/genomes/Broadhs37/hs37d5.fa     -input $2".vcf"   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /data/test_gatk/bundle/hapmap_3.3.b37.vcf     -resource:omni,known=false,training=true,truth=true,prior=12.0 /data/test_gatk/bundle/1000G_omni2.5.b37.vcf     -resource:1000G,known=false,training=true,truth=false,prior=10.0 /data/test_gatk/bundle/1000G_phase1.snps.high_confidence.b37.vcf     -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /data/test_gatk/bundle/dbsnp_138.b37.vcf     -an DP     -an QD     -an FS     -an SOR     -an MQ    -an MQRankSum     -an ReadPosRankSum    -mode SNP     -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0     -recalFile recalibrate_$2".recal"     -tranchesFile recalibrate_$2".tranches"     -rscriptFile recalibrate_$2".R" -dt NONE 

#java -jar -Djava.io.tmpdir=tmp /data/software/GenomeAnalysisTK.jar -T ApplyRecalibration     -R /data/genomes/Broadhs37/hs37d5.fa -input $2".vcf"      -mode SNP     --ts_filter_level 99.0     -recalFile /data/test_gatk/testing_bams/recalibrate_$2".recal"    -tranchesFile /data/test_gatk/testing_bams/recalibrate_$2".tranches"     -o recalibrated_snps_raw_indels_$2".vcf" 

#We cannot recalibrate variants because we have too few, since we have a small targeted experiment
#java -jar -Djava.io.tmpdir=tmp /data/software/GenomeAnalysisTK.jar -T VariantRecalibrator     -R /data/genomes/Broadhs37/hs37d5.fa    -input /data/test_gatk/testing_bams/recalibrated_snps_raw_indels_$2".vcf"     -resource:mills,known=true,training=true,truth=true,prior=12.0 /data/test_gatk/bundle/Mills_and_1000G_gold_standard.indels.b37.vcf    -an QD     -an DP     -an FS     -an SOR     -an MQRankSum     -an ReadPosRankSum     -mode INDEL     -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0     --maxGaussians 4     -recalFile recalibrate_INDEL_$2".recal"     -tranchesFile recalibrate_INDEL_$2".tranches"     -rscriptFile recalibrate_INDEL_plots_$2".R" 

#java -jar -Djava.io.tmpdir=tmp /data/software/GenomeAnalysisTK.jar -T ApplyRecalibration     -R /data/genomes/Broadhs37/hs37d5.fa     -input /data/test_gatk/testing_bams/recalibrated_snps_raw_indels_$2".vcf"     -mode INDEL     --ts_filter_level 99.0     -recalFile recalibrate_INDEL_$2".recal"     -tranchesFile recalibrate_INDEL_$2".tranches"     -o GATK_$2".vcf" 


bgzip GATK_$2".vcf"
tabix GATK_$2".vcf.gz"

