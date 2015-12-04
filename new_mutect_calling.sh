#!/bin/bash
<<"COMMENT"
The code should do the following :
Given a properly pre-processed file :
- Run mutect on the input file
- select only the variants in amplicon regions
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
  echo 'Usage : ./mutect_calling.sh <bam_file.bam> <library_number>'
  exit
}

if [ "$#" -ne 2 ]
then
  usage
fi

echo "Running MuTect pipeline on" $1 " "to create Mutect_$2".vcf.gz"

# Choose the proper amplicons file
bedfile=/data/TruSight_analysis/trusight-myeloid-amplicon-track.bed

cat $bedfile |awk '{print $1":"$2"-"$3;}' >$2".intervals"

java -Xmx4g -jar -Djava.io.tmpdir=tmp /data/software/mutect-src/mutect/target/mutect-1.1.7.jar --analysis_type MuTect --reference_sequence  /data/genomes/Broadhs37/hs37d5.fa --cosmic /data/test_mutect/b37_cosmic_v54_120711.vcf --dbsnp /data/test_mutect/dbsnp_132_b37.leftAligned.vcf.gz --input_file:tumor $1 --out $2"_call_stats.out" --coverage_file coverage.$2".wig.txt" -vcf Mutect_$2".vcf" -rf BadCigar -dt NONE --intervals $2".intervals"

bgzip Mutect_$2".vcf"
tabix Mutect_$2".vcf.gz"
