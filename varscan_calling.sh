#!/bin/bash
<<"COMMENT"
The code should do the following :
Given a properly pre-processed file :
- Run samtools mpileup on the bam file
- Run varscan2 on the resulting mpileup
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
  echo 'Usage : ./varscan_calling.sh <bam_file.bam> <library_number>'
  exit
}

if [ "$#" -ne 2 ]
then
  usage
fi
      
echo "Running varscan pipeline on" $1 " "to create Varscan_$2".vcf.gz"

# Choose the proper amplicons file
bedfile=/data/TruSight_analysis/trusight-myeloid-amplicon-track.bed

samtools mpileup -f /data/genomes/Broadhs37/hs37d5.fa  -l $bedfile $1 >$2".mpileup"

java  -Xmx4g -jar -Djava.io.tmpdir=tmp /data/software/VarScan.v2.3.9.jar mpileup2snp $2".mpileup" --min-coverage 100 --min-var-freq 0.05 --min-avg-qual 20 --p-value 0.1 --output-vcf >Varscan_$2".vcf"  --strand-filter 0 --min-strand 0

bgzip Varscan_$2".vcf"
tabix Varscan_$2".vcf.gz"

