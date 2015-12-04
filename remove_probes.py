#!/usr/bin/python
'''
The code should do the following :
	- read in a bam alignment file
	- read in a list of probe regions
	- check if the read begins at the probe, if reverse if it ends at the probe
	- calculate the length of bases to remove, and then clip the read..
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
'''

import sys
import pysam
import bbcflib
from bbcflib.track import track
import clip_functions
from clip_functions import *
import argparse
import pybedtools
from pybedtools import BedTool

parser = argparse.ArgumentParser(
    description='Clip ends of reads according to probes structure')
parser.add_argument('Probes file', metavar='probes_file.bed', type=str, nargs=1,
                   help='A bed file containing the coordinates of the probes used in the design')  #In this file the bp lengths to remove are in the 11th column, from left to right
parser.add_argument('Alignments bam', metavar='in_file.bam', type=str, nargs=1,
                   help='A bam file containing the alignments to process')
parser.add_argument('Output file', metavar='out_file.bam', type=str, nargs=1,
                   help='A bam file containing the output')

args=parser.parse_args()
myargs=vars(args)

samfile = pysam.AlignmentFile(myargs["Alignments bam"][0], "rb")
pairedreads = pysam.AlignmentFile(myargs["Output file"][0], "wb", template=samfile)

amplicons=BedTool(myargs["Probes file"][0])

maxdist=6

seqnames=list()

for a in amplicons:
        seqnames.append(a.chrom)

def chrom_filter(feature,chrom):
	return feature.chrom==chrom

for chr in set(seqnames):
	left_lengths=dict()
	right_lengths=dict()
	a_starts=dict()
	a_ends=dict()
	amplicons_chrom=amplicons.filter(chrom_filter, chrom=chr)
	for a in amplicons_chrom:
		left_lengths[a.name]=map(int,a.fields[10].split(","))[0]
		right_lengths[a.name]=map(int,a.fields[10].split(","))[1]
		a_starts[a.name]=a.start
		a_ends[a.name]=a.stop

	chrom=chr[3:]
	print chrom
	for read in samfile.fetch(str(chrom)):
		if read.is_reverse:
			if read.is_unmapped==False:
				dists=dict()
				abs_dists=dict()
				for k, v in a_ends.items():
					dists[k]=read.reference_end-v
					abs_dists[k]=abs(dists[k])	#This should return which amplicon the read belongs to. Tolerance for errors should be written in, and the dist returned
				if min(abs_dists.values())<maxdist:
					ampliconID=min(abs_dists,key=abs_dists.get)
					clip_right(read,right_lengths[ampliconID]+dists[ampliconID])
		else:
			if read.is_unmapped==False:
				dists=dict()
				abs_dists=dict()
				for k, v in a_starts.items():
					dists[k]=read.reference_start-v
                                        abs_dists[k]=abs(dists[k])      #This should return which amplicon the read belongs to. Tolerance for errors should be written in, and the dist returned
                                if min(abs_dists.values())<maxdist:
                                        ampliconID=min(abs_dists,key=abs_dists.get)
                                        clip_left(read,left_lengths[ampliconID]-dists[ampliconID])	#Same but based on the start
		#print read unless it's now empty of any alignment then we have to set it to unmapped :
		if read.is_unmapped==False:
			if len(read.cigar)==0 or 0 not in (item[0] for item in read.cigar):
				read.is_unmapped=True
				read.cigar=""
				read.mapq=0 
				read.is_secondary=False
		if read.rlen>0:
			pairedreads.write(read)

#Let's add something to retain the rest of the reads :
samfile.close()
samfile = pysam.AlignmentFile(myargs["Alignments bam"][0], "rb")

for read in samfile.fetch(until_eof=True):
	if read.tid==-1:
		pairedreads.write(read)
	elif "chr"+samfile.getrname(read.tid) not in set(seqnames):
		pairedreads.write(read)


pairedreads.close()
samfile.close()

