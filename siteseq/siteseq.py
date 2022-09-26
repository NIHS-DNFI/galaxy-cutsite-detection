# -*- coding: utf-8 -*-

"""
This python script was modified from a published research as below:
Cameron, P., Fuller, C., Donohoue, P. et al.
Mapping the genomic landscape of CRISPR–Cas9 cleavage.
Nat Methods 14, 600–606 (2017). https://doi.org/10.1038/nmeth.4284
"""

from __future__ import division
import os
import sys
import pysam
import re
from pyfaidx import Fasta
from collections import Counter
import argparse

### MAIN FUNCTIONS ###
# find_initial_read_pileups_modc #
"""
Identifies any pileups in aligned sequencing reads that exceed
a user specified threshold. Returns the genomic locus (coordinates)
of contiguous sequence above the pileup threshold.

Args:
    bam_file (str):
        The full path to a bam file that has been sorted and
        indexed.

    outputfile (str):
        The full path of output file name.

    depth_tresh (int):
        The minimum number of reads required to call a pileup.
        Default: 5
    banned_chroms (list):
        A list of strings that cannot exist within the chromosome
        name for a given pileup location.
        Default: ["alt", "random", "Un", "chrM"]
Returns:
    peak_locations (list):
        A list of genomic coordinates in the form chr1:1234-5678
        that exceed the read pileup threshold.
"""


def find_initial_read_pileups_modc(bam_file,
                                outputfile_one,
                                banned_chroms=["alt", "random", "Un"],
                                depth_thresh=5):
                f = open(outputfile_one, 'w')

                global peak_locations
                pysam.index(bam_file)
                bam = pysam.AlignmentFile(bam_file, 'rb')
                pileup = bam.pileup()
                in_peak = False
                # Various parameters for the peak identification.
                peak_start = -1
                peak_end = -1
                last_chrom = -1
                peak_start = -1
                peak_end = -1
                last_pos = -1
                for col in pileup:
                    # Ignore undesired chromosomes.
                    bad_chroms = [b for b in banned_chroms if b in col.reference_name]
                    if any(bad_chroms):
                          continue
                    if in_peak:
                        # If any of these are true, terminate the peak.
                        if(int(col.nsegments) < int(depth_thresh)) or \
                                (int(col.reference_pos) != (last_pos + 1)) or \
                                (col.reference_name != last_chrom):
                            in_peak = False
                            peak_end = last_pos
                            # Ignore 1 length peaks.
                            if peak_end == peak_start:
                                continue
                            coord = "{c}:{s}-{e}".format(c=last_chrom,
                                                        s=peak_start,
                                                        e=peak_end)
                            peak_locations.append(coord)
                    else:
                        if(int(col.nsegments) >= int(depth_thresh)):
                            in_peak = True
                            peak_start = col.reference_pos
                    last_pos = col.reference_pos
                    last_chrom = col.reference_name
                f.write("\n".join(peak_locations))
                f.close

# call_site_seq_features #
"""
    This function takes peaks and determines whether or
    not they are site-seq features. The steps for feature
    detection are as follows:

    1) Remove peaks that lie within undesired chromosomes.
    2) For each peak, identify the count for locations where sequencing
       reads terminate (cliffs, either 5' or 3') at the same position.
       Reads that are not long enough (determined by min_read_length)
       are removed from the counts.
    3) Ensure that the number of cliff reads for that peak is high enough
       (determined by min_read_length). If the read count for the cliff is
       too low the potential feature is ignored.
    4) Ensure that the ratio of non-cliff reads to cliff reads is not
       too high. If the ratio is too high (determined by required_cliff_ratio)
       then the potential feature is ignored.
    5) All potential features that pass the above filters are reported.

    Args:
        initial_peaks (list):
            A list of peaks already found. These are strings formatted like
            'chr#:#-#'
        bam_file (str):
            The path and name of a sorted and indexed bam file.
        pyfaidx_fasta_obj (pyFaidx Fasta Object):
            A pyFaidx object created from a reference genome.
        out_name (str):
            The name/path of the output fasta file. The fasta file contains
            information about the peak in the header as well as the sequence
            around
        min_dup (int):
            The number of reads that terminate at the exact same location
            (cliff_reads) needed for a potential feature to be considered.
            Default: 5 reads
        min_read_length (int):
            The minimum read length required. Reads shorter than this will
            be removed.
            Default: 60.
        flank_from_feature (int):
            Sets the sequence around the SITE-Seq feature to return in the
            fasta file output by the function.
            Default: 20.
        required_cliff_ratio (float):
            The minimum ratio of cliff reads to cliff spanning reads
            required in order for a feature to be considered.
            Default: 0.9.
        banned_chroms (list):
            A list of chromosome names that should not be considered when
            attempting to identify features. This may need to be
            changed depending on the reference used.

    Returns:
        out_name (str):
            The path to the fasta file that contains the features.
"""

def call_site_seq_features(initial_peaks, bam_file,
                              pyfaidx_fasta_obj, out_name,
                              min_dup=5, min_read_length=60,
                              flank_from_feature=20,
                              required_cliff_ratio=0.9,
                              banned_chroms=["alt", "random", "Un", "chrM"]):
            
            bam = pysam.AlignmentFile(bam_file, 'rb')
            features = []

            for i, peak in enumerate(initial_peaks):
                r1_starts = []
                chrom, start, end = parse_coordinate(peak)
                if any(p for p in banned_chroms if p in chrom):
                    continue
                try:
                    reads = bam.fetch(chrom, start - 5, end + 5)
                except:
                    continue
                for read in reads:
                    read_start = read.reference_start
                    read_end = read.reference_end
                    read_on_minus = read.is_reverse
                    if read_start is None or read_end is None:
                        continue
                    aligned_length = abs(int(read_start) - int(read_end))
                    if aligned_length < min_read_length:
                        continue
                    read_append = read_end if read_on_minus else read_start
                    r1_starts.append(read_append)
                terminations = Counter(r1_starts)
                if len(terminations) == 0:
                    continue
                r1_start = terminations.most_common(1)[0][0]
                r1_count = terminations.most_common(1)[0][1]
                if r1_count < min_dup:
                    continue
                non_cliff_reads_at_cliff = 0
                if not r1_start:
                    continue
                for read in bam.fetch(chrom, r1_start - 1, r1_start + 1):
                    if read.reference_start != r1_start and \
                            read.reference_end != r1_start:
                        non_cliff_reads_at_cliff += 1
                try:
                    cliff_ratio = non_cliff_reads_at_cliff / r1_count
                except ZeroDivisionError:
                    cliff_ratio = 0
                if cliff_ratio >= required_cliff_ratio:
                    continue
                fasta_dict = {}
                fasta_dict["feature_region"] = "{c}:{s}-{e}".format(
                    c=chrom, s=r1_start - flank_from_feature,
                    e=r1_start + flank_from_feature)
                fasta_dict["feature_cut"] = "{c}:{s}".format(c=chrom,
                                                              s=r1_start)
                fasta_dict["r1_start_count"] = r1_count
                fasta_dict["reads_near_r1"] = non_cliff_reads_at_cliff
                features.append(fasta_dict)
            fasta_from_features(out_name, features, pyfaidx_fasta_obj)

## HELPER FUNCTIONS refered in MAINs ##
def fasta_from_features(out_name, features, pyfaidx_fasta_obj):
    """
    This function requires a reference genome in a fasta file.
    See https://github.com/mdshw5/pyfaidx for more information.
    """
    with open(out_name, "w+") as f:
        fileheader = "region\tcut_site\tnoncliff_reads_at_cliff\tr1_start_count\tsequence\n"
        f.write(fileheader)
        for feature in features:
            sequence = retrieve_sequence(feature["feature_region"],
                                         pyfaidx_fasta_obj)
            header = ">{region}\t{cut}\t{r1_start}\t{start_count}".format(
                region=feature["feature_region"], cut=feature["feature_cut"],
                r1_start=feature["reads_near_r1"],
                start_count=feature["r1_start_count"])
            f.write(header + "\t")
            f.write(sequence + "\n")
    return out_name

def retrieve_sequence(coordinate, pyfaidx_fasta_obj):
    """
    Obtain sequence from a pyfaidx fasta object given a
    coordinate.
    """
    chrom, start, stop = parse_coordinate(coordinate)
    return pyfaidx_fasta_obj[chrom][start - 1:stop].seq


def parse_coordinate(coordinate):
    """
    Parse a coordinate formatted like :
    chrZ:#-# OR chrZ:# where Z and # are placeholders.
    """
    coord_range_pattern = r'(.*):\d+-\d+$'
    coord_site_pattern = r'(.*):\d+$'
    if re.match(coord_range_pattern, coordinate):
        chrom = coordinate.split(':')[0]
        start = int(coordinate.split(':')[1].split('-')[0])
        stop = int(coordinate.split(':')[1].split('-')[1])
        return chrom, start, stop
    elif re.match(coord_site_pattern, coordinate):
        chrom = coordinate.split(':')[0]
        start = int(coordinate.split(':')[1])
        return chrom, start

### Main script codes ###
"""
 1st function:
    Extracts rougly all of the peaks above the read number threshold and returns the peak location list as a list named 'peal_locations'. In this script, the output is also saved as a tab-delimited txt file in a designated file path.
 2nd function:
    Filters the result list from the 1st def with the SITEseq-peak features and returns the final SITE-seq target candidate list in the designated path.
"""
## argument for 1st def defined in the global scope
peak_locations = []

def main():
    parser = argparse.ArgumentParser(
	    description='This is a program to call the SITE-seq feature and suggest the possible off-target site of CRISPR reaction.',
	    prog='SITEseq_peak_call_test.py',
	    usage='%(prog)s -i [input file path] -R [reference fasta file path] -o [output file path] [other options]',
	    epilog='Make sure you include .txt extension in your input/output file path.',
	    add_help=True,
	    )

#===== parser arguments for 1st def =====#
#input bam file
    parser.add_argument('-i', '--input', action="store", dest="inputbam", required=True, help="Path to the input bam file; This is REQUIRED argument") #= siteseq_bam
#output file for 1st def
    parser.add_argument('-p', '--output1', dest="allpileup", action="store", help="Path to the output to the 1st function which lists all pileup sites") #= outputfile_allpileup
#ng_chrom = ["alt", "random", "Un"]
#minimum read depth for pileup call
    parser.add_argument('-t', '--minread', action="store", dest="minread", type=int, default=5, help="Read depth threshold for initial peak call") #= read_thresh
#===== parser arguments for 2nd def =====#
#reference genome fasta file
    parser.add_argument('-R', '--reference', action="store", dest="reference", required=True, help="Path to the reference genome fasta file; This is REQUIRED argument") #= reference_genome
#final output file of the sites with SITEseq feature (output of 2nd def)
    parser.add_argument('-o', '--output2', action="store", dest="SITEseqpeak", required=True, help="Path to the final output file which lists SITEseq cut locations; This is REQUIRED argument") #= outputfile_SITEseqpeak
#read depth threshold for SITEseq feature
    parser.add_argument('-m', '--mindup', action="store", dest="mindup", type=int, default=5, help="Minimum read depth which ends at the cliff required to be considered as the SITEseq feature") #= min_dup
#minimum read length for SITEseq feature
    parser.add_argument('-l', '--minlen', action="store", dest="minlen", type=int, default=60, help="Minimum read length required to be considered as the SITEseq feature") #= min_read_length
#sequence length to be returned around the SITEseq cut position
    parser.add_argument('-f', '--flank', action="store", dest="flank", type=int, default=20, help="The sequence length around the detected SITEseq feature to be returned in the output file") #= flank_from_feature
#minimum ratio of cliff reads for SITEseq feature
    parser.add_argument('-r', '--cliffratio', action="store", dest="cliffratio", type=float, default=0.9, help="The minimum ratio of cliff reads to cliff spanning reads required in order for a feature to be considered") #= required_cliff_ratio

    args = parser.parse_args()
    siteseq_bam = args.inputbam
    outputfile_allpileup = args.allpileup
    ng_chrom = ["alt", "random", "Un"]
    read_thresh = args.minread

    reference_genome = args.reference
    outputfile_SITEseqpeak = args.SITEseqpeak
    pyfaidx_fasta_obj = Fasta(reference_genome)
    min_dup = args.mindup
    min_read_length = args.minlen
    flank_from_feature = args.flank
    required_cliff_ratio = args.cliffratio

    print('SITEseq bam file is', siteseq_bam)
    print('Reference genome file is', reference_genome)
    print('Read depth threshold for peak call is', read_thresh)
    print('Minimum read depth which ends at the cliff required to be considered as the SITEseq feature is', min_dup)
    print('Minimum read length required to be considered as the SITEseq feature is', min_read_length)
    print('The sequence length around the detected SITEseq feature to be returned in the output file is',  flank_from_feature)
    print('The minimum ratio of cliff reads to cliff spanning reads required in order for a feature to be considered is',  required_cliff_ratio)
    print('First output file will be', outputfile_allpileup)
    print('Second output file will be', outputfile_SITEseqpeak)
    print("")

    find_initial_read_pileups_modc(siteseq_bam, outputfile_allpileup, ng_chrom, read_thresh)
    call_site_seq_features(peak_locations, siteseq_bam, pyfaidx_fasta_obj, outputfile_SITEseqpeak)

if __name__ == '__main__':
    main()
