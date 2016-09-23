#!/usr/bin/python
#######################################################################################
# Analysis of 2D signal in mutational profiling sequencing data (MaP2D)
#######################################################################################
#
# Simple to RDAT
#
# This script generates a 2D correlated mutational profiling dataset using data in 'simple' binary format that records the mutations in sequenced cDNAs
#
# The outputs of the script are:
#        1. An RDAT file containing the MaP2D (two-dimensional mutational profiling) data
#        2. A .log file recording minor output from the analysis, e.g. input files
#
# Clarence Cheng, 2015-2016
# with advice and code from Pablo Cordero and Joe Yesselman
#

import os, sys, time
import argparse
import numpy as np
import nwalign as nw
import sys
import string
import re
from itertools import izip_longest
from rdatkit.datahandlers import RDATFile
from matplotlib.pylab import *

parser = argparse.ArgumentParser( description='Generates a 2D correlated mutational profiling dataset using demultiplexed Read1 and Read 2 files, aligned Read1 and Read2 text files, or a .simple file')

parser.add_argument('sequencefile', type=argparse.FileType('r'), nargs='+', help='name of sequence file in .fasta format')
parser.add_argument('--simplefile', type=argparse.FileType('r'), help='name of .simple file')
parser.add_argument('--name',       type=str,  help='name of RNA; by default take from sequencefile', default='PLACEHOLDER')
parser.add_argument('--offset',     type=int,  help='integer to add to residue numbers to get conventional numbers', default=0)
parser.add_argument('--outprefix',  type=str,  default='out')
parser.add_argument('--num_hits_cutoff',        type=int,  help='quality filter: maximum number of hits to allow before recording', default=10)
parser.add_argument('--num_hits_cutoff_lower',  type=int,  help='filter minimum number of hits, for testing signal', default=-1)
parser.add_argument('--start_pos_cutoff',       type=int,  help='full-length filter: minimal nucleotide position to which alignment must extend', default=10)
parser.add_argument('--end_pos_cutoff',         type=int,  help='full-length filter: minimal number of nucleotides from 3''-end at which alignment must start', default=10)        ############################################################################################################
parser.add_argument('--simplelinelengths',      type=bool,  help='Output file with both length of read and num mutations', default=False)

args = parser.parse_args()

i = 0
seqs = []
names = []
for seqfile in args.sequencefile:
    seqfile_lines = seqfile.readlines()
    seqs.append( seqfile_lines[1].strip().upper().replace("U","T") )
    names.append( seqfile_lines[0].strip()[2:] )
    print names[i]
    print seqs[i]
    i += 1

# if args.name == 'PLACEHOLDER':
#     args.name = sequencefile_lines[0].strip()[1:]


def timeStamp():        # From ShapeMapper
    t = time.localtime()
    month = t.tm_mon
    day = t.tm_mday
    hour = t.tm_hour
    minute = t.tm_min
    second = t.tm_sec
    return '%i:%i:%i, %i/%i'%(hour, minute, second, month, day)


def output_rdat( filename, args, sequence, name, row_indices, seqpos, WTdata, WTdata_err, data2d, data2d_err, f_log ):
        # Output data and annotations to RDAT file
        construct = name
        sequence_RNA = string.replace(sequence,'T','U')
        structure = ('.'*len(sequence))
        offset = args.offset
        version = 0.34
        comments = 'Created from %s using simple_to_rdat.py\n' % (args.simplefile.name)
                
        annotations = {'experimentType:MutationalProfiling'}
        data_annotations = []
        data_annotations.append({'mutation':'WT'})
        if len(row_indices)> 0:
            for idx, row_idx in enumerate(row_indices):
                position = seqpos[row_idx] - offset
                data_annotations.append({'mutation':[string.replace(sequence[position-1],'T','U') + str(position + offset) + 'X']})
        data = np.concatenate([WTdata, data2d], axis=0)                        # combine WT data and 2D data

        r = RDATFile()
        r.save_a_construct(construct, data, sequence_RNA, structure, offset, annotations, data_annotations, filename, comments, version)
        # Note: Should include # of simple file lines, filters used, and # mutants passing filters in RDAT file
        #       and edit read_rdat_file.m in HiTRACE to read these figures; will streamline plotting

        output_string = 'RDAT file created: ' + filename
        sys.stdout.write( output_string+'\n' )
        f_log.write( output_string+'\n' )


def simple_to_rdat( args, sequence, name ):
        print 'Starting analysis from simple file at ' + timeStamp()
        f_log.write( 'Starting analysis from simple file at ' + timeStamp() + '\n\n' )
        
        #### Get length and seqpos
        WTlen = len(sequence)
        seqpos = arange(1, WTlen+1) + args.offset
        
        #### Get mutation indices from simple data
        print 'Reading file: '+str(args.simplefile.name)

        # #### Get data for distribution of start and end positions of reads, length and # muts per read
        # f_startendpos = open(currdir + '/' + args.outprefix + '_startendpos.csv', 'w')
        # if args.simplelinelengths:
        #     f_lengthandmuts = open(currdir + '/' + args.outprefix + '_readlenandmuts.csv', 'w')
        #     f_lengthandmuts.write( 'Read_length,Number_of_muts\n' )

        #### Set up arrays for 1D and 2D data
        count = 0
        fltnm = 0
        usenm = 0
        stpnm = 0
        WTdata = np.zeros((1, WTlen))
        data2d = np.zeros(([len(seqpos),len(seqpos)]))
        mut_idxs   = []
        mut_counts = []
        mut_mat    = np.zeros((1,WTlen+1))
        start_pos_counts = np.zeros((1,WTlen))
        start_pos_for_norm = np.zeros(WTlen)
        stop_reactivity  = np.zeros((1,WTlen))

        #### Read simple file and generate 2D data
        for line in args.simplefile:
            # Read from file
            fields = line.strip().split('\t')
            if len( fields ) < 3: continue
            start_pos = int(fields[0])
            end_pos = int(fields[1])
            # f_startendpos.write(str(start_pos) + ',' + str(end_pos) + '\n')
            simple_string = fields[2].strip()   # Remove endline
            assert( len(simple_string) == (int( fields[1]) - int(fields[0]) + 1 ) )
            simple_line = [ int(x) for x in simple_string ]
            count += 1                          # Record total sequences
            if count % 50000 == 0: print 'Reading line number ', count

            # Simple-to-RDAT stop and mutation assignment
            mut_counts = 0
            start_pos_to_use = start_pos
            # Record mutation indices
            mut_idx_init = array( [ (idx + start_pos_to_use - 1) for idx, val in enumerate(simple_line) if val == 1])
            mut_idx = [idx for idx in mut_idx_init if idx < WTlen]
            mut_counts = len( mut_idx )
            # if args.simplelinelengths:
            #     f_lengthandmuts.write( str(len(simple_string)) + ',' + str(mut_counts) + '\n' )
            if ( start_pos <= args.start_pos_cutoff ): mut_mat[0,mut_counts] += 1       # record mutations per read
            if ( mut_counts > args.num_hits_cutoff ): continue
            if ( mut_counts <= args.num_hits_cutoff_lower ): continue
            # Build data arrays
            if ( start_pos <= args.start_pos_cutoff ):
                if ( end_pos >= WTlen - args.end_pos_cutoff ):
                    if mut_counts > 0:
                        WTdata[ 0, mut_idx ] += 1               # Build 1D profile
                        usenm += 1
                    for mutpos in mut_idx:
                        data2d[mutpos, mut_idx] += 1                # Build 2D dataset
                    fltnm += 1                                      # count sequences that pass filters 
            # Record reactivity associated with reverse transcriptase stopping (probably should look at a 2D version of this, like in MOHCA-seq. -- rhiju)
            assert ( start_pos_to_use  <= (WTlen + 1) )
            if ( start_pos_to_use > WTlen ): continue
            mod_pos = start_pos_to_use - 1                      # "start pos" actually records where reverse transcriptase stopped.
            # note further offset: mod_pos=0 means *no modifications*; mod_pos = 1 means modification at position 1
            # if we swtich to to 2D readout, maybe should decrement (add -1). Then no mod would actually *wrap* to WTlen (?)
            start_pos_counts[ 0, mod_pos ] += 1
            stpnm += 1                                          # count sequences used for stop reactivities

        # if args.simplelinelengths:
        #     f_lengthandmuts.close()

        f_dist = open(currdir + '/' + args.outprefix + '_mutsperread.csv','w')
        f_dist.write( 'Number of mutations,Number of reads\n')
        for (x,y), nummuts in np.ndenumerate(mut_mat): f_dist.write( str(y) + ',' + str(nummuts) + '\n')
        f_dist.close()

        f_log.write( '\nTotal number of sequences: ' + str(count) )
        f_log.write( '\nNumber of sequences used for RT stop reactivities: ' + str(stpnm) )
        f_log.write( '\nFull-length filter: alignment must extend to at least position ' + str(args.start_pos_cutoff) )
        f_log.write( '\nQuality filter (use for 2D data if N or fewer mismatches to WT): N = ' + str(args.num_hits_cutoff) )
        f_log.write( '\nNumber of sequences passing full-length and quality filters: ' + str(fltnm) )
        f_log.write( '\nNumber of sequences used for 2D data (at least 1 mutation): ' + str(usenm) + '\n\n')

        # #### Get stop reactivities, using Fi/[F0 + F1 ... + Fi] expression.
        # sum_counts = start_pos_counts[ 0, 0 ]
        # for (idx,counts) in enumerate(start_pos_counts[0,1:]):
        #     sum_counts += counts
        #     stop_reactivity[0, idx ] = ( float(counts)/sum_counts )
        # filename = currdir + '/' + args.outprefix + '_' + name + '.stop_reactivity.rdat'
        # output_rdat( filename, args, sequence, name, [], seqpos, stop_reactivity, np.zeros((0,WTlen)), f_log )

        if fltnm == 0:
            f_log.write( '\nNo sequences passing both full-length and quality filters! Setting normalization factor for "modification fraction" to 1\n' )
            fltnm = 1
        
        #### Get reactivity for WT data
        # normalizes to total number of reads -- should be a true 'modification fraction'
        WTdata_reactivity       = WTdata/fltnm
        WTdata_err              = np.sqrt(WTdata)
        WTdata_reactivity_err   = WTdata_err/fltnm
        # # normalizes to total number of hits -- original choice
        # WTdata_norm = WTdata/WTdata.sum()

        #### Get reactivity for 2D data
        # Normalize by total signal in row and store row indices for building RDAT file
        row_indices = []
        data2d_reactivity       = 0 * data2d
        data2d_err              = np.sqrt(data2d)    # estimate error for each 2D position as square root of number of reads with 2D signal at that position
        data2d_reactivity_err   = 0 * data2d_err
        # data2d_norm             = 0 * data2d
        for row_idx in xrange(data2d.shape[0]):
            if data2d[row_idx,:].sum() > 0*data2d.shape[1]:
                # normalizes to total number of reads with modification at row_idx position -- should be a true 'modification fraction'
                data2d_reactivity[row_idx, :]       = data2d[row_idx,:]     / data2d[ row_idx, row_idx ]      # [row_idx, row_idx] contains total reads with mod at row_idx position
                data2d_reactivity_err[row_idx, :]   = data2d_err[row_idx,:] / data2d[ row_idx, row_idx ]      # normalize errors by same scale factor as reactivities
                # # normalizes to total number of hits -- original choice
                # data2d_norm[row_idx, :] = data2d[ row_idx, : ] / data2d[row_idx, :].sum()
            row_indices.append(row_idx)

        #### Output RDATs
        # reactivity RDAT (normalize by total reads with a mutation at row_idx)
        filename = currdir + '/' + args.outprefix + '_' + name +  '.reactivity.rdat'
        output_rdat( filename, args, sequence, name, row_indices, seqpos, WTdata_reactivity, WTdata_reactivity_err, data2d_reactivity, data2d_reactivity_err, f_log )

        # raw RDAT (# reads not normalized; can subsequently normalize by e.g. total reads per barcode, total aligned reads, total reads with >=1 mutation)
        filename = currdir + '/' + args.outprefix + '_' + name +  '.raw.rdat'
        output_rdat( filename, args, sequence, name, row_indices, seqpos, WTdata, WTdata_err, data2d, data2d_err, f_log )

        # filename = currdir + '/' + args.outprefix + '_' + name +  '.rdat'
        # output_rdat( filename, args, sequence, name, row_indices, seqpos, WTdata_norm, data2d_norm, f_log )

        f_log.write( '\nFinished analysis from simple file at ' + timeStamp() + '\n\n' )

        return


if __name__ == '__main__':

        # get target directory ready.
        targetdir = os.path.dirname( args.outprefix )
        if ( len( targetdir ) > 0 and not os.path.exists( targetdir ) ): os.makedirs( targetdir )

        # Get current directory and open log file
        currdir = os.getcwd()
        f_log = open(currdir + '/' + args.outprefix + '.log', 'w')
        f_log.write( 'Current directory: ' + currdir + '\n\n' )

        # See if user input .simple data file and run commands
        if args.simplefile is not None:
            f_log.write( 'Simple data file provided:\t' + args.simplefile.name + '\nGenerating 2D data from simple data file.\n\n' )
            for i in xrange(len(names)):
                simple_to_rdat( args, seqs[i], names[i] )

        else:
            f_log.write( 'No data files provided.\n' )
            exit()

        f_log.close()



