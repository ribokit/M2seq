#!/usr/bin/env python
#######################################################################################
# Analysis of 2D signal in mutational profiling sequencing data (MaP2D)
#######################################################################################
#
# FASTQ to RDAT
#
# This script generates a 2D correlated mutational profiling dataset using one of two sources:
# 	1. Demultiplexed and quality-filtered Read 1 and Read 2 FASTQ files
# 	2. Data in 'simple' binary format that records the mutations in sequenced cDNAs
# 	If a 'simple' file is input, the script does not use FASTQ files.
#
# The outputs of the script are:
# 	1. An RDAT file containing the 2D MaP data
# 	2. A 'simple' file that records the mutations in each cDNA in binary format
# 	3. A text file for each input FASTQ showing all sequences after alignment ('...Read1Aligned.txt' and '...Read2Aligned.txt')
#	4. A '...mutations.csv' file quantifying the numbers of mutations and mutation rates
#	5. A .log file recording minor output from the analysis, e.g. input files
#	If a 'simple' file is input, only the RDAT and log files are generated.
#
# Clarence Cheng, 2015
# with advice and code from Pablo Cordero and Joe Yesselman
#

import os, sys, time
import argparse
import numpy as np
import nwalign as nw
import string
from rdatkit.datahandlers import RDATFile
from matplotlib.pylab import *

parser = argparse.ArgumentParser()

parser.add_argument('sequencefile', type=argparse.FileType('r'))
parser.add_argument('--read1fastq', type=argparse.FileType('r'))
parser.add_argument('--read2fastq', type=argparse.FileType('r'))
parser.add_argument('--simplefile', type=argparse.FileType('r'))
parser.add_argument('--name', type=str, default='PLACEHOLDER')
parser.add_argument('--offset', type=int, default=0)
parser.add_argument('--outprefix', type=str, default='out')

args = parser.parse_args()
sequencefile_lines = args.sequencefile.readlines()
sequence = sequencefile_lines[1].strip().upper()

if args.name == 'PLACEHOLDER':
	args.name = sequencefile_lines[0].strip()[1:]


def timeStamp():										# From ShapeMapper
    t = time.localtime()
    month = t.tm_mon
    day = t.tm_mday
    hour = t.tm_hour
    minute = t.tm_min
    second = t.tm_sec
    return '%i:%i:%i, %i/%i'%(hour, minute, second, month, day)


def reverse_complement( sequence ):
	sequence_dict = {'A':'T','T':'A','C':'G','G':'C'}
	return ''.join([sequence_dict[nt] for nt in reversed(sequence)])


def fastq_to_simple( args ):			###### Build 'simple' array with 1 at each mutation and 0 at each WT nt

	# Get WT sequence, rev comp, and length
	f_log.write( 'Starting analysis at ' + timeStamp() )
	f_log.write( '\n\nRNA name: ' + args.name )
	WTfwd = sequence
	WTrev = reverse_complement(WTfwd)
	WTlen = len(WTfwd)
	seqpos = arange(0, WTlen) + args.offset
	f_log.write( '\nLength of WT sequence: ' + str(WTlen) + '\n' )
	f_log.write( '\n\nLength of seqpos: ' + str(len(seqpos)) )
	f_log.write( '\nWild-type sequence for matching to Read 2:\n' + WTfwd + '\n' )
	f_log.write( '\nReverse complement for matching to Read 1:\n' + WTrev + '\n' )

	# Grab cDNA sequences from fastq files
	seqs_read1 = []
	for i, line in enumerate(args.read1fastq):
		if i % 4 == 1:
			seqs_read1.append( line.strip() )
	seqs_read2 = []
	for i, line in enumerate(args.read2fastq):
		if i % 4 == 1:
			seqs_read2.append( line.strip() )

	# Truncate WT reverse complement to match length of read
	WTrev_trunc = WTrev[0:len(seqs_read1[1])]
	f_log.write( '\nTruncated WT rev sequence: ' + WTrev_trunc )
	f_log.write( '\nFirst read 1 sequence:     ' + seqs_read1[1] + '\n' )
	WTfwd_trunc = WTfwd[0:len(seqs_read2[1])]
	f_log.write( '\nTruncated WT fwd sequence: ' + WTfwd_trunc )
	f_log.write( '\nFirst read 2 sequence:     ' + seqs_read2[1] + '\n' )

	if len(seqs_read2) != len(seqs_read1):
		print 'Number of sequences in read 1 and read 2 FASTQs are unequal!'
		exit()

	# Align WT sequence to read1 and read2
	f_log.write( '\nStarting alignment at ' + timeStamp() )
	seqs_read1_align = []
	seqs_read2_align = []
	f_seqs_read1 = open(currdir + '/' + args.outprefix + '_Read1aligned.txt','w')
	f_seqs_read2 = open(currdir + '/' + args.outprefix + '_Read2aligned.txt','w')
        GAP_OPEN   = -5 # originally -5
        GAP_EXTEND = -1
	for line, seq in enumerate(seqs_read1):
                seq_read1 = seqs_read1[ line ]
                pos = seq_read1.find( 'AGATCGGAAGAGC' ) # position of ligation adapter. Easy to recognize.
                if ( pos > -1 ): seq_read1 = seq_read1[ :pos];
                WTrev_trunc = WTrev[:pos]
		aligned_read1 = nw.global_align( WTrev_trunc, seq_read1, gap_open=GAP_OPEN, gap_extend=GAP_EXTEND )

                pos = aligned_read1[1].find( '----' )
                if ( pos > -1 ): aligned_read1 = (aligned_read1[0][ :pos], aligned_read1[1][ :pos ] )
		seqs_read1_align.append( aligned_read1 )
		f_seqs_read1.write( aligned_read1[1]+' '+aligned_read1[0] )
		f_seqs_read1.write('\n')

		aligned_read2 = nw.global_align( WTfwd_trunc, seqs_read2[line], gap_open=GAP_OPEN, gap_extend=GAP_EXTEND )
		seqs_read2_align.append( aligned_read2 )
		f_seqs_read2.write( aligned_read2[1]+' '+aligned_read2[0])
		f_seqs_read2.write('\n')
	f_log.write( '\nAlignment finished at ' + timeStamp() )
	f_seqs_read1.close()
	f_seqs_read2.close()

	# Compare WT sequence to both read1 and read2 and build simple file
	simple = [[0 for col in range(WTlen)] for row in range(len(seqs_read1_align))]			# initialize simple array for recording mutations per position per read in binary
	mutations = np.zeros([len(seqs_read1_align),20,WTlen])
	mutdict = {'AT':0 ,'AG':1 ,'AC':2 ,'AN':3 ,'A-':4,
                   'TA':5 ,'TG':6 ,'TC':7 ,'TN':8 ,'T-':9,
                   'GA':10,'GT':11,'GC':12,'GN':13,'G-':14,
                   'CA':15,'CT':16,'CG':17,'CN':18,'C-':19}
			# Dictionary includes all potential mutations, e.g. A to T, G, C, del, or N (ambiguous read, not counted); del = '-'
	mut_projection = np.zeros([20,WTlen])

        print_out_problem = False
        for line, seq_align in enumerate(seqs_read1_align):
                idx = -1
                if line == 1 : print seq_align[0]+'\n'+seq_align[1]+'\n'
                for count,wt_nt in enumerate( seq_align[0] ): # should be wild type sequence
                        if wt_nt != '-': idx += 1
                        read_nt = seq_align[1][count]
                        if read_nt == wt_nt:
                                simple[line][WTlen-idx-1] = 0
                        else:
                                if read_nt == 'N':										# Ignore reads with N because not necessarily mutated
                                        simple[line][WTlen-idx-1] = 0
                                else:
                                        if ( wt_nt != '-' ): # note recording insertions -- fix this!
                                                if ( line == 1 ): print 'mismatch: ', idx, wt_nt, '-->', read_nt
                                                simple[line][WTlen-idx-1] = 1
                                                pair = wt_nt+read_nt
                                                mutations[line][mutdict[pair]][WTlen-idx-1] += 1	# Record 2D array of mutation frequencies for each sequence in Read 1

                if not print_out_problem and simple[line].count(1) > 20:
                     print_out_problem = True
                     print 'Problem line: ', line+1, ' has ', simple[line].count(1), ' mutations'
                     print seq_align[0]+'\n'+seq_align[1]+'\n'

        if False: # off for now.
                for idx, nt in enumerate(WTfwd_trunc):
                        for line, seq in enumerate(seqs_read2_align):
                                if seq[idx] == nt:
                                        simple[line][idx] = 0
                                elif seq[idx] != nt:
                                        if seq[idx] == 'N':
                                                simple[line][idx] = 0
                                        else:
                                                simple[line][idx] = 1
                                                pair = nt+seq[idx]
                                                mutations[line][mutdict[pair]][idx] += 1			# Record 2D array of mutation frequencies for each sequence in Read 2

	# Project mutation frequencies across sequences, filtering out reads with 10 or more
	count = 0
	mut_count = np.sum(mutations, axis=(1,2))
	for line, seq in enumerate(seqs_read2_align):
		# print np.sum(mutations, axis=(1,2))
		if mut_count[line] < 10:
			mut_projection = mut_projection + mutations[line]
			count += 1
	# print str(count) # should match line "Number of sequences with fewer than 10 mismatches to WT" in .log file (assessed in simple_to_rdat fxn)

	mut_projection = np.concatenate([mut_projection, np.sum(mut_projection, axis=1, keepdims=True)], axis=1)			# Total number of each type of mutation across all positions
	mut_projection = np.concatenate([mut_projection, np.sum(mut_projection, axis=0, keepdims=True) / count], axis=0)	# Mutation frequency at each position

	f_muts = open(currdir + '/' + args.outprefix + '_mutations.csv','w')
	f_muts.write('Nucleotide:,' + ','.join(sequence) + '\n')
	f_muts.write('Sequence position:,' + ','.join(map(str,seqpos+1)) + '\n')
	mut_count = 0
	for line in mut_projection:
		if mut_count <= 19:
			revdict = dict((v,k) for k,v in mutdict.iteritems())
			f_muts.write(revdict[mut_count][0] + ' to ' + string.replace(revdict[mut_count][1],'-','del') + ',' + ','.join(map(str,line)) + '\n')
			mut_count += 1
		else:
			f_muts.write('Mutation rate ([total mutations] / [total reads with < 10 mutations])' + ',' + ','.join(map(str,line)) + '\n')
	f_muts.write('Total reads with < 10 mutations:,' + str(count))

	f_simple = open(currdir + '/' + args.outprefix + '.simple','w')
	for line in simple:
		f_simple.write('1' + '\t' + str(WTlen) + '\t')
		f_simple.write(''.join(map(str,line)) + '\n')
	f_simple.close()

	return simple


def simple_to_rdat( args, sequence, sfilein=0, simple=[] ):

	# Setup: get length and seqpos
	WTlen = len(sequence)
	seqpos = arange(0, WTlen) + args.offset
	numreads = 0

	# Get mutation indices from simple data
	if sfilein == 0:											# If simple array was generated by fastq_to_simple
		mut_idxs = []
		for line, dat in enumerate(simple):
			numreads = numreads + 1
			mut_idxs.append(array([idx for idx, val in enumerate(dat) if val == 1]))

	else:														# If .simple file was input
		max_end_pos = -inf
		min_start_pos =  inf
		mut_idxs = []
		for line in args.simplefile.readlines():
			numreads = numreads + 1
			fields = line.split('\t')
			start_pos = int(fields[0])
			end_pos = int(fields[1])
			if end_pos > max_end_pos:
				max_end_pos = end_pos
			if start_pos < min_start_pos:
				min_start_pos = start_pos
			mut_idxs.append(array([idx + start_pos - 1 for idx, char in enumerate(fields[2]) if char == '1']))
			simple.append(fields[2])

	# Set up arrays for 1D and 2D data
	WTdata = np.zeros((1, WTlen))
	data2d = np.zeros(([len(seqpos),len(seqpos)]))

	# Build data arrays
	count = 0
	for indices in mut_idxs:
		# indices -= seqpos[0]									# Adjust sequence position by starting sequence position
		if len(indices) >= 10:
			continue
		elif len(indices) < 10:
			count += 1
			for mutpos in indices:
				data2d[mutpos, indices] += 1 					# Build 2D dataset
				WTdata[0, indices] += 1 						# Build 1D profile
	f_log.write( '\n\nTotal number of sequences: ' + str(numreads) )
	f_log.write( '\nNumber of sequences with fewer than 10 mismatches to WT (used to get 2D data): ' + str(count) + '\n' )

	# Normalize by total signal in row and store row indices for building RDAT file
	row_indices = []
	for row_idx in xrange(data2d.shape[0]):
		if data2d[row_idx,:].sum() > 0*data2d.shape[1]:
			data2d[row_idx, :] /= data2d[row_idx, :].sum()
		row_indices.append(row_idx)
	WTdata /= WTdata.sum()

	# Output data and annotations to RDAT file
	construct = args.name
	sequence_RNA = string.replace(sequence,'T','U')
	structure = ('.'*len(sequence))
	offset = args.offset - 1
	version = 0.34
	filename = currdir + '/' + args.outprefix + '.rdat'
	if sfilein == 0:
		comments = 'Created from %s and %s using fastq_to_rdat.py\n' % (args.read1fastq.name, args.read2fastq.name)
	else:
		comments = 'Created from %s using fastq_to_rdat.py\n' % (args.simplefile.name)
	annotations = {'experimentType:MutationalProfiling'}
	data_annotations = []
	data_annotations.append({'mutation':'WT'})
	for idx, row_idx in enumerate(row_indices):
		position = seqpos[row_idx] - args.offset
		data_annotations.append({'mutation':[string.replace(sequence[position],'T','U') + str(position + args.offset) + string.replace(reverse_complement(sequence[position]),'T','U')]})
	data = np.concatenate([WTdata, data2d], axis=0)			# combine WT data and 2D data

	# np.set_printoptions(threshold=np.nan)
	# print str(data[0,:])

	r = RDATFile()
	r.save_a_construct(construct, data, sequence_RNA, structure, offset, annotations, data_annotations, filename, comments, version)


	f_log.write( '\n\nRDAT file created: ' + currdir + '/' + args.outprefix + '.rdat')
	f_log.close()

	return



if __name__ == '__main__':

	# Get current directory and open log file
	currdir = os.getcwd()
	f_log = open(currdir + '/' + args.outprefix + '.log', 'w')
	f_log.write( 'Current directory: ' + currdir + '\n\n' )

	# See if user input .simple data file or FASTQ files and run appropriate commands
	if args.simplefile is not None:
		f_log.write( 'Simple data file provided: ' + args.simplefile.name + '\nGenerating 2D data from simple data file.\n\n' )
		simple_to_rdat( args, sequence, 1 )

	elif args.read1fastq is not None and args.read2fastq is not None:
		f_log.write( 'FASTQs provided:\t(Read 1) ' + args.read1fastq.name + '\t(Read 2) ' + args.read2fastq.name + '\nGenerating simple data file from FASTQs.\n\n')
		simple = fastq_to_simple( args )
		simple_to_rdat( args, sequence, 0, simple )

	else:
		f_log.write( 'No data files provided.\n' )
		exit()







