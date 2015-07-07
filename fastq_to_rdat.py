#!/usr/bin/python
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
import sys
import string
from rdatkit.datahandlers import RDATFile
from matplotlib.pylab import *

parser = argparse.ArgumentParser( description='Generates a 2D correlated mutational profiling dataset using demultiplexed Read1 and Read 2 files or a .simple file')

parser.add_argument('sequencefile', type=argparse.FileType('r'), help='name of sequence file in .fasta format')
parser.add_argument('--read1fastq', type=argparse.FileType('r'), help='name of Illumina R1 fastq file')
parser.add_argument('--read2fastq', type=argparse.FileType('r'), help='name of Illumina R2 fastq file')
parser.add_argument('--simplefile', type=argparse.FileType('r'), help='name of .simple file if you have it; otherwise must specify fastq files.')
parser.add_argument('--name', type=str, help='name of RNA; by default take from sequencefile',default='PLACEHOLDER')
parser.add_argument('--offset', type=int, help='integer to add to residue numbers to get conventional numbers', default=0)
parser.add_argument('--outprefix', type=str, default='out')
parser.add_argument('--num_hits_cutoff', type=int, help='quality filter: maximum number of hits to allow before recording.', default=10)
parser.add_argument('--start_pos_cutoff', type=int, help='full-length filter: minimal nucleotide position to which alignment must extend',default=10)
parser.add_argument('--collapse_adjacent_mutations', type=bool, help='if several mutations/indels occur in a row, count them as a single mutation occurring at last position.',default=True)
parser.add_argument('--ignore_dels', type=bool, help='do not record deletions in simple/rdat files.',    default=False)
parser.add_argument('--ignore_inserts', type=bool, help='do not record insertions in simple/rdat files.',default=True)
parser.add_argument('--max_reads', type=int, help='maximum number of reads to take from FASTQ',default=0)

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


# Dictionary includes all potential mutations, e.g. A to T, G, C, del, or N (ambiguous read, not counted); del = '-'
muts = [ 'AT','AG','AC','AN','A-',
         'TA','TG','TC','TN','T-',
         'GA','GT','GC','GN','G-',
         'CA','CT','CG','CN','C-',
         '-A','-C','-G','-T','-N', # insertions
         'sequencing_error1', 'sequencing_error2' ]

mutdict = {}
for i,mut in enumerate(muts): mutdict[ mut ] = i

def record_mutations( seq_align, line, simple, mutations, match,  direction = 'reverse'):
        idx = -1
        for count,(wt_nt,read_nt) in enumerate( zip( seq_align[0], seq_align[1] ) ):
                if wt_nt != '-': idx += 1
                idx_to_use = idx
                if (direction == 'reverse'): idx_to_use = len( simple[0] ) - idx -1
                if ( idx_to_use >= len( simple[0] ) ): continue
                if read_nt == wt_nt:
                        # check for errors -- but skip indels -- this is really awful. Should really align READ1 and READ2 to each other *first*, then align to reference sequence.
                        if simple[line][idx_to_use] == 1 and \
                           not( mutations[line][mutdict['A-']][idx_to_use] or mutations[line][mutdict['C-']][idx_to_use] or mutations[line][mutdict['T-']][idx_to_use] or mutations[line][mutdict['G-']][idx_to_use]) and \
                           not( mutations[line][mutdict['-A']][idx_to_use] or mutations[line][mutdict['-C']][idx_to_use] or mutations[line][mutdict['-T']][idx_to_use] or mutations[line][mutdict['-G']][idx_to_use]):
                                mutations[line][mutdict['sequencing_error1']][idx_to_use] += 1	# Error in read1 if read2 gets matched
                        simple[line][idx_to_use] = 0
                        match [line][idx_to_use] = 1	# Record a match.
                else:
                        if read_nt == 'N':										# Ignore reads with N because not necessarily mutated
                                simple[line][idx_to_use] = 0
                        else:
                                pair = wt_nt+read_nt
                                # note that mutation count are currently crazy -- read1 and read2 are separate -- indels always recorded in mutation, but not in simple
                                mutations[line][mutdict[pair]][idx_to_use] += 1	# Record 2D array of mutation frequencies for each sequence in Read 1
                                if args.ignore_dels    and read_nt == '-': continue
                                if args.ignore_inserts and wt_nt   == '-': continue
                                # check if matched before -- record as a sequencer error. currently not reverting indels -- see note above.
                                if match[line][idx_to_use] == 1 and (read_nt != '-') and (wt_nt != '-'):
                                        #if ( line == 1 ) :  print line,idx_to_use, wt_nt, read_nt
                                        assert( simple[line][idx_to_use] == 0 )
                                        mutations[line][mutdict['sequencing_error2']][idx_to_use] += 1	# Error in read2 if read1 was matched
                                else:
                                        simple[line][idx_to_use] = 1
                                        match [line][idx_to_use] = 0


def truncate_junk_at_end( aligned_read1 ):
        q = aligned_read1[ 0 ]
        r = aligned_read1[ 1 ]
        MIN_MATCH = 4 # number of consecutive nts that must match to determine if we are in a 'good' region.
        pos = len( q )
        while pos > (MIN_MATCH-1):
                if ( q[ pos-MIN_MATCH : pos ] == r[ pos-MIN_MATCH : pos ] ):
                        break
                pos -= 1
        aligned_read_truncate = ( q[:pos],r[:pos] )
        num_junk_nts = 0
        for char in r[pos:]:
                if char != '-': num_junk_nts += 1

        maxpos_wt = 0
        for char in q[:pos]:
                if char != '-': maxpos_wt += 1

        return ( aligned_read_truncate, num_junk_nts, maxpos_wt )

def collapse_adjacent_mutations( simple ):
        for line, simple_line in enumerate(simple):
                last_mismatch = -1
                for pos, simple_symbol in enumerate(simple_line):
                        if simple_symbol == 1:
                                if last_mismatch == -1:
                                        last_mismatch = pos # entering a mismatch
                        else:
                                if last_mismatch > -1:
                                        for idx in range( last_mismatch, pos-1 ): simple[ line ][ idx ] = 0 # put mutation on right-most nucleotide -- that's where RT messed up.
                                        last_mismatch = -1 # not in a mismatch anymore

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
        # Probably should not load these all into memory at once.
	seqs_read1 = []
	for i, line in enumerate(args.read1fastq):
                if ( i % 40000 == 0 ): print 'Read in line ', i/4
		if i % 4 == 1:
			seqs_read1.append( line.strip() )
                if args.max_reads > 0 and len( seqs_read1 ) >= args.max_reads: break
        print 'Finished reading in file for Read 1. Number of seqs:', len( seqs_read1)

	seqs_read2 = []
        for i, line in enumerate(args.read2fastq):
		if i % 4 == 1:
			seqs_read2.append( line.strip() )
                if args.max_reads > 0 and len( seqs_read2 ) >= args.max_reads: break

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
        start_pos        = []
	f_seqs_read1 = open(currdir + '/' + args.outprefix + '_Read1aligned.txt','w')
	f_seqs_read2 = open(currdir + '/' + args.outprefix + '_Read2aligned.txt','w')
        NW_GAP_OPEN   = -5 # originally -5
        NW_GAP_EXTEND = -1 # originally -1
	for line, (seq_read1,seq_read2) in enumerate(zip(seqs_read1,seqs_read2)):
                if ( line % 10000 == 0): print 'Doing alignment for line ', line, ' out of ', len( seqs_read1 )
                # A somewhat better workflow here might be to align READ1 and READ2 to each other first,
                # detect and remove sequencing errors, then do a local alignment of this high-quality sequence into WT.
                # do READ1 (reverse read)
                maxpos1     = seq_read1.find( 'AGATCGGAAGAGC' ) # position of ligation adapter. Easy to recognize. Should not hardcode this in, though.
                if ( maxpos1 == -1 ): maxpos1 = WTlen  # maybe read was not long enough to see adapter
                # truncation would not be necessary if we used a local align:
                seq_read1   = seq_read1[:maxpos1];
                WTrev_trunc = WTrev[    :maxpos1]
                aligned_read1 = nw.global_align( WTrev_trunc, seq_read1, gap_open=NW_GAP_OPEN, gap_extend=NW_GAP_EXTEND )

                # Following edits aligned_read1 to remove 'junk' nts that SSIII tacks onto the ends of transcripts
                # note that maxpos1 is updated (position at which reverse read stops, in WTrev numbering)
                ( aligned_read1, num_junk_nts, maxpos1 ) = truncate_junk_at_end( aligned_read1 )

		seqs_read1_align.append( aligned_read1 )
		f_seqs_read1.write( aligned_read1[0]+'\n'+aligned_read1[1]+'\n\n' )

                # do forward read. Have some information for where it starts based on where ligation site showed up in read1.
                seq_read2   = seq_read2[ num_junk_nts : ] # remove the junk (was 3' in read1; is 5' in read2)
                maxpos2     = WTlen - maxpos1 + len( seq_read2 ) # where read2 should end in WTfwd
                WTfwd_trunc = WTfwd[      WTlen - maxpos1    : maxpos2           ] # sometimes maxpos2 is beyond WTlen -- then this is slightly shorter than seq_read2
                seq_read2   = seq_read2[  : len( WTfwd_trunc) ]
                aligned_read2 = nw.global_align( WTfwd_trunc, seq_read2, gap_open=NW_GAP_OPEN, gap_extend=NW_GAP_EXTEND )

                # In case RT stopped before 5' end of RNA (maxpos1 < WTlen):
                # need to pad? actually should get rid of this, and just keep track of start_pos
                pad_sequence = ''
                for k in range( WTlen - maxpos1 ): pad_sequence += 'N'
                aligned_read2 = ( pad_sequence+aligned_read2[0],  pad_sequence+aligned_read2[1] )
                # should be under an option? these are not actually junk nts, but errors arising at edges of global alignment
                # just used here for trimming aligned_read2
                ( aligned_read2, num_junk_nts2, maxpos2 ) = truncate_junk_at_end( aligned_read2 )

		seqs_read2_align.append( aligned_read2 )
		f_seqs_read2.write( aligned_read2[0]+'\n'+aligned_read2[1]+'\n\n')
                start_pos.append( WTlen - maxpos1 + 1 ) # for 'simple' output.

	f_log.write( '\nAlignment finished at ' + timeStamp() )
	f_seqs_read1.close()
	f_seqs_read2.close()

	# Compare WT sequence to both read1 and read2 and build simple file
	simple = [[0 for col in range(WTlen)] for row in range(len(seqs_read1_align))]			# initialize simple array for recording mutations per position per read in binary
	mutations = np.zeros([len(seqs_read1_align),len(mutdict),WTlen])
        match = np.zeros([len(seqs_read1_align),WTlen])
	mut_projection = np.zeros([len(muts),WTlen])

        print_out_problem = 0
        for line, seq_align in enumerate(seqs_read1_align):
                if ( line % 10000 == 0): print 'Doing mutation assignment for line ', line, ' out of ', len( seqs_read1_align )
                record_mutations( seq_align, line, simple, mutations, match, direction = 'reverse' )

        for line, seq_align in enumerate(seqs_read2_align):
                if ( line % 10000 == 0): print 'Doing mutation assignment for line ', line, ' out of ', len( seqs_read2_align )
                record_mutations( seq_align, line, simple, mutations, match, direction = 'forward' )

        # If several mutations/indels occur in a row, count them as a single mutation occurring at last position.
        if args.collapse_adjacent_mutations: collapse_adjacent_mutations( simple )

	# Project mutation frequencies across sequences, filtering out reads with 10 or more
	count = 0
	mut_count = np.sum(mutations[:,:-1,:], axis=(1,2) )
	for line, seq in enumerate(seqs_read1_align):
		# print np.sum(mutations, axis=(1,2))
		if mut_count[line] <= args.num_hits_cutoff and start_pos[line] <= args.start_pos_cutoff:
			mut_projection = mut_projection + mutations[line]
			count += 1

        # print str(count) # should match line "Number of sequences with fewer than 10 mismatches to WT" in .log file (assessed in simple_to_rdat fxn)

	mut_projection = np.concatenate([mut_projection, np.sum(mut_projection, axis=1, keepdims=True)], axis=1)	        # Total number of each type of mutation across all positions
        mut_projection = np.concatenate([mut_projection, np.sum(mut_projection, axis=0, keepdims=True) / count], axis=0)	# Mutation frequency at each position

	f_muts = open(currdir + '/' + args.outprefix + '_mutations.csv','w')
	f_muts.write('Nucleotide:,' + ','.join(sequence) + '\n')
	f_muts.write('Sequence position:,' + ','.join(map(str,seqpos+1)) + '\n')
        revdict = dict((v,k) for k,v in mutdict.iteritems())
	for mut_count,line in enumerate( mut_projection ):
		if mut_count < len( muts ):
                        tag = revdict[mut_count]
                        if len( tag ) == 2: tag = string.replace(revdict[mut_count][0],'-','del') + ' to ' + string.replace(revdict[mut_count][1],'-','del')
			f_muts.write( tag + ',' + ','.join(map(str,line)) + '\n')
		else:
			f_muts.write('Mutation rate ([total mutations] / [total reads with < '+str(args.num_hits_cutoff)+' mutations])' + ',' + ','.join(map(str,line)) + '\n')
	f_muts.write('Total reads with < '+str(args.num_hits_cutoff)+' mutations:,' + str(count))

        # could pre-filter by num_hits or start_pos here, to reduce disk output.
	f_simple = open(currdir + '/' + args.outprefix + '.simple','w')
	for line, simple_line in enumerate(simple):
		f_simple.write(str(start_pos[ line ]) + '\t' + str(WTlen) + '\t')
		f_simple.write(''.join(map(str,simple_line[ start_pos[line]-1 : ] )) + '\n')
	f_simple.close()

        output_string = '\n\nsimple file created: ' + currdir + '/' + args.outprefix + '.simple'
        print output_string

	return ( simple, start_pos )



def output_rdat( filename, args, sequence, row_indices, seqpos, WTdata, data2d, f_log, sfilein ):
        # Output data and annotations to RDAT file
	construct = args.name
	sequence_RNA = string.replace(sequence,'T','U')
	structure = ('.'*len(sequence))
	offset = args.offset
	version = 0.34
	if sfilein == 0:
		comments = 'Created from %s and %s using fastq_to_rdat.py\n' % (args.read1fastq.name, args.read2fastq.name)
	else:
		comments = 'Created from %s using fastq_to_rdat.py\n' % (args.simplefile.name)
	annotations = {'experimentType:MutationalProfiling'}
	data_annotations = []
	data_annotations.append({'mutation':'WT'})
	for idx, row_idx in enumerate(row_indices):
		position = seqpos[row_idx] - offset
		data_annotations.append({'mutation':[string.replace(sequence[position-1],'T','U') + str(position + offset) + 'X']})
	data = np.concatenate([WTdata, data2d], axis=0)			# combine WT data and 2D data
	# np.set_printoptions(threshold=np.nan)
	# print str(data[0,:])

	r = RDATFile()
        # what's the deal with offset? MATLAB scripts add this number to 1, 2, ... N. But Python scripts add this number to 0, 1, ... N-1?
	r.save_a_construct(construct, data, sequence_RNA, structure, offset, annotations, data_annotations, filename, comments, version)

        output_string = '\n\nRDAT file created: ' + currdir + '/' + args.outprefix + '.rdat'
        sys.stdout.write( output_string + '\n' )
	f_log.write( output_string )

def simple_to_rdat( args, sequence, sfilein=0, simple=[], start_pos=[] ):

	# Setup: get length and seqpos
	WTlen = len(sequence)
	seqpos = arange(1, WTlen+1) + args.offset
	# Get mutation indices from simple data
	if sfilein != 0:
                assert( len(simple) == 0 )
                assert( len(start_pos) == 0 )
                print 'Reading file: '+str(args.simplefile.name)
		for line in args.simplefile.readlines():
			fields = line.split('\t')
                        if len( fields ) < 3: continue
			start_pos.append( int(fields[0]) )
                        simple_string = fields[2][:-1] # no endline
                        assert( len(simple_string) == (int( fields[1]) - int(fields[0]) + 1 ) )
			simple.append( [ int(x) for x in simple_string ] )

        mut_idxs = []
        for line, dat in enumerate(simple):
                if ( line % 10000 == 0): print 'Doing simple-to-RDAT mutation assignment for line ', line, ' out of ', len( simple )
                start_pos_to_use = 1
                if ( len( start_pos ) > 0 ):
                        start_pos_to_use = start_pos[ line ]
                        if ( start_pos[ line ] > args.start_pos_cutoff ): continue
                mut_idx = array( [ (idx + start_pos_to_use - 1) for idx, val in enumerate(dat) if val == 1])
                if ( len( mut_idx ) > args.num_hits_cutoff ): continue
                mut_idxs.append(mut_idx)
        numreads = len( mut_idxs )

	# Set up arrays for 1D and 2D data
	WTdata = np.zeros((1, WTlen))
	data2d = np.zeros(([len(seqpos),len(seqpos)]))

	# Build data arrays
	count = 0
	for (line,indices) in enumerate( mut_idxs ):
                if ( line % 10000 == 0): print 'Doing data2d compilation for line  ', line, ' out of ', len( mut_idxs )
		# indices -= seqpos[0]									# Adjust sequence position by starting sequence position
                assert( len(indices) <= args.num_hits_cutoff )
                count += 1
                if len( indices ) > 0:
                        WTdata[ 0, indices] += 1 					# Build 1D profile
                for mutpos in indices:
                        data2d[mutpos, indices] += 1 					# Build 2D dataset
	f_log.write( '\n\nTotal number of sequences: ' + str(numreads) )
	f_log.write( '\nNumber of sequences with fewer than 10 mismatches to WT (used to get 2D data): ' + str(count) + '\n' )

        # normalizes to total number of reads -- should be a true 'modification fraction'
        WTdata_reactivity = WTdata/count
        # normalizes to total number of hits -- original choice
	WTdata_norm = WTdata/WTdata.sum()

	# Normalize by total signal in row and store row indices for building RDAT file
	row_indices = []
        data2d_reactivity = 0 * data2d
        data2d_norm       = 0 * data2d
	for row_idx in xrange(data2d.shape[0]):
		if data2d[row_idx,:].sum() > 0*data2d.shape[1]:
                        # normalizes to total number of reads with modification at row_idx position -- should be a true 'modification fraction'
                        data2d_reactivity[row_idx, :] = data2d[row_idx,:] / data2d[ row_idx, row_idx ]
                        # normalizes to total number of hits -- original choice
			data2d_norm[row_idx, :] = data2d[ row_idx, : ] / data2d[row_idx, :].sum()
		row_indices.append(row_idx)

	filename = currdir + '/' + args.outprefix + '.reactivity.rdat'
        output_rdat( filename, args, sequence, row_indices, seqpos, WTdata_reactivity, data2d_reactivity, f_log, sfilein )

	filename = currdir + '/' + args.outprefix + '.rdat'
        output_rdat( filename, args, sequence, row_indices, seqpos, WTdata_norm, data2d_norm, f_log, sfilein )

	f_log.close()

	return


if __name__ == '__main__':

        # get target directory ready.
        targetdir = os.path.dirname( args.outprefix )
        if ( len( targetdir ) > 0 and not os.path.exists( targetdir ) ): os.makedirs( targetdir )

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
                ( simple, start_pos ) = fastq_to_simple( args )
		simple_to_rdat( args, sequence, 0, simple, start_pos )

	else:
		f_log.write( 'No data files provided.\n' )
		exit()







