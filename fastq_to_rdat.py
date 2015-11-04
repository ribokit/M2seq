#!/usr/bin/python
#######################################################################################
# Analysis of 2D signal in mutational profiling sequencing data (MaP2D)
#######################################################################################
#
# FASTQ to RDAT
#
# This script generates a 2D correlated mutational profiling dataset using one of three sources:
#        1. Demultiplexed and quality-filtered Read 1 and Read 2 FASTQ files
#        2. Text files containing Read 1 and Read 2 sequences aligned to WT sequence of interest
#        3. Data in 'simple' binary format that records the mutations in sequenced cDNAs
#
# The outputs of the script are:
#        1. An RDAT file containing the MaP2D (two-dimensional mutational profiling) data
#        2. A 'simple' file that records the mutations in each cDNA in binary format
#        3. A text file for each input FASTQ showing all sequences after alignment ('...Read1Aligned.txt' and '...Read2Aligned.txt'). These aligned text files can be used as inputs to this script.
#        4. A '...mutations.csv' file quantifying the numbers of mutations and mutation rates
#        5. A .log file recording minor output from the analysis, e.g. input files
#        NOTE: If aligned reads files are input, only 1, 2, 4, and 5 are generated.
#        NOTE: If a 'simple' file is input, only 1 and 5 are generated.
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
from itertools import izip_longest
from rdatkit.datahandlers import RDATFile
from matplotlib.pylab import *

parser = argparse.ArgumentParser( description='Generates a 2D correlated mutational profiling dataset using demultiplexed Read1 and Read 2 files or a .simple file')

parser.add_argument('sequencefile', type=argparse.FileType('r'), help='name of sequence file in .fasta format')
parser.add_argument('--read1fastq', type=argparse.FileType('r'), help='name of Illumina R1 fastq file')
parser.add_argument('--read2fastq', type=argparse.FileType('r'), help='name of Illumina R2 fastq file')
parser.add_argument('--read1align', type=argparse.FileType('r'), help='name of aligned reads from Illumina R1 fastq file')
parser.add_argument('--read2align', type=argparse.FileType('r'), help='name of aligned reads from Illumina R2 fastq file')
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


def timeStamp():        # From ShapeMapper
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


def record_mutations( seq_align, line, simple_line, mutations_line, match_line,  direction = 'reverse'):
        idx = -1
        for count,(wt_nt,read_nt) in enumerate( zip( seq_align[0], seq_align[1] ) ):
                if wt_nt != '-': idx += 1
                idx_to_use = idx
                if (direction == 'reverse'): idx_to_use = len( simple_line ) - idx -1
                if ( idx_to_use >= len( simple_line ) ): continue
                if read_nt == wt_nt:
                        # check for errors -- but skip indels -- this is really awful. Should really align READ1 and READ2 to each other *first*, then align to reference sequence.
                        if simple_line[idx_to_use] == 1 and \
                           not( mutations_line[mutdict['A-']][idx_to_use] or mutations_line[mutdict['C-']][idx_to_use] or mutations_line[mutdict['T-']][idx_to_use] or mutations_line[mutdict['G-']][idx_to_use]) and \
                           not( mutations_line[mutdict['-A']][idx_to_use] or mutations_line[mutdict['-C']][idx_to_use] or mutations_line[mutdict['-T']][idx_to_use] or mutations_line[mutdict['-G']][idx_to_use]):
                                mutations_line[mutdict['sequencing_error1']][idx_to_use] += 1               # Error in read1 if read2 gets matched
                        simple_line[idx_to_use] = 0
                        match_line[idx_to_use]  = 1         # Record a match.
                else:
                        if read_nt == 'N':                                                                  # Ignore reads with N because not necessarily mutated
                                simple_line[idx_to_use] = 0
                        else:
                                pair = wt_nt+read_nt
                                # note that mutation count are currently crazy -- read1 and read2 are separate -- indels always recorded in mutation, but not in simple
                                mutations_line[mutdict[pair]][idx_to_use] += 1                              # Record 2D array of mutation frequencies for each sequence in Read 1
                                if args.ignore_dels    and read_nt == '-': continue
                                if args.ignore_inserts and wt_nt   == '-': continue
                                # check if matched before -- record as a sequencer error. currently not reverting indels -- see note above.
                                if match_line[idx_to_use] == 1 and (read_nt != '-') and (wt_nt != '-'):
                                        #if ( line == 1 ) :  print line,idx_to_use, wt_nt, read_nt
                                        assert( simple_line[idx_to_use] == 0 )
                                        mutations_line[mutdict['sequencing_error2']][idx_to_use] += 1       # Error in read2 if read1 was matched
                                else:
                                        simple_line[idx_to_use] = 1
                                        match_line[idx_to_use]  = 0


def record_muts_generate_simple( aligned_read1, aligned_read2, f_simple, mut_projection, seqnm, fltnm, WTlen, start_pos ):
        
        # Initialize variables for mutation counting
        simple_line    = [0 for col in range(WTlen)]                # initialize simple line for recording mutations per position, in binary
        mutations_line = np.zeros([len(mutdict),WTlen])
        match_line     = np.zeros([WTlen])

        # Record mutations
        record_mutations( aligned_read1, seqnm, simple_line, mutations_line, match_line, direction = 'reverse' )
        record_mutations( aligned_read2, seqnm, simple_line, mutations_line, match_line, direction = 'forward' )
        
        # If several mutations/indels occur in a row, count them as a single mutation occurring at last position.
        if args.collapse_adjacent_mutations: collapse_adjacent_mutations( simple_line )
        
        # Truncate simple line to start at start_pos
        simple_line = simple_line[ (start_pos - 1): ]

        # Write simple line to file
        f_simple.write(str(start_pos) + '\t' + str(WTlen) + '\t')
        f_simple.write(''.join(map( str,simple_line )) + '\n')

        # Project mutation frequencies across sequences, filtering out reads with 10 or more
        mut_count = np.sum(mutations_line[:-1,:], axis=(0,1) )      # gets total mutations for each sequence (?)
        if mut_count <= args.num_hits_cutoff and start_pos <= args.start_pos_cutoff:
            mut_projection = mut_projection + mutations_line
            fltnm += 1

        return ( fltnm, mut_projection )


def truncate_junk_at_end( aligned_read1 ):
        q = aligned_read1[ 0 ]
        r = aligned_read1[ 1 ]
        MIN_MATCH = 4       # number of consecutive nts that must match to determine if we are in a 'good' region.
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


def collapse_adjacent_mutations( simple_line ):
        last_mismatch = -1
        for pos, simple_symbol in enumerate(simple_line):
                if simple_symbol == 1:
                        if last_mismatch == -1:
                                last_mismatch = pos # entering a mismatch
                else:
                        if last_mismatch > -1:
                                for idx in range( last_mismatch, pos-1 ): simple_line[ idx ] = 0    # put mutation on right-most nucleotide -- that's where RT messed up.
                                last_mismatch = -1  # not in a mismatch anymore


def fastq_to_simple( args ):            # Build 'simple' array with 1 at each mutation and 0 at each WT nt, starting from fastq files
        print 'Starting analysis from FASTQs at ' + timeStamp()
        f_log.write( 'Starting analysis from FASTQs at ' + timeStamp() )
        f_log.write( '\n\nRNA name: ' + args.name )
        
        # Get WT sequence, rev comp, and length
        WTfwd = sequence
        WTrev = reverse_complement(WTfwd)
        WTlen = len(WTfwd)
        seqpos = arange(0, WTlen) + args.offset + 1
        f_log.write( '\nLength of WT sequence: ' + str(WTlen) + '\n' )
        f_log.write( '\n\nLength of seqpos: ' + str(len(seqpos)) )
        f_log.write( '\nWild-type sequence for matching to Read 2:\n' + WTfwd + '\n' )
        f_log.write( '\nReverse complement for matching to Read 1:\n' + WTrev + '\n' )
        
        # INITIALIZE VARIABLES
        start_pos = []
        NW_GAP_OPEN   = -5
        NW_GAP_EXTEND = -1
        f_seqs_read1 = open(currdir + '/' + args.outprefix + '_Read1aligned.txt','w')
        f_seqs_read2 = open(currdir + '/' + args.outprefix + '_Read2aligned.txt','w')
        f_simple     = open(currdir + '/' + args.outprefix + '.simple','w')

        read1 = args.read1fastq
        read2 = args.read2fastq
        count = 0
        seqnm = 0
        fltnm = 0
        mut_projection = np.zeros([len(muts),WTlen])

        for x,y in izip_longest(read1, read2):      # iterate through read1 and read2 (with both sequences read in first)

            if x is None or y is None: print 'Number of sequences in read 1 and read 2 FASTQs are unequal!'; exit()     # check if files are different length
            if args.max_reads > 0 and count >= args.max_reads: break                                                    # check if gone beyond max reads, if max reads input by user

            if count % 4 == 1:
                seq_read1 = x.strip()
                seq_read2 = y.strip()
                
                # Truncate WT fwd and rev sequences to length of read, in case read is shorter than probed sequence
                if count == 1: WTrev_trunc = WTrev[0:len(seq_read1)]; WTfwd_trunc = WTfwd[0:len(seq_read2)];

                ############################################################
                # ALIGN READ1 AND READ2 TO WT (in future, should ALIGN READ1 and READ2 to each other as FULLREAD (detecting/removing sequencing errors), then LOCALLY ALIGN FULLREAD to WT - but for now, still do each read separately)
                maxpos1     = seq_read1.find( 'AGATCGGAAGAGC' )         # position of ligation adapter. Easy to recognize. Should not hardcode this in, though.
                if ( maxpos1 == -1 ): maxpos1 = WTlen - 1               # maybe read was not long enough to see adapter

                # Align read1 (reverse read)
                # truncation would not be necessary if we used a local align:
                seq_read1   = seq_read1[ :maxpos1 ]
                WTrev_trunc = WTrev[ :maxpos1 ]
                aligned_read1 = nw.global_align( WTrev_trunc, seq_read1, gap_open=NW_GAP_OPEN, gap_extend=NW_GAP_EXTEND )
                ( aligned_read1, num_junk_nts, maxpos1 ) = truncate_junk_at_end( aligned_read1 )        # Remove 'junk' nts that SSIII tacks onto ends of transcripts; note that maxpos1 is updated (position at which reverse read stops, in WTrev numbering)
                f_seqs_read1.write( aligned_read1[0]+'\n'+aligned_read1[1]+'\n\n' )

                # Align read2 (forward read)
                # Have some information for where it starts based on where ligation site showed up in read1.
                seq_read2   = seq_read2[ num_junk_nts: ]                # remove the junk (was 3' in read1; is 5' in read2)
                maxpos2     = WTlen - maxpos1 + len( seq_read2 )        # where read2 should end in WTfwd
                WTfwd_trunc = WTfwd[ WTlen - maxpos1:maxpos2 ]          # sometimes maxpos2 is beyond WTlen -- then this is slightly shorter than seq_read2
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
                f_seqs_read2.write( aligned_read2[0]+'\n'+aligned_read2[1]+'\n\n')
                
                start_pos = WTlen - maxpos1 + 1 # for 'simple' output.
                
                ############################################################
                # RECORD MUTATIONS AND GENERATE SIMPLE FILE
                ( fltnm, mut_projection ) = record_muts_generate_simple( aligned_read1, aligned_read2, f_simple, mut_projection, seqnm, fltnm, WTlen, start_pos )

                seqnm += 1
                if seqnm % 50000 == 0: print 'Reading sequence number ', seqnm

            # Advance count
            count += 1

        mut_projection = np.concatenate([mut_projection, np.sum(mut_projection, axis=1, keepdims=True)], axis=1)                # Total number of each type of mutation across all positions
        mut_projection = np.concatenate([mut_projection, np.sum(mut_projection, axis=0, keepdims=True) / fltnm], axis=0)        # Mutation frequency at each position

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
        
        f_muts.write('Total reads with < '+str(args.num_hits_cutoff)+' mutations:,' + str(fltnm))

        simple_path = f_simple.name
        output_string = '\n\nSimple file created: ' + simple_path
        print output_string

        f_seqs_read1.close()
        f_seqs_read2.close()
        f_simple.close()

        f_log.write( '\nTotal number of sequences: ' + str(seqnm) )
        # print str(fltnm) # should match line "Number of sequences with fewer than 10 mismatches to WT" in .log file (assessed in simple_to_rdat fxn)
        f_log.write( '\nFinished analysis from FASTQs at ' + timeStamp() + '\n\n' )

        return ( simple_path )


def aligned_to_simple( args ):          # Build 'simple' array with 1 at each mutation and 0 at each WT nt, starting from aligned sequence text files
        print 'Starting analysis from aligned sequence files at ' + timeStamp()
        f_log.write( 'Starting analysis from aligned sequence files at ' + timeStamp() )
        f_log.write( '\n\nRNA name: ' + args.name )
        
        # Get WT sequence, rev comp, and length
        WTfwd = sequence
        WTrev = reverse_complement(WTfwd)
        WTlen = len(WTfwd)
        seqpos = arange(1, WTlen+1) + args.offset
        f_log.write( '\nLength of WT sequence: ' + str(WTlen) + '\n' )
        f_log.write( '\n\nLength of seqpos: ' + str(len(seqpos)) )
        f_log.write( '\nWild-type sequence for matching to Read 2:\n' + WTfwd + '\n' )
        f_log.write( '\nReverse complement for matching to Read 1:\n' + WTrev + '\n' )
        
        # INITIALIZE VARIABLES
        start_pos = []
        NW_GAP_OPEN   = -5
        NW_GAP_EXTEND = -1
        f_simple     = open(currdir + '/' + args.outprefix + '.simple','w')

        read1 = args.read1align
        read2 = args.read2align
        count = 0
        seqnm = 0
        fltnm = 0
        temp1 = []
        temp2 = []
        mut_projection = np.zeros([len(muts),WTlen])

        for x,y in izip_longest(read1, read2):      # iterate through read1 and read2 (with both sequences read in first)

            if x is None or y is None: print 'Number of sequences in read 1 and read 2 aligned sequence files are unequal!'; exit()     # check if files are different length
            if args.max_reads > 0 and count >= args.max_reads: break                                                                    # check if gone beyond max reads, if max reads input by user

            aligned_read1 = []
            aligned_read2 = []
            
            if count % 3 == 0:
                temp1.append( x.strip() )
                temp2.append( y.strip() )
                maxpos_wt = 0
                for char in temp1[0]:
                    if char != '-': maxpos_wt += 1
                start_pos = WTlen - maxpos_wt + 1   # for 'simple' output
            elif count % 3 == 1:
                temp1.append( x.strip() )
                temp2.append( y.strip() )
            elif count % 3 == 2:
                aligned_read1 = temp1
                aligned_read2 = temp2
                temp1 = []
                temp2 = []

                ############################################################
                # RECORD MUTATIONS AND GENERATE SIMPLE FILE
                ( fltnm, mut_projection ) = record_muts_generate_simple( aligned_read1, aligned_read2, f_simple, mut_projection, seqnm, fltnm, WTlen, start_pos )
                
                seqnm += 1
                if seqnm % 50000 == 0: print 'Reading sequence number ', seqnm

            # Advance count
            count += 1

        mut_projection = np.concatenate([mut_projection, np.sum(mut_projection, axis=1, keepdims=True)], axis=1)                # Total number of each type of mutation across all positions
        mut_projection = np.concatenate([mut_projection, np.sum(mut_projection, axis=0, keepdims=True) / fltnm], axis=0)        # Mutation frequency at each position

        f_muts = open(currdir + '/' + args.outprefix + '_mutations.csv','w')
        f_muts.write('Nucleotide:,' + ','.join(sequence) + '\n')
        f_muts.write('Sequence position:,' + ','.join(map(str,seqpos)) + '\n')
        revdict = dict((v,k) for k,v in mutdict.iteritems())
        for mut_count,line in enumerate( mut_projection ):
                if mut_count < len( muts ):
                        tag = revdict[mut_count]
                        if len( tag ) == 2: tag = string.replace(revdict[mut_count][0],'-','del') + ' to ' + string.replace(revdict[mut_count][1],'-','del')
                        f_muts.write( tag + ',' + ','.join(map(str,line)) + '\n')
                else:
                        f_muts.write('Mutation rate ([total mutations] / [total reads with < '+str(args.num_hits_cutoff)+' mutations])' + ',' + ','.join(map(str,line)) + '\n')
        
        f_muts.write('Total reads with < '+str(args.num_hits_cutoff)+' mutations:,' + str(fltnm))

        simple_path = f_simple.name
        output_string = '\n\nSimple file created: ' + simple_path
        print output_string

        f_simple.close()

        f_log.write( '\nTotal number of sequences: ' + str(seqnm) )
        # print str(fltnm) # should match line "Number of sequences with fewer than 10 mismatches to WT" in .log file (assessed in simple_to_rdat fxn)
        f_log.write( '\nFinished analysis from aligned sequence files at ' + timeStamp() + '\n\n' )

        return ( simple_path )


def output_rdat( filename, args, sequence, row_indices, seqpos, WTdata, data2d, f_log, filesin ):
        # Output data and annotations to RDAT file
        construct = args.name
        sequence_RNA = string.replace(sequence,'T','U')
        structure = ('.'*len(sequence))
        offset = args.offset
        version = 0.34
        if filesin == 0:
            comments = 'Created from %s using fastq_to_rdat.py\n' % (args.simplefile.name)
        elif filesin == 1:
            comments = 'Created from %s and %s using fastq_to_rdat.py\n' % (args.read1align.name, args.read2align.name)
        elif filesin == 2:
            comments = 'Created from %s and %s using fastq_to_rdat.py\n' % (args.read1fastq.name, args.read2fastq.name)
                
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

        output_string = 'RDAT file created: ' + filename
        sys.stdout.write( output_string+'\n' )
        f_log.write( output_string+'\n' )


def simple_to_rdat( args, sequence, filesin=2 ):
        print 'Starting analysis from simple file at ' + timeStamp()
        f_log.write( 'Starting analysis from simple file at ' + timeStamp() + '\n\n' )
        
        # Get length and seqpos
        WTlen = len(sequence)
        seqpos = arange(1, WTlen+1) + args.offset
        
        # Get mutation indices from simple data
        print 'Reading file: '+str(args.simplefile.name)

        # Set up arrays for 1D and 2D data
        count = 0
        fltnm = 0
        WTdata = np.zeros((1, WTlen))
        data2d = np.zeros(([len(seqpos),len(seqpos)]))
        mut_idxs   = []
        mut_counts = []
        start_pos_counts = np.zeros((1,WTlen))  
        stop_reactivity = np.zeros((1,WTlen))
        for line in args.simplefile.readlines():                
                # Read from file
                fields = line.split('\t')
                if len( fields ) < 3: continue
                start_pos = int(fields[0])
                simple_string = fields[2][:-1]  # Remove endline
                assert( len(simple_string) == (int( fields[1]) - int(fields[0]) + 1 ) )
                simple_line = [ int(x) for x in simple_string ]
                count += 1                      # Record total sequences
                if count % 50000 == 0: print 'Reading line number ', count

                # Simple-to-RDAT stop and mutation assignment
                mut_counts = 0
                start_pos_to_use = start_pos
                # record reactivity associated with reverse transcriptase stopping (probably should look at a 2D version of this, like in MOHCA-seq. -- rhiju)
                mod_pos = start_pos_to_use - 1                  # "start pos" actually records where reverse transcriptase stopped.
                # note further offset: mod_pos=0 means *no modifications*; mod_pos = 1 means modification at position 1
                # if we swtich to to 2D readout, maybe should decrement (add -1). Then no mod would actually *wrap* to WTlen (?)
                start_pos_counts[ 0, mod_pos ] += 1
                # record mutation indices
                if ( start_pos > args.start_pos_cutoff ): continue
                mut_idx = array( [ (idx + start_pos_to_use - 1) for idx, val in enumerate(simple_line) if val == 1])
                mut_counts = len( mut_idx )
                if ( mut_counts > args.num_hits_cutoff ): continue

                # Build data arrays
                if mut_counts > 0:
                        WTdata[ 0, mut_idx ] += 1               # Build 1D profile
                for mutpos in mut_idx:
                        data2d[mutpos, mut_idx] += 1            # Build 2D dataset

                fltnm += 1                                      # count sequences that pass filters above

        f_log.write( '\nTotal number of sequences: ' + str(count) )
        f_log.write( '\nNumber of sequences with fewer than 10 mismatches to WT (used to get 2D data): ' + str(fltnm) + '\n\n' )

        # Get reactivities, using Fi/[F0 + F1 ... + Fi] expression.
        sum_counts = start_pos_counts[ 0, 0 ]
        for (idx,counts) in enumerate(start_pos_counts[0,1:]):
                sum_counts += counts
                stop_reactivity[0, idx ] = ( float(counts)/sum_counts )
        filename = currdir + '/' + args.outprefix + '.stop_reactivity.rdat'
        output_rdat( filename, args, sequence, [], seqpos, stop_reactivity, np.zeros((0,WTlen)), f_log, filesin )

        # normalizes to total number of reads -- should be a true 'modification fraction'
        WTdata_reactivity = WTdata/fltnm
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
        output_rdat( filename, args, sequence, row_indices, seqpos, WTdata_reactivity, data2d_reactivity, f_log, filesin )

        filename = currdir + '/' + args.outprefix + '.rdat'
        output_rdat( filename, args, sequence, row_indices, seqpos, WTdata_norm, data2d_norm, f_log, filesin )

        f_log.write( '\nFinished analysis from simple file at ' + timeStamp() )
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

        # See if user input .simple data file, aligned sequence text files, or FASTQ files and run appropriate commands
        if args.simplefile is not None:
                f_log.write( 'Simple data file provided: ' + args.simplefile.name + '\nGenerating 2D data from simple data file.\n\n' )
                filesin = 0
                simple_to_rdat( args, sequence, filesin )

        elif args.read1align is not None and args.read2align is not None:
                f_log.write( 'Aligned sequences provided:\t(Read1) ' + args.read1align.name + '\t(Read 2) ' + args.read2align.name + '\nGenerating simple data file from aligned sequence files.\n\n')
                filesin = 1
                simple_path = aligned_to_simple( args )
                args.simplefile = open(simple_path, 'r')
                simple_to_rdat( args, sequence, filesin )

        elif args.read1fastq is not None and args.read2fastq is not None:
                f_log.write( 'FASTQs provided:\t(Read 1) ' + args.read1fastq.name + '\t(Read 2) ' + args.read2fastq.name + '\nGenerating simple data file from FASTQs.\n\n')
                filesin = 2
                simple_path = fastq_to_simple( args )
                args.simplefile = open(simple_path, 'r')
                simple_to_rdat( args, sequence, filesin )

        else:
                f_log.write( 'No data files provided.\n' )
                exit()



