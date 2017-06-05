#!/usr/bin/env python
###############################################################################
# Analysis of 2D signal in mutational profiling sequencing data (MaP2D)
###############################################################################
#
# MaP2D
#
# This script takes raw FASTQ files from a sequencing run and performs the following:
#
# 1. Demultiplexes the raw FASTQs using user-provided barcodes using Novobarcode
# 2. Uses the demultiplexed FASTQs and the WT sequence to generate 2D mutational profiling data
#       Use ShapeMapper for read alignment to reference sequence and generation of mutation strings
# 3. Calculates 2D datasets and outputs RDAT files using simple_to_rdat.py
#
# Clarence Cheng, 2015-2016
#

import os, sys, time
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('sequencefile', type=argparse.FileType('r'))
parser.add_argument('barcodes', type=argparse.FileType('r'))
parser.add_argument('read1fastq', type=argparse.FileType('r'))
parser.add_argument('read2fastq', type=argparse.FileType('r'))
parser.add_argument('--config', type=argparse.FileType('r'))
parser.add_argument('--name', type=str, default='PLACEHOLDER')
parser.add_argument('--offset', type=int, default=0)
parser.add_argument('--outprefix', type=str, default='out')

args = parser.parse_args()

sequencefile_lines = args.sequencefile.readlines()
sequence = sequencefile_lines[1].strip().upper()

if args.name == 'PLACEHOLDER':
    args.name = sequencefile_lines[0].strip()[1:]
    print args.name

currdir = os.getcwd()


def timeStamp():
    t = time.localtime()
    month = t.tm_mon
    day = t.tm_mday
    hour = t.tm_hour
    minute = t.tm_min
    second = t.tm_sec
    return '%i:%i:%i, %i/%i'%(hour, minute, second, month, day)


def make_dir(path):
    try:
        os.mkdir(path)
    except OSError:
        if not os.path.isdir(path):
            raise


######################## Demultiplex using Novobarcode ########################
f_log = open(currdir + '/' + 'AnalysisLog.txt', 'w')

make_dir( currdir + '/1_Demultiplex' )
if not os.path.exists( '1_Demultiplex/novobarcode_log_Distance4.txt' ):
    f_log.write( 'Starting Novobarcode demultiplexing at: ' + timeStamp() )
    print 'Starting Novobarcode demultiplexing'
    os.system('novobarcode -b ' + args.barcodes.name + ' -f ' + args.read1fastq.name + ' ' + args.read2fastq.name + ' -d 1_Demultiplex > 1_Demultiplex/novobarcode_log_Distance4.txt')
    f_log.write( '\nFinished demultiplexing at: ' + timeStamp() + '\n')

print 'Reading barcodes from: '+args.barcodes.name
f_log.write( 'Primers:\n' )
lines = open( args.barcodes.name ).readlines()
primer_tags = []
barcode_sequences = {}
for line in lines:
    if len( line ) < 2: continue
    col = line.rstrip( '\n' )
    cols = col.split( '\t' )
    if len( cols ) != 2: continue
    if len( cols[1] ) < 2: continue
    primer_tags.append( cols[0] )
    barcode_sequences[ cols[0] ] = cols[1]
    f_log.write( line )
    print cols[0]+'\t'+cols[1]


######################## Analyze using ShapeMapper.py ########################
# NOTE: .cfg config file input to ShapeMapper should have the following settings:
#  Under trimReads options, minPhred = 0
#  Under countMutations options, minPhredtoCount = 20
#  Under countMutations options, makeOldMutationStrings = on (write mutation strings in format that can be converted to binary .simple files)

if args.config is not None:
    make_dir( currdir + '/2_ShapeMapper' )

    os.chdir( currdir )

    # Move demultiplexed fastq files to ShapeMapper folder
    old_fastq_names = [os.path.basename(args.read1fastq.name), os.path.basename(args.read2fastq.name)]
    for primer_tag in primer_tags:
        new_fastq_names = [ primer_tag+'_S1_L001_R1_001.fastq', primer_tag+'_S1_L001_R2_001.fastq' ]
        # print new_fastq_names
        for (old_fastq_name,new_fastq_name) in zip(old_fastq_names,new_fastq_names):
            os.system('mv 1_Demultiplex/'+barcode_sequences[primer_tag]+'/'+old_fastq_name+' 2_ShapeMapper/'+new_fastq_name)
            print 'mv 1_Demultiplex/'+barcode_sequences[primer_tag]+'/'+old_fastq_name+' 2_ShapeMapper/'+new_fastq_name

    # Run ShapeMapper
    f_log.write( '\nStarting ShapeMapper analysis at: ' + timeStamp() )
    print 'Starting ShapeMapper analysis'
    print 'cp ' + args.sequencefile.name + ' "' + currdir + '"/2_ShapeMapper/' + args.name + '.fa'
    print 'cp ' + args.config.name + ' "' + currdir + '"/2_ShapeMapper/' + args.name + '.cfg'
    os.system('cp ' + args.sequencefile.name + ' "' + currdir + '"/2_ShapeMapper/' + args.name + '.fa')
    os.system('cp ' + args.config.name + ' "' + currdir + '"/2_ShapeMapper/' + args.name + '.cfg')
    os.chdir( currdir + '/2_ShapeMapper' )
    command_ShapeMapper = 'ShapeMapper.py ' + args.name + '.cfg'
    f_log.write( '\nShapeMapper command: ' + command_ShapeMapper )
    os.system( command_ShapeMapper )
    f_log.write( '\nFinished ShapeMapper analysis at: ' + timeStamp() )

    # Generate simple files
    f_log.write( '\nGenerating simple files at: ' + timeStamp() )
    os.chdir( currdir + '/2_ShapeMapper/output/mutation_strings_oldstyle/' )
    for file in os.listdir(currdir + '/2_ShapeMapper/output/mutation_strings_oldstyle/'):
        if file.endswith('.txt'):
            os.system('muts_to_simple.py ' + file)
    f_log.write( '\nFinished generating simple files at: ' + timeStamp() )

    os.chdir( currdir )

    # Move demultiplexed fastq files back to Demultiplex folder
    for primer_tag in primer_tags:
        new_fastq_names = [ primer_tag+'_S1_L001_R1_001.fastq', primer_tag+'_S1_L001_R2_001.fastq' ]
        new_fastq_names_a = [ '1_Demultiplex/'+barcode_sequences[primer_tag]+'/'+name for name in new_fastq_names]
        new_fastq_names_b = [ '2_ShapeMapper/'+name for name in new_fastq_names]
        command1 = 'mv '+new_fastq_names_b[0]+' '+new_fastq_names_a[0]
        command2 = 'mv '+new_fastq_names_b[1]+' '+new_fastq_names_a[1]
        print command1
        print command2
        os.system( command1 )
        os.system( command2 )


######################## Generate RDAT of 2D data using simple_to_rdat.py ########################
if args.config is not None:
    make_dir( currdir + '/3_MaP2D' )
    make_dir( currdir + '/3_MaP2D/simple_files')
    os.chdir( currdir )

    print 'Starting MaP2D analysis'
    f_log.write( '\nStarting MaP2D analysis at: ' + timeStamp() + '\n' )

    # Move simple format files and counted mutations to MaP2D folder
    # os.chdir( currdir )
    for file in os.listdir(currdir + '/2_ShapeMapper/output/mutation_strings_oldstyle/'):
        if file.endswith('.simple'):
            os.system('mv "' + currdir + '"/2_ShapeMapper/output/mutation_strings_oldstyle/' + file + ' "' + currdir + '"/3_MaP2D/simple_files/')

    # Run simple_to_rdat.py
    os.system('cp ' + args.sequencefile.name + ' "' + currdir + '"/3_MaP2D/' + args.name + '.fa')
    os.chdir( currdir + '/3_MaP2D/simple_files')
    for file in os.listdir(currdir + '/3_MaP2D/simple_files'):
        if file.endswith('.simple'):
            command_simple2rdat = 'simple_to_rdat.py ' + '../' + args.sequencefile.name + ' --simplefile ' + file + ' --name ' + args.name + ' --offset ' + str(args.offset) + ' --outprefix ' + file.split('.')[0]
            # note: getting different .fa files for a single pair of FASTQs (multiple RNAs per sequencing run) is currently unsupported
            print command_simple2rdat
            f_log.write( '\nMaP2D command: ' + command_simple2rdat )
            os.system( command_simple2rdat )

    f_log.write( '\nFinished MaP2D analysis at: ' + timeStamp() + '\n' )

os.chdir( currdir )

f_log.close()

