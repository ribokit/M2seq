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
#       (note: initially used ShapeMapper for pre-alignment and read trimming by Phred score, but
#       now uses nwalign without trimming; results are similar with and without)
# 3. Calculates 1D reactivity profiles from the demultiplexed FASTQs using ShapeMapper
#       (note: may write in-house scripts for 1D reactivity calculations in future)
#
# Each analysis is performed in a separate folder.
#
# Clarence Cheng, 2015
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
    f_log.write( '\nFinished demultiplexing at: ' + timeStamp() )

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
    f_log.write( line + '\n' )
    print cols[0]+'\t'+cols[1]


######################## Generate RDAT of 2D data using fastq_to_rdat.py ########################
make_dir( currdir + '/2_MaP2D' )

os.chdir( currdir )

print 'Starting MaP2D analysis'
f_log.write( '\nStarting MaP2D analysis at: ' + timeStamp() + '\n' )

# Move demultiplexed fastq files to MaP2D folder
old_fastq_names = [os.path.basename(args.read1fastq.name), os.path.basename(args.read2fastq.name)]
os.system('cp %s 2_MaP2D' % (args.sequencefile.name) )
for primer_tag in primer_tags:

    new_fastq_names = [ primer_tag+'_S1_L001_R1_001.fastq', primer_tag+'_S1_L001_R2_001.fastq' ]
    for (old_fastq_name,new_fastq_name) in zip(old_fastq_names,new_fastq_names):
        os.system('mv 1_Demultiplex/'+barcode_sequences[primer_tag]+'/'+old_fastq_name+' 1_Demultiplex/'+barcode_sequences[primer_tag]+'/'+new_fastq_name)
        os.system('mkdir -p 2_MaP2D/%s' % (primer_tag) )

    os.chdir( currdir + '/2_MaP2D/' + primer_tag )
    print 'Starting '+primer_tag
    fastq_names = [ '../../1_Demultiplex/'+barcode_sequences[primer_tag]+'/'+name for name in new_fastq_names]
    commandline = 'fastq_to_rdat.py ../' + args.sequencefile.name + ' --read1fastq ' + fastq_names[0] + ' --read2fastq ' + fastq_names[1] + ' --name ' + args.name.lstrip() + ' --offset ' + str(args.offset) + ' --outprefix ' + primer_tag
    print 'Command: '+commandline
    f_log.write( 'Command: ' + commandline + '\n' )
    os.system( commandline )
    os.chdir( currdir )

f_log.write( '\nFinished MaP2D analysis at: ' + timeStamp() )


os.chdir( currdir )


######################## Calculate 1D reactivities using ShapeMapper.py ########################
if args.config is not None:
    make_dir( currdir + '/3_ShapeMapper' )

    os.chdir( currdir )

    # Move demultiplexed fastq files to ShapeMapper folder
    for primer_tag in primer_tags:
        new_fastq_names = [ primer_tag+'_S1_L001_R1_001.fastq', primer_tag+'_S1_L001_R2_001.fastq' ]
        new_fastq_names_a = [ '1_Demultiplex/'+barcode_sequences[primer_tag]+'/'+name for name in new_fastq_names]
        new_fastq_names_b = [ '3_ShapeMapper/'+name for name in new_fastq_names]
        command1 = 'mv '+new_fastq_names_a[0]+' '+new_fastq_names_b[0]
        command2 = 'mv '+new_fastq_names_a[1]+' '+new_fastq_names_b[1]
        print command1
        print command2
        os.system( command1 )
        os.system( command2 )

    f_log.write( '\nStarting ShapeMapper analysis at: ' + timeStamp() )
    os.system('cp ' + args.sequencefile.name + ' ' + currdir + '/3_ShapeMapper/' + args.name + '.fa')
    os.system('cp ' + args.config.name + ' ' + currdir + '/3_ShapeMapper/' + args.name + '.cfg')
    os.chdir( currdir + '/3_ShapeMapper' )
    os.system('ShapeMapper.py ' + args.name + '.cfg')
    f_log.write( '\nFinished ShapeMapper analysis at: ' + timeStamp() )

    os.chdir( currdir )

    # Move demultiplexed fastq files back to Demultiplex folder
    for primer_tag in primer_tags:
        new_fastq_names = [ primer_tag+'_S1_L001_R1_001.fastq', primer_tag+'_S1_L001_R2_001.fastq' ]
        new_fastq_names_a = [ '1_Demultiplex/'+barcode_sequences[primer_tag]+'/'+name for name in new_fastq_names]
        new_fastq_names_b = [ '3_ShapeMapper/'+name for name in new_fastq_names]
        command1 = 'mv '+new_fastq_names_b[0]+' '+new_fastq_names_a[0]
        command2 = 'mv '+new_fastq_names_b[1]+' '+new_fastq_names_a[1]
        print command1
        print command2
        os.system( command1 )
        os.system( command2 )

f_log.close()
