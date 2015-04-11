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
# To do:
#       -Parse specific barcodes provided in 'barcodes' input file, which may be a
#        subset of the 16 RTB primers used by Ann in the pilot experiment.
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

f_log.write( 'Starting Novobarcode demultiplexing at: ' + timeStamp() )
os.system('novobarcode -b ' + args.barcodes.name + ' -f ' + args.read1fastq.name + ' ' + args.read2fastq.name + ' -d 1_Demultiplex > 1_Demultiplex/novobarcode_log_Distance4.txt')
f_log.write( '\nFinished demultiplexing at: ' + timeStamp() )




######################## Generate RDAT of 2D data using fastq_to_rdat.py ########################
make_dir( currdir + '/2_MaP2D' )

# Move demultiplexed fastq files to MaP2D folder
os.chdir( currdir )

# Here and below, should parse 'barcodes' input file to get actual barcode sequences (may be different barcodes or a subset of barcodes)
# Read 1
# print 'mv 1_Demultiplex/ACCAGGCGCTGG/' + args.read1fastq.name + ' 2_MaP2D/RTB000_S1_L001_R1_001.fastq'
os.system('mv 1_Demultiplex/ACCAGGCGCTGG/' + args.read1fastq.name + ' 2_MaP2D/RTB000_S1_L001_R1_001.fastq')
os.system('mv 1_Demultiplex/GAGGCCTTGGCC/' + args.read1fastq.name + ' 2_MaP2D/RTB001_S1_L001_R1_001.fastq')
os.system('mv 1_Demultiplex/CTTTAAAATATA/' + args.read1fastq.name + ' 2_MaP2D/RTB002_S1_L001_R1_001.fastq')
os.system('mv 1_Demultiplex/TGACTTGCACAT/' + args.read1fastq.name + ' 2_MaP2D/RTB003_S1_L001_R1_001.fastq')
os.system('mv 1_Demultiplex/TGCGCCATTGCT/' + args.read1fastq.name + ' 2_MaP2D/RTB004_S1_L001_R1_001.fastq')
os.system('mv 1_Demultiplex/ACAAAATGGTGG/' + args.read1fastq.name + ' 2_MaP2D/RTB005_S1_L001_R1_001.fastq')
os.system('mv 1_Demultiplex/CTGCGTGCAAAC/' + args.read1fastq.name + ' 2_MaP2D/RTB006_S1_L001_R1_001.fastq')
os.system('mv 1_Demultiplex/GATTTGCACCTA/' + args.read1fastq.name + ' 2_MaP2D/RTB007_S1_L001_R1_001.fastq')
os.system('mv 1_Demultiplex/GGTATATGTACA/' + args.read1fastq.name + ' 2_MaP2D/RTB008_S1_L001_R1_001.fastq')
os.system('mv 1_Demultiplex/CCCGCGCTGGGT/' + args.read1fastq.name + ' 2_MaP2D/RTB009_S1_L001_R1_001.fastq')
os.system('mv 1_Demultiplex/ATGCATGCACAG/' + args.read1fastq.name + ' 2_MaP2D/RTB010_S1_L001_R1_001.fastq')
os.system('mv 1_Demultiplex/TAATGCAACTTC/' + args.read1fastq.name + ' 2_MaP2D/RTB011_S1_L001_R1_001.fastq')
os.system('mv 1_Demultiplex/GCAAATGTGCTA/' + args.read1fastq.name + ' 2_MaP2D/RTB012_S1_L001_R1_001.fastq')
os.system('mv 1_Demultiplex/TGGCGAACATGG/' + args.read1fastq.name + ' 2_MaP2D/RTB013_S1_L001_R1_001.fastq')
os.system('mv 1_Demultiplex/CTTTCCCACACT/' + args.read1fastq.name + ' 2_MaP2D/RTB014_S1_L001_R1_001.fastq')
os.system('mv 1_Demultiplex/AACGTGTGTGAC/' + args.read1fastq.name + ' 2_MaP2D/RTB015_S1_L001_R1_001.fastq')
# Read 2
os.system('mv 1_Demultiplex/ACCAGGCGCTGG/' + args.read2fastq.name + ' 2_MaP2D/RTB000_S1_L001_R2_001.fastq')
os.system('mv 1_Demultiplex/GAGGCCTTGGCC/' + args.read2fastq.name + ' 2_MaP2D/RTB001_S1_L001_R2_001.fastq')
os.system('mv 1_Demultiplex/CTTTAAAATATA/' + args.read2fastq.name + ' 2_MaP2D/RTB002_S1_L001_R2_001.fastq')
os.system('mv 1_Demultiplex/TGACTTGCACAT/' + args.read2fastq.name + ' 2_MaP2D/RTB003_S1_L001_R2_001.fastq')
os.system('mv 1_Demultiplex/TGCGCCATTGCT/' + args.read2fastq.name + ' 2_MaP2D/RTB004_S1_L001_R2_001.fastq')
os.system('mv 1_Demultiplex/ACAAAATGGTGG/' + args.read2fastq.name + ' 2_MaP2D/RTB005_S1_L001_R2_001.fastq')
os.system('mv 1_Demultiplex/CTGCGTGCAAAC/' + args.read2fastq.name + ' 2_MaP2D/RTB006_S1_L001_R2_001.fastq')
os.system('mv 1_Demultiplex/GATTTGCACCTA/' + args.read2fastq.name + ' 2_MaP2D/RTB007_S1_L001_R2_001.fastq')
os.system('mv 1_Demultiplex/GGTATATGTACA/' + args.read2fastq.name + ' 2_MaP2D/RTB008_S1_L001_R2_001.fastq')
os.system('mv 1_Demultiplex/CCCGCGCTGGGT/' + args.read2fastq.name + ' 2_MaP2D/RTB009_S1_L001_R2_001.fastq')
os.system('mv 1_Demultiplex/ATGCATGCACAG/' + args.read2fastq.name + ' 2_MaP2D/RTB010_S1_L001_R2_001.fastq')
os.system('mv 1_Demultiplex/TAATGCAACTTC/' + args.read2fastq.name + ' 2_MaP2D/RTB011_S1_L001_R2_001.fastq')
os.system('mv 1_Demultiplex/GCAAATGTGCTA/' + args.read2fastq.name + ' 2_MaP2D/RTB012_S1_L001_R2_001.fastq')
os.system('mv 1_Demultiplex/TGGCGAACATGG/' + args.read2fastq.name + ' 2_MaP2D/RTB013_S1_L001_R2_001.fastq')
os.system('mv 1_Demultiplex/CTTTCCCACACT/' + args.read2fastq.name + ' 2_MaP2D/RTB014_S1_L001_R2_001.fastq')
os.system('mv 1_Demultiplex/AACGTGTGTGAC/' + args.read2fastq.name + ' 2_MaP2D/RTB015_S1_L001_R2_001.fastq')


os.system('cp ' + args.sequencefile.name + ' ' + currdir + '/2_MaP2D/' + args.name + '.fa')

for i in range(0, 16):
    if i < 10:
        make_dir( currdir + '/2_MaP2D/RTB00' + str(i) )
    else:
        make_dir( currdir + '/2_MaP2D/RTB0' + str(i) )

os.chdir( currdir + '/2_MaP2D' )
f_log.write( '\nStarting MaP2D analysis at: ' + timeStamp() )
for i in range(0, 16):
    if i < 10:
        os.system('fastq_to_rdat.py ' + args.name + '.fa' + ' --read1fastq RTB00' + str(i) + '_S1_L001_R1_001.fastq --read2fastq RTB00' + str(i) + '_S1_L001_R2_001.fastq --name ' + args.name + ' --offset ' + str(args.offset) + ' --outprefix RTB00' + str(i) + '/out')
    else:
        os.system('fastq_to_rdat.py ' + args.name + '.fa' + ' --read1fastq RTB0' + str(i) + '_S1_L001_R1_001.fastq --read2fastq RTB0' + str(i) + '_S1_L001_R2_001.fastq --name ' + args.name + ' --offset ' + str(args.offset) + ' --outprefix RTB0' + str(i) + '/out')
f_log.write( '\nFinished MaP2D analysis at: ' + timeStamp() )

os.chdir( currdir )


# Move demultiplexed fastq files back to demultiplex folder
# Read 1
os.system('mv 2_MaP2D/RTB000_S1_L001_R1_001.fastq 1_Demultiplex/ACCAGGCGCTGG/RTB000_S1_L001_R1_001.fastq')
os.system('mv 2_MaP2D/RTB001_S1_L001_R1_001.fastq 1_Demultiplex/GAGGCCTTGGCC/RTB001_S1_L001_R1_001.fastq')
os.system('mv 2_MaP2D/RTB002_S1_L001_R1_001.fastq 1_Demultiplex/CTTTAAAATATA/RTB002_S1_L001_R1_001.fastq')
os.system('mv 2_MaP2D/RTB003_S1_L001_R1_001.fastq 1_Demultiplex/TGACTTGCACAT/RTB003_S1_L001_R1_001.fastq')
os.system('mv 2_MaP2D/RTB004_S1_L001_R1_001.fastq 1_Demultiplex/TGCGCCATTGCT/RTB004_S1_L001_R1_001.fastq')
os.system('mv 2_MaP2D/RTB005_S1_L001_R1_001.fastq 1_Demultiplex/ACAAAATGGTGG/RTB005_S1_L001_R1_001.fastq')
os.system('mv 2_MaP2D/RTB006_S1_L001_R1_001.fastq 1_Demultiplex/CTGCGTGCAAAC/RTB006_S1_L001_R1_001.fastq')
os.system('mv 2_MaP2D/RTB007_S1_L001_R1_001.fastq 1_Demultiplex/GATTTGCACCTA/RTB007_S1_L001_R1_001.fastq')
os.system('mv 2_MaP2D/RTB008_S1_L001_R1_001.fastq 1_Demultiplex/GGTATATGTACA/RTB008_S1_L001_R1_001.fastq')
os.system('mv 2_MaP2D/RTB009_S1_L001_R1_001.fastq 1_Demultiplex/CCCGCGCTGGGT/RTB009_S1_L001_R1_001.fastq')
os.system('mv 2_MaP2D/RTB010_S1_L001_R1_001.fastq 1_Demultiplex/ATGCATGCACAG/RTB010_S1_L001_R1_001.fastq')
os.system('mv 2_MaP2D/RTB011_S1_L001_R1_001.fastq 1_Demultiplex/TAATGCAACTTC/RTB011_S1_L001_R1_001.fastq')
os.system('mv 2_MaP2D/RTB012_S1_L001_R1_001.fastq 1_Demultiplex/GCAAATGTGCTA/RTB012_S1_L001_R1_001.fastq')
os.system('mv 2_MaP2D/RTB013_S1_L001_R1_001.fastq 1_Demultiplex/TGGCGAACATGG/RTB013_S1_L001_R1_001.fastq')
os.system('mv 2_MaP2D/RTB014_S1_L001_R1_001.fastq 1_Demultiplex/CTTTCCCACACT/RTB014_S1_L001_R1_001.fastq')
os.system('mv 2_MaP2D/RTB015_S1_L001_R1_001.fastq 1_Demultiplex/AACGTGTGTGAC/RTB015_S1_L001_R1_001.fastq')
# Read 2
os.system('mv 2_MaP2D/RTB000_S1_L001_R2_001.fastq 1_Demultiplex/ACCAGGCGCTGG/RTB000_S1_L001_R2_001.fastq')
os.system('mv 2_MaP2D/RTB001_S1_L001_R2_001.fastq 1_Demultiplex/GAGGCCTTGGCC/RTB001_S1_L001_R2_001.fastq')
os.system('mv 2_MaP2D/RTB002_S1_L001_R2_001.fastq 1_Demultiplex/CTTTAAAATATA/RTB002_S1_L001_R2_001.fastq')
os.system('mv 2_MaP2D/RTB003_S1_L001_R2_001.fastq 1_Demultiplex/TGACTTGCACAT/RTB003_S1_L001_R2_001.fastq')
os.system('mv 2_MaP2D/RTB004_S1_L001_R2_001.fastq 1_Demultiplex/TGCGCCATTGCT/RTB004_S1_L001_R2_001.fastq')
os.system('mv 2_MaP2D/RTB005_S1_L001_R2_001.fastq 1_Demultiplex/ACAAAATGGTGG/RTB005_S1_L001_R2_001.fastq')
os.system('mv 2_MaP2D/RTB006_S1_L001_R2_001.fastq 1_Demultiplex/CTGCGTGCAAAC/RTB006_S1_L001_R2_001.fastq')
os.system('mv 2_MaP2D/RTB007_S1_L001_R2_001.fastq 1_Demultiplex/GATTTGCACCTA/RTB007_S1_L001_R2_001.fastq')
os.system('mv 2_MaP2D/RTB008_S1_L001_R2_001.fastq 1_Demultiplex/GGTATATGTACA/RTB008_S1_L001_R2_001.fastq')
os.system('mv 2_MaP2D/RTB009_S1_L001_R2_001.fastq 1_Demultiplex/CCCGCGCTGGGT/RTB009_S1_L001_R2_001.fastq')
os.system('mv 2_MaP2D/RTB010_S1_L001_R2_001.fastq 1_Demultiplex/ATGCATGCACAG/RTB010_S1_L001_R2_001.fastq')
os.system('mv 2_MaP2D/RTB011_S1_L001_R2_001.fastq 1_Demultiplex/TAATGCAACTTC/RTB011_S1_L001_R2_001.fastq')
os.system('mv 2_MaP2D/RTB012_S1_L001_R2_001.fastq 1_Demultiplex/GCAAATGTGCTA/RTB012_S1_L001_R2_001.fastq')
os.system('mv 2_MaP2D/RTB013_S1_L001_R2_001.fastq 1_Demultiplex/TGGCGAACATGG/RTB013_S1_L001_R2_001.fastq')
os.system('mv 2_MaP2D/RTB014_S1_L001_R2_001.fastq 1_Demultiplex/CTTTCCCACACT/RTB014_S1_L001_R2_001.fastq')
os.system('mv 2_MaP2D/RTB015_S1_L001_R2_001.fastq 1_Demultiplex/AACGTGTGTGAC/RTB015_S1_L001_R2_001.fastq')



######################## Calculate 1D reactivities using ShapeMapper.py ########################
if args.config is not None:
    make_dir( currdir + '/3_ShapeMapper' )

    os.chdir( currdir )

    # Move demultiplexed fastq files to ShapeMapper folder
    # Read 1
    os.system('mv 1_Demultiplex/ACCAGGCGCTGG/RTB000_S1_L001_R1_001.fastq 3_ShapeMapper/RTB000_S1_L001_R1_001.fastq')
    os.system('mv 1_Demultiplex/GAGGCCTTGGCC/RTB001_S1_L001_R1_001.fastq 3_ShapeMapper/RTB001_S1_L001_R1_001.fastq')
    os.system('mv 1_Demultiplex/CTTTAAAATATA/RTB002_S1_L001_R1_001.fastq 3_ShapeMapper/RTB002_S1_L001_R1_001.fastq')
    os.system('mv 1_Demultiplex/TGACTTGCACAT/RTB003_S1_L001_R1_001.fastq 3_ShapeMapper/RTB003_S1_L001_R1_001.fastq')
    os.system('mv 1_Demultiplex/TGCGCCATTGCT/RTB004_S1_L001_R1_001.fastq 3_ShapeMapper/RTB004_S1_L001_R1_001.fastq')
    os.system('mv 1_Demultiplex/ACAAAATGGTGG/RTB005_S1_L001_R1_001.fastq 3_ShapeMapper/RTB005_S1_L001_R1_001.fastq')
    os.system('mv 1_Demultiplex/CTGCGTGCAAAC/RTB006_S1_L001_R1_001.fastq 3_ShapeMapper/RTB006_S1_L001_R1_001.fastq')
    os.system('mv 1_Demultiplex/GATTTGCACCTA/RTB007_S1_L001_R1_001.fastq 3_ShapeMapper/RTB007_S1_L001_R1_001.fastq')
    os.system('mv 1_Demultiplex/GGTATATGTACA/RTB008_S1_L001_R1_001.fastq 3_ShapeMapper/RTB008_S1_L001_R1_001.fastq')
    os.system('mv 1_Demultiplex/CCCGCGCTGGGT/RTB009_S1_L001_R1_001.fastq 3_ShapeMapper/RTB009_S1_L001_R1_001.fastq')
    os.system('mv 1_Demultiplex/ATGCATGCACAG/RTB010_S1_L001_R1_001.fastq 3_ShapeMapper/RTB010_S1_L001_R1_001.fastq')
    os.system('mv 1_Demultiplex/TAATGCAACTTC/RTB011_S1_L001_R1_001.fastq 3_ShapeMapper/RTB011_S1_L001_R1_001.fastq')
    os.system('mv 1_Demultiplex/GCAAATGTGCTA/RTB012_S1_L001_R1_001.fastq 3_ShapeMapper/RTB012_S1_L001_R1_001.fastq')
    os.system('mv 1_Demultiplex/TGGCGAACATGG/RTB013_S1_L001_R1_001.fastq 3_ShapeMapper/RTB013_S1_L001_R1_001.fastq')
    os.system('mv 1_Demultiplex/CTTTCCCACACT/RTB014_S1_L001_R1_001.fastq 3_ShapeMapper/RTB014_S1_L001_R1_001.fastq')
    os.system('mv 1_Demultiplex/AACGTGTGTGAC/RTB015_S1_L001_R1_001.fastq 3_ShapeMapper/RTB015_S1_L001_R1_001.fastq')
    # Read 2
    os.system('mv 1_Demultiplex/ACCAGGCGCTGG/RTB000_S1_L001_R2_001.fastq 3_ShapeMapper/RTB000_S1_L001_R2_001.fastq')
    os.system('mv 1_Demultiplex/GAGGCCTTGGCC/RTB001_S1_L001_R2_001.fastq 3_ShapeMapper/RTB001_S1_L001_R2_001.fastq')
    os.system('mv 1_Demultiplex/CTTTAAAATATA/RTB002_S1_L001_R2_001.fastq 3_ShapeMapper/RTB002_S1_L001_R2_001.fastq')
    os.system('mv 1_Demultiplex/TGACTTGCACAT/RTB003_S1_L001_R2_001.fastq 3_ShapeMapper/RTB003_S1_L001_R2_001.fastq')
    os.system('mv 1_Demultiplex/TGCGCCATTGCT/RTB004_S1_L001_R2_001.fastq 3_ShapeMapper/RTB004_S1_L001_R2_001.fastq')
    os.system('mv 1_Demultiplex/ACAAAATGGTGG/RTB005_S1_L001_R2_001.fastq 3_ShapeMapper/RTB005_S1_L001_R2_001.fastq')
    os.system('mv 1_Demultiplex/CTGCGTGCAAAC/RTB006_S1_L001_R2_001.fastq 3_ShapeMapper/RTB006_S1_L001_R2_001.fastq')
    os.system('mv 1_Demultiplex/GATTTGCACCTA/RTB007_S1_L001_R2_001.fastq 3_ShapeMapper/RTB007_S1_L001_R2_001.fastq')
    os.system('mv 1_Demultiplex/GGTATATGTACA/RTB008_S1_L001_R2_001.fastq 3_ShapeMapper/RTB008_S1_L001_R2_001.fastq')
    os.system('mv 1_Demultiplex/CCCGCGCTGGGT/RTB009_S1_L001_R2_001.fastq 3_ShapeMapper/RTB009_S1_L001_R2_001.fastq')
    os.system('mv 1_Demultiplex/ATGCATGCACAG/RTB010_S1_L001_R2_001.fastq 3_ShapeMapper/RTB010_S1_L001_R2_001.fastq')
    os.system('mv 1_Demultiplex/TAATGCAACTTC/RTB011_S1_L001_R2_001.fastq 3_ShapeMapper/RTB011_S1_L001_R2_001.fastq')
    os.system('mv 1_Demultiplex/GCAAATGTGCTA/RTB012_S1_L001_R2_001.fastq 3_ShapeMapper/RTB012_S1_L001_R2_001.fastq')
    os.system('mv 1_Demultiplex/TGGCGAACATGG/RTB013_S1_L001_R2_001.fastq 3_ShapeMapper/RTB013_S1_L001_R2_001.fastq')
    os.system('mv 1_Demultiplex/CTTTCCCACACT/RTB014_S1_L001_R2_001.fastq 3_ShapeMapper/RTB014_S1_L001_R2_001.fastq')
    os.system('mv 1_Demultiplex/AACGTGTGTGAC/RTB015_S1_L001_R2_001.fastq 3_ShapeMapper/RTB015_S1_L001_R2_001.fastq')


    f_log.write( '\nStarting ShapeMapper analysis at: ' + timeStamp() )
    # print 'cp ' + args.sequencefile.name + ' ' + currdir + '/3_ShapeMapper/' + args.name + '.fa'
    # print 'cp ' + args.config.name + ' ' + currdir + '/3_ShapeMapper/' + args.name + '.cfg'
    os.system('cp ' + args.sequencefile.name + ' ' + currdir + '/3_ShapeMapper/' + args.name + '.fa')
    os.system('cp ' + args.config.name + ' ' + currdir + '/3_ShapeMapper/' + args.name + '.cfg')
    os.chdir( currdir + '/3_ShapeMapper' )
    os.system('ShapeMapper.py ' + args.name + '.cfg')
    f_log.write( '\nFinished ShapeMapper analysis at: ' + timeStamp() )


    os.chdir( currdir )

    # Move demultiplexed fastq files back to Demultiplex folder
    # Read 1
    os.system('mv 3_ShapeMapper/RTB000_S1_L001_R1_001.fastq 1_Demultiplex/ACCAGGCGCTGG/RTB000_S1_L001_R1_001.fastq')
    os.system('mv 3_ShapeMapper/RTB001_S1_L001_R1_001.fastq 1_Demultiplex/GAGGCCTTGGCC/RTB001_S1_L001_R1_001.fastq')
    os.system('mv 3_ShapeMapper/RTB002_S1_L001_R1_001.fastq 1_Demultiplex/CTTTAAAATATA/RTB002_S1_L001_R1_001.fastq')
    os.system('mv 3_ShapeMapper/RTB003_S1_L001_R1_001.fastq 1_Demultiplex/TGACTTGCACAT/RTB003_S1_L001_R1_001.fastq')
    os.system('mv 3_ShapeMapper/RTB004_S1_L001_R1_001.fastq 1_Demultiplex/TGCGCCATTGCT/RTB004_S1_L001_R1_001.fastq')
    os.system('mv 3_ShapeMapper/RTB005_S1_L001_R1_001.fastq 1_Demultiplex/ACAAAATGGTGG/RTB005_S1_L001_R1_001.fastq')
    os.system('mv 3_ShapeMapper/RTB006_S1_L001_R1_001.fastq 1_Demultiplex/CTGCGTGCAAAC/RTB006_S1_L001_R1_001.fastq')
    os.system('mv 3_ShapeMapper/RTB007_S1_L001_R1_001.fastq 1_Demultiplex/GATTTGCACCTA/RTB007_S1_L001_R1_001.fastq')
    os.system('mv 3_ShapeMapper/RTB008_S1_L001_R1_001.fastq 1_Demultiplex/GGTATATGTACA/RTB008_S1_L001_R1_001.fastq')
    os.system('mv 3_ShapeMapper/RTB009_S1_L001_R1_001.fastq 1_Demultiplex/CCCGCGCTGGGT/RTB009_S1_L001_R1_001.fastq')
    os.system('mv 3_ShapeMapper/RTB010_S1_L001_R1_001.fastq 1_Demultiplex/ATGCATGCACAG/RTB010_S1_L001_R1_001.fastq')
    os.system('mv 3_ShapeMapper/RTB011_S1_L001_R1_001.fastq 1_Demultiplex/TAATGCAACTTC/RTB011_S1_L001_R1_001.fastq')
    os.system('mv 3_ShapeMapper/RTB012_S1_L001_R1_001.fastq 1_Demultiplex/GCAAATGTGCTA/RTB012_S1_L001_R1_001.fastq')
    os.system('mv 3_ShapeMapper/RTB013_S1_L001_R1_001.fastq 1_Demultiplex/TGGCGAACATGG/RTB013_S1_L001_R1_001.fastq')
    os.system('mv 3_ShapeMapper/RTB014_S1_L001_R1_001.fastq 1_Demultiplex/CTTTCCCACACT/RTB014_S1_L001_R1_001.fastq')
    os.system('mv 3_ShapeMapper/RTB015_S1_L001_R1_001.fastq 1_Demultiplex/AACGTGTGTGAC/RTB015_S1_L001_R1_001.fastq')
    # Read 2
    os.system('mv 3_ShapeMapper/RTB000_S1_L001_R2_001.fastq 1_Demultiplex/ACCAGGCGCTGG/RTB000_S1_L001_R2_001.fastq')
    os.system('mv 3_ShapeMapper/RTB001_S1_L001_R2_001.fastq 1_Demultiplex/GAGGCCTTGGCC/RTB001_S1_L001_R2_001.fastq')
    os.system('mv 3_ShapeMapper/RTB002_S1_L001_R2_001.fastq 1_Demultiplex/CTTTAAAATATA/RTB002_S1_L001_R2_001.fastq')
    os.system('mv 3_ShapeMapper/RTB003_S1_L001_R2_001.fastq 1_Demultiplex/TGACTTGCACAT/RTB003_S1_L001_R2_001.fastq')
    os.system('mv 3_ShapeMapper/RTB004_S1_L001_R2_001.fastq 1_Demultiplex/TGCGCCATTGCT/RTB004_S1_L001_R2_001.fastq')
    os.system('mv 3_ShapeMapper/RTB005_S1_L001_R2_001.fastq 1_Demultiplex/ACAAAATGGTGG/RTB005_S1_L001_R2_001.fastq')
    os.system('mv 3_ShapeMapper/RTB006_S1_L001_R2_001.fastq 1_Demultiplex/CTGCGTGCAAAC/RTB006_S1_L001_R2_001.fastq')
    os.system('mv 3_ShapeMapper/RTB007_S1_L001_R2_001.fastq 1_Demultiplex/GATTTGCACCTA/RTB007_S1_L001_R2_001.fastq')
    os.system('mv 3_ShapeMapper/RTB008_S1_L001_R2_001.fastq 1_Demultiplex/GGTATATGTACA/RTB008_S1_L001_R2_001.fastq')
    os.system('mv 3_ShapeMapper/RTB009_S1_L001_R2_001.fastq 1_Demultiplex/CCCGCGCTGGGT/RTB009_S1_L001_R2_001.fastq')
    os.system('mv 3_ShapeMapper/RTB010_S1_L001_R2_001.fastq 1_Demultiplex/ATGCATGCACAG/RTB010_S1_L001_R2_001.fastq')
    os.system('mv 3_ShapeMapper/RTB011_S1_L001_R2_001.fastq 1_Demultiplex/TAATGCAACTTC/RTB011_S1_L001_R2_001.fastq')
    os.system('mv 3_ShapeMapper/RTB012_S1_L001_R2_001.fastq 1_Demultiplex/GCAAATGTGCTA/RTB012_S1_L001_R2_001.fastq')
    os.system('mv 3_ShapeMapper/RTB013_S1_L001_R2_001.fastq 1_Demultiplex/TGGCGAACATGG/RTB013_S1_L001_R2_001.fastq')
    os.system('mv 3_ShapeMapper/RTB014_S1_L001_R2_001.fastq 1_Demultiplex/CTTTCCCACACT/RTB014_S1_L001_R2_001.fastq')
    os.system('mv 3_ShapeMapper/RTB015_S1_L001_R2_001.fastq 1_Demultiplex/AACGTGTGTGAC/RTB015_S1_L001_R2_001.fastq')


f_log.close()
