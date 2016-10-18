#!/usr/bin/python

# Convert old-style mutation strings to simple strings, write to simple file

import os, sys
import argparse

parser = argparse.ArgumentParser( description='Input old-style mutation strings file')
parser.add_argument('mutsfile', type=argparse.FileType('r'), help='Name of old-style mutation strings file')
args = parser.parse_args()

simplename = args.mutsfile.name + '.simple'
f_simple = open(simplename,'w')

start_pos = []
count = 0

for line in args.mutsfile:
    fields = line.strip().split('\t')
    if len( fields ) < 3: continue

    start_pos   = int(fields[0])
    end_pos     = int(fields[1])
    mut_string  = fields[2]
    assert( len(mut_string) == (int( fields[1]) - int(fields[0]) + 1 ) )

    count += 1      # record total sequences
    if count % 50000 == 0: print 'Reading line number ', count


    # ignore non-overlapping reads
    temp = mut_string.strip('s~')
    non_overlap = 0
    for i in xrange( len(temp) ):
        if temp[i] == 's' or temp[i] == '~':
            if temp[i+1] == '|':
                non_overlap = 1
                break

    # adjust start and end positions
    if non_overlap == 0:
        for i in xrange( len(mut_string) ):
            if mut_string[i] == 's' or mut_string[i] == '~':
                start_pos += 1
            else:
                break

        for i in reversed( xrange( len(mut_string) ) ):
            if mut_string[i] == 's' or mut_string[i] == '~':
                end_pos -= 1
            else:
                break

        simple_line = mut_string.strip('s~').replace('|','0').replace('~','0').replace('A','1').replace('T','1').replace('G','1').replace('C','1').replace('-','1')

        f_simple.write(str(start_pos) + '\t' + str(end_pos) + '\t' + simple_line + '\n')

print '\nSimple file created: ' + f_simple.name
print '\nTotal number of sequences: ' + str(count)

f_simple.close()


