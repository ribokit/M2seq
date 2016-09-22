#!/usr/bin/python

# Convert culled events text file to simplified format with "0" indicating no mutation and "1" indicating mutation
# Read start and end are 1-based indices
# Mate pairs that do not overlap are excluded

# Usage: python simplify_mutations.py <culled_events_file.txt>
# output: <culled_events_file.txt>.simple

import sys, os

fileIn = open(sys.argv[1], "rU")
lines = fileIn.readlines()
fileIn.close()

fileOut = open(sys.argv[1]+".simple","w")


class Read:
    def __init__(self):
        self.start = 0
        self.end =0
        self.simpleEvents = []
        self.inRead = False
        self.afterRead = False

lineCount = 0
for line in lines:
    lineCount += 1
    #if lineCount > 200:
    #    break
    #print line.strip()
    splitLine = line.strip().split()
    events = splitLine[2]
    start = int(splitLine[0])
    end = int(splitLine[1])
    r1 = Read()
    r1.start = int(start)
    r1.end = int(end)

    # identify reads with non-overlapping mate pairs (exclude these reads for now)
    noOverlap = False
    left = 0
    for i in range(0, len(events)):
        if events[i] not in ["~", "s"]:
            left = i
            break
    for i in range(left, len(events)):
        nuc = events[i]
        if nuc in ["~","s"]:
            try:
                if events[i+1] in ["|"]:
                   noOverlap = True
                   break
            except IndexError:
                pass

    if noOverlap == False:
        # update read start and end with clipped and uncovered nucs removed
        startIndex = 0
        for i in range(len(events)):
            if events[i] not in ["~","s"]:
                break
            else:
                r1.start += 1
                startIndex += 1
        endIndex = len(events)-1
        for i in range(len(events))[::-1]:
            if events[i] not in ["~","s"]:
                break
            else:
                r1.end -= 1
                endIndex -= 1

        for i in range(startIndex,endIndex+1):
            nuc = events[i]
            if nuc in ["|","~"]:
                r1.simpleEvents.append("0")
            else:
                r1.simpleEvents.append("1")
 
    if len(r1.simpleEvents) > 0:
        fileOut.write("%i\t%i\t%s\n"%(r1.start, r1.end, "".join(r1.simpleEvents)))        
    for i in range(len(r1.simpleEvents)-1):
        if r1.simpleEvents[i]=="1" and r1.simpleEvents[i+1]=="1":
            print line.strip().split()[2]
            print "".join(r1.simpleEvents)+"\n"
            break

