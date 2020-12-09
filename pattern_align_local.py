#! /usr/bin/env python

"""
Given query sequence, subject sequence, and subject annotations,
align query to subject, then transfer annotations to query;
gapped positions in the query will be ignored.

Example:
q = AAAAGTTTT
s = AAAATTTTC
p = 012345679

Alignment:
q = AAAAGTTTT-
s = AAAA-TTTTC
p = 0123-45678

Output:
q = AAAAGTTTT
s = AAAA-TTTT
p = 0123-4567 ***

===============================================
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""

import os, sys, re, glob, argparse, subprocess

def align ( query, subject, pattern ): 

    # process pattern
    if type( pattern ) is not list:
        pattern = re.split( "\s+", pattern ) if " " in pattern else [k for k in pattern]
    pgap = "-" * max( [len(item) for item in pattern] )
    if len( subject ) != len( pattern ):
        sys.exit( "subject and pattern are not the same length" )

    # write the query and subject sequences to memory
    with open( "/tmp/q.temp", "w" ) as fh:
        print >>fh, ">q.temp"
        print >>fh, query
    with open( "/tmp/s.temp", "w" ) as fh:
        print >>fh, ">s.temp"
        print >>fh, subject

    # execute blast alignment (bl2seq equivalent)
    items = [
        "blastp",
        "-query /tmp/q.temp",
        "-subject /tmp/s.temp",
        "-outfmt \"6 qlen slen qstart qend sstart send qseq sseq\"",
        ]
    cmd = subprocess.Popen( " ".join( items ), shell=True, stdout=subprocess.PIPE )
    results = cmd.stdout.readline().strip().split( "\t" )

    # format blast result
    qlen, slen, qstart, qend, sstart, send = [int( k ) for k in results[0:-2]]
    qseq, sseq = results[-2:]
    alen = len( qseq )

    # use alignment result to align pattern
    pseq = []
    qindex = qstart - 1
    sindex = sstart - 1
    for aindex in range( alen ):
        # no query, skip subject/pattern info
        if qseq[aindex] == "-":
            qindex += 0
            sindex += 1
        # no subject, gap
        elif sseq[aindex] == "-":
            pseq.append( pgap )
            qindex += 1
            sindex += 0
        # transfer pattern
        else:
            pseq.append( pattern[sindex] )
            qindex += 1
            sindex += 1

    # add opening / closing gaps
    opening_gaps = [pgap] * (qstart - 1)
    closing_gaps = [pgap] * (qlen - qend)
    pseq = opening_gaps + pseq + closing_gaps

    # check / return
    if len( pseq ) != len( query ):
        sys.exit( "query and aligned pattern not of the same length" )
    return pseq

def main():
    assert os.path.exists( sys.argv[0] ), \
        "if running as script arg1 is file 'aaseq1<\n>aaseq2<\n>pattern2'"
    print " ".join( align( *[k.strip() for k in open( sys.argv[1] ).readlines() ] ) )

if __name__ == "__main__":
    main()
