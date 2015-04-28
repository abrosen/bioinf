#emulate on Ruffus

import os

def ebola_conservative_regions (infile, outfile):
    current_region =0
    for line in open (infile):
        if line[0] == '>':
            current_region += 1
            current_file = open ("%d.%s" % (current_region, outfile), "w")
        current_file.write(line)

def runBlast (infile, outfile):
    os.system ("blastn -query %s -db /usr/local/blastdata/human_genomic -out %s.output -outfmt 0" % (infile,outfile))

#ebola_conservative_regions ('teflon_segs.txt', 'teflon_segment')

runBlast("Ebolatest.fas", "testfile")

#no_segments =11
#for segment in range(no_segments):
#    runBlast ('%d.teflon_segment' % segment, '%d.teflon_segment' % segment)

