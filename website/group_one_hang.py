# !/usr/bin/python

# Part 1 Step 1 : Hang Huynh
# pillomavirus_type_63_uid15486 This assume that the "virus genome files" is saved manually file by file
# grabbing files with common  ".fas" files in the same directories 
# input "common name" by user

import os, glob
from Bio import SeqIO
import fileinput

# Finding the common name of the virus genome in the file to be scan against
# the first letter should be capitalized (based on how i saved the file)
def fileDirectory(name):
    	file1 = '*.fna'
    	file2 = name    # this part can be change depends on how the data is save
    	file3 = (file2) + (file1)
    	file4 = '*' + (file3)                      # and if the user select "none" then it assume that the user 
    	filelist = glob.glob(file4)                # is having all the virus family genome in one file
    	with open ("virus.fas", "w") as outfile:   # or do not have the common name such as new virus type ect.
        	for f in filelist:
            		with open(f, "rU") as infile:
			 	seq_records = SeqIO.parse(infile, "fasta")
				SeqIO.write(seq_records, outfile, "fasta")
fileDirectory()
"""# This part create the file called "mergeall.fas" which merge all the file...
# with the common name "user input common name"

# This part merge the "Mergeall.fas" file with "userinput.fas" file to... 
# create the merging.fas file
def mergeDirectory():
	file5 = ["mergeall.fas","test.fas"]           # ".fas" is user input file
	with open ("virus.fas", "w") as outfile:    
		for f in file5: 
			with open(f, "rU") as infile:
				seq_records = SeqIO.parse(infile,"fasta")
				SeqIO.write(seq_records, outfile,"fasta")
mergeDirectory()
# the final file ouput called "virus.fas" as the input for the ClustalW alignment
"""

# Part 1 Step 2 : Andres Wong
# This will take the user input merged with deposited database genomes
# As a fasta file
# It will perform a multiple sequencfe alignment with ClustalW
# Then look at conserved regions
# From the user input, it will make a list of conserved seqences
# As long as they meet the minimum requirements

def alineacion():
	from Bio import SeqIO
        from Bio.Align.Applications import ClustalwCommandline
        from Bio import AlignIO
# This is the file from previous step, named Ebola.fas
        MERGED = "virus.fas"
        ALIGNED = "virus.aln"
        INPUT = ClustalwCommandline("clustalw", infile = MERGED)
        INPUT()

#ALIGNED = AlignIO.read(ALIGNED, "clustal")
# print ALIGNED
# End of Andrew R modified code
# This program will scan through the ClustalW output
# It will tell you the length of the longest continuous sequence
def uninterrupted():
        zapato = open('virus.aln', 'r')
        longest   = 0
        das_list  = []
        blank     = 0
        nada      = 0
        howmuch   = 0
        hone_it   = []
        max_me    = 0
        min_me    = 0
        range_it  = ()
        map_it    = []
        location  = 0
        for line in zapato:
                if line[0] == ' ':
                        blank += 36
                        for item in line:
                                howmuch += 1
                                if item == '*':
                                        longest += 1
                                        if longest >= 17:
                                                location = howmuch - blank
                                                hone_it.append(location)
                                                max_me = max(hone_it)
                                                min_me = min(hone_it)
                                if item != '*' and longest > 0:
                                        longest = 0
                                        range_it = (min_me,max_me)
                                        hone_it = []
                                if range_it not in map_it and range_it != (0, 0) and range_it != ():
                                        map_it.append(range_it)

        zapato.close()
#               print map_it
#               print howmuch
#               print location
        return map_it

# The minimum length criteria can vary, for now, 17 nucleotide minimum shRNA
def place():
        map_it = uninterrupted()
        primero   = 0
        segundo   = 0
        where_is  = ()
        WALDO     = []
        criteria  = 16
        for item in map_it:
                primero  = item[0] - criteria
                segundo  = item[1]
                where_is = (primero,segundo)
                WALDO.append(where_is)
#               print WALDO
        return WALDO

# For ceviche to work the program has to assign a name to user input variable
# It can be called EBO for now
# It cannot have any aronym of nucleotide bases
# Just change the Tai in code below to name given to user sequence

def ceviche():
        zapato = open('virus.aln', 'r')
        picante = ''
        valid   = ['C','T','G','A','-']
        for line in zapato:
                if line[0] == 'g' and line[1] == 'i' and line[2] == '|' and line[3] == '1':
                        for item in line:
                                if item in valid:
                                        picante += str(item)
        zapato.close()
        return picante

def candidates():
        picante  = ceviche()
        WALDO = place()
        ready    = 0
        go       = 0
        CONSERVED = ''
        OUTPUT   = []
        for item in WALDO:
                ready      = item[0]
                go         = item[1] + 1
                CONSERVED  = picante[ready:go]
                OUTPUT.append(CONSERVED)
        return OUTPUT

# Below is in FASTA format
# Gives te sequences to BLAST against human genome

def split_em():
        OUTPUT = candidates()
        TERMINATOR = []
        original = ''
        for item in OUTPUT:
	# print 'This is the original long sequence %s and is %s long' % (item,len(item))
                while len(item) > 16:
                        TERMINATOR.append(item)
                        original = item
                        item = item[1:]
# Code below compensates for item being appended earlier
# This way, the length is one more than the original item
                        while len(original) > 17:
                                original = original[0:-1]
#                               print original
                                TERMINATOR.append(original)
        return TERMINATOR

def to_blast():
        alineacion()
        TERMINATOR = split_em()
        blast_me = open('andres.fas','a')
        fasta_id = 0
        for item in TERMINATOR:
                fasta_id += 1
                blast_me.write('>pre_candidate_%s\n' % (fasta_id))
                blast_me.write('%s\n' % (item))
        blast_me.close()

# Read the file andres.fas for the next step for BLAST
# User interface below
to_blast()
"""
#IIIIIIII
# parser of alignment from blastn output
# based on the currently deprecated biopython/NCBIStandalone.py
# Biopython's old parser of blastall output
# Small parts of the original codes are modified to fit our needs


   # blastn_scanner: input a blastn_output file, looking for info we want to parse.



class blastn_scanner (object):
    def __init__ (me,
                  version,                  #what version of blastn used?
                  database_name,             #what db used for blastn?
                  database_sequences,       #number of sequences in the target db
                  database_letters,        #total size (in basepairs) of the db
                  counts_hits,              #number of targets hit
                  results):                 #holder for all the info we want parsed for next step
                                            #results contains info produced from class blastn_result
        me.version = version
        me.database_name = database_name
        me.database_sequences = database_sequences
        me.database_letters = database_letters
        me.counts_hits = counts_hits
        me.results = results

	# blastn_results: holder of all info from HSP

class blastn_results (object):
    def __init__ (me,
                  identities,                #% identities; 100% for siRNA
                  strand,                    #Plus or Minus
                  query_sequence,            #sequences matched to target, saved for next step
                  query_start,
                  query_end):                 #coordinates of the query HSP - not much meaningful for us
        me.identities = identities
        me.strand = strand
        me.query_sequence = query_sequence
        me.query_start = query_start
        me.query_end = query_end
def scan_blastn (blastn_output_file, identity_cutoff = 90):
    	line = blastn_output_file.readline()
    	if not line.startswith("BLASTN"):
        	raise TypeError ('Not a BLASTN output file')
    	version = line.rstrip()                 #store "BLATN 2.2.30+"
    	while True:
        	line = blastn_output_file.readline()
        	if line.startswith("Database:"):
            	database_name = line.rstrip().split()[1] #store the name of HG
            	break
    	while True:
        	line = blastn_output_file.readline()
        	if "sequences" in line or "total letters" in line:
            	words = line.rstrip().split()
            	database_sequences = int (words[0].replace (",","") )    #how many sequences in db?
            	database_letters = int (words[2].replace (",","") )       #how big is the db?
            	break
    	while True:                                             #locate the benchmark of hits or no hits
        	line = blastn_output_file.readline()
        	if "No hits found" in line:
            		break                           #exit()? to terminate
        	elif line.startswith ("Sequences producing significant alignments:"):
            		break
        	elif not line:
            		raise TypeError ("No description section in blastn output")
    		# Move to each entry of HSP
    
    	results = []
    	counts_hits =0
    	parse_flag = False
    	while True:
        	hsp_line = blastn_output_file.readline()
        	if not hsp_line:
            		break
        if "Identities" in hsp_line:                        #extract %identities
            	parse_flag = False
            	counts_hits += 1
            	temp = hsp_line.rstrip().split()[3].replace("(","")
            	identities = int (temp.replace("%),", ""))
        	if "Strand" in hsp_line:                            #extract strand polarity
            temp = hsp_line.rstrip().replace("/", " ")
            strand = temp.split()[-1]
        if hsp_line.startswith ("Query"):                   #reach Query line
            temp = hsp_line.rstrip().split()
            query_start = int(temp[1])
            query_sequence = temp[2]
            query_end = int(temp[3])
            while True:
                hsp_line = blastn_output_file.readline()
                if not hsp_line:
                    break
                if hsp_line.startswith("Query"):
                    temp = hsp_line.rstrip().split()
                    start = int (temp[1])
                    #if abs (start - query_end) !=1:
                    #raise TypeError ("something wrong with the Query sequence positions")
                    if start < query_start:
                        query_start = start
                        raise TypeError ("start: minus Query or something wrong")
                    query_sequence = query_sequence + temp[2]

                    end = int (temp[3])
                    if end > query_end:
                        query_end = end
                    else:
                        raise TypeError ("end: minus Query or something wrong")
                    #print ("test while T ", counts_hits, "  start: ", query_start, "  end: ", query_end)
                elif "Score" in hsp_line:
                    parse_flag = True
                    break
                elif "Lambda" in hsp_line:
                    parse_flag = True
                    break
            print ("test while T ", counts_hits, "  start: ", query_start, "  end: ", query_end)
            #print ("\n")

       # record the info extracted from each hsp entry to results

        if parse_flag:
            result = blastn_results (identities, strand, query_sequence, query_start, query_end)
            results.append (result)
    print ("end of program")
    return blastn_scanner (version, database_name, database_sequences, database_letters, counts_hits, results)

def fasta_output (records, fasta_output_file):
    target = open (fasta_output_file, 'w')
    for i in range (records):
        target.write (">hits%d \n" % i)
        target.write (test.results[i].query_sequence.replace("-",""))
        target.write ("\n")

handle = open ("Ebola_blast_output.txt")
#handle = open ("blastp.txt")
test= scan_blastn(handle, identity_cutoff = 90)
fasta_output (test.counts_hits, 'Ebola_hits')

# IIIIIIII
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
    	os.system ("blastn -query %s -db NT_033778.fna -out %s.output -outfmt 0" % (infile,outfile))
ebola_conservative_regions ('teflon_segs.txt', 'teflon_segment')
no_segments =11
for segment in range(no_segments):
    	runBlast ('%d.teflon_segment' % segment, '%d.teflon_segment' % segment)



# IIIIIIII
# This program was written to handle step 5
# It searches the ebola genome input from user
# Looking for the hits from the BLAST search
# And eliminates the human hits from the user viral genome
# So this kind of precedes my original stepG
# The output is an ebola genome without regions of human alignments
# In FASTA format
# Written by Andres W

def separalo():
        pre_file = open('Ebola_Human_Hits.fas','r')
        DNA = list('GACT')
        segment_list = []
        temp_seq = ''
        temp_seq2 = ''
        for line in pre_file:
                if line[0] != ">":
                        temp_seq = temp_seq + str(line)
                        for item in temp_seq:
                                if item in DNA:
                                        temp_seq2 = temp_seq2 + str(item)
                        segment_list.append(temp_seq2)
                        temp_seq = ''
                        temp_seq2 = ''
        pre_file.close()
        return segment_list

def virus():
        pre_file = open('Ebola_genome.fas','r')
        DNA = list('GACT')
        pre_viral_genome = ''
        viral_genome = ''
        for line in pre_file:
                if line[0] != ">":
                        pre_viral_genome += str(line)
                        for item in pre_viral_genome:
                                if item in DNA:
                                        viral_genome += str(item)
        pre_file.close()
        return viral_genome

def comparalo_via_clustal():
        from Bio import SeqIO
        from Bio.Align.Applications import ClustalwCommandline
        from Bio import AlignIO
        MERGED = "Ebola_virus_plus_hits.fas"
        ALIGNED = "Ebola_virus_plus_hits.aln"
        INPUT     = ClustalwCommandline("clustalw", infile = MERGED)
        INPUT()
        ALIGNED = AlignIO.read(ALIGNED, "clustal")

def comparalo_via_re():
        import re
        segment_list = separalo()
        viral_genome = virus()
        edited_genome = ''
        removed = ''
        for segment in segment_list:
                if segment in viral_genome:
                        edited_genome = re.sub(segment, '\n', viral_genome)
                        viral_genome = edited_genome
        return edited_genome

def double_check():
        edited_genome = comparalo_via_re()
        segment_list = separalo()
        count = 0
        for segment in segment_list:
                if segment not in edited_genome:
                        count += 1
        print count

def to_group2():
        edited_genome = comparalo_via_re()
        group1_output = open('dre.fas','a')
        group1_output.write('>Ebola_genome_minus_human_blastn_hits\n')
        cuenta = 0
        for base in edited_genome:
                cuenta += 1
                group1_output.write('%s' % (base))
                if cuenta == 70:
                        cuenta = 0
                        group1_output.write('\n')
        group1_output.close()
to_group2()
"""
