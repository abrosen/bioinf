# !/usr/bin/python

# Finalized group 1 code
# COMPILED BY HANG and ANDRES
#
#
# STEP 1: YUN TAO
# This step is independent from STEP 1 and STEP 2
# Blasts the ebola genome against the human genome
# Return hits for STEP 4


import os

def ebola_conservative_regions (infile, outfile):
    current_region =0
    for line in open (infile):
        if line[0] == '>':
            current_region += 1
            current_file = open ("%d.%s" % (current_region, outfile), "w")
        current_file.write(line)

def runBlast (infile, outfile):
    os.system ("blastn -word_size 11 -query %s -db blastdata/human_genomic -out %s.output -outfmt 0" % (infile,outfile))

#ebola_conservative_regions ('teflon_segs.txt', 'teflon_segment')

#no_segments =11
#for segment in range(no_segments):
runBlast ('username.fas', 'ebola')



# parser of alignment from blastn output
# based on the currently deprecated biopython/NCBIStandalone.py
# Biopython's old parser of blastall output
# Small parts of the original codes are modified to fit our needs


"""
    blastn_scanner: input a blastn_output file, looking for info we want to parse.

"""


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
"""
    blastn_results: holder of all info from HSP
    
"""

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

    """
    Move to each entry of HSP
    """

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


                elif "Score" in hsp_line:
                    parse_flag = True
                    break

                elif "Lambda" in hsp_line:
                    parse_flag = True
                    break



        """
        record the info extracted from each hsp entry to results
        """
        if parse_flag:
            result = blastn_results (identities, strand, query_sequence, query_start, query_end)
            results.append (result)


    return blastn_scanner (version, database_name, database_sequences, database_letters, counts_hits, results)

def fasta_output (records, fasta_output_file):
    target = open (fasta_output_file, 'w')
    for i in range (records):
        target.write (">hits%d \n" % i)
        target.write (test.results[i].query_sequence.replace("-",""))
        target.write ("\n")



"""
main program
"""
handle = open ("ebola.output")
#handle = open ("blastp.txt")
test= scan_blastn(handle, identity_cutoff = 90)
fasta_output (test.counts_hits, 'ebola_blastn.fas')

"""
check results:
test.counts_hits
test.results[0].query_sequence



"""

# Part 1 Step 1 : Hang Huynh
# pillomavirus_type_63_uid15486 This assume that the "virus genome files" is saved manually file by file
# grabbing files with common  ".fas" files in the same directories 
# input "common name" by user


def blastin():
	
	pre_file = open('ebola_blastn.fas','r')	
	
	DNA = list('GACT')
        blastned = []
        temp_seq = ''
	temp_seq2 = ''

        for line in pre_file:
                if line[0] != ">":
			temp_seq = temp_seq + str(line)
			for item in temp_seq:
				if item in DNA: 
					temp_seq2 = temp_seq2 + str(item)
			blastned.append(temp_seq2)
	      	 	temp_seq = ''
			temp_seq2 = ''

	pre_file.close()

        return blastned

def editin():

        pre_file = open('username.fas','r')

        DNA = list('GACT')
        pre_viral_genome = ''
        viral_genome = ''

        for line in pre_file:
                if line[0] != ">":
                        pre_viral_genome += str(line)
                        for item in pre_viral_genome:
                                if item in DNA:
                                        viral_genome += str(item)
			pre_viral_genome = ''

        pre_file.close()

        return viral_genome

def reducir():


        segment_list = blastin()
        viral_genome = editin()

        genoma_editado = ''

        for segment in segment_list:
                if segment in viral_genome:
                        genoma_editado = re.sub(segment, '\n', viral_genome)
                        viral_genome = genoma_editado


        return viral_genome


def fasta():

	modified = reducir()		
	merge = open('ebola_edited.fas','w')

	merge.write('>user_input\n')	
	for line in modified:
		merge.write(line)
	
	merge.close()

# Finding the common name of the virus genome in the file to be scan against
# the first letter should be capitalized (based on how i saved the file)
def fileDirectory():
        file1 = '*.fna'
#        file2 = raw_input('input name of virus in lowercase:')    #ANDREW LOOK HERE PLEASE 
	file2 = 'ebola'
        file3 = (file2) + (file1)
        file4 = '*' + (file3)                      # and if the user select "none" then it assume that the user 
        filelist = glob.glob(file4)                # is having all the virus family genome in one file
        with open ("mergeall.fas", "w") as outfile:   # or do not have the common name such as new virus type ect.
                for f in filelist:
                        with open(f, "rU") as infile:
                                seq_records = SeqIO.parse(infile, "fasta")
                                SeqIO.write(seq_records, outfile, "fasta")
		

# This part create the file called "mergeall.fas" which merge all the file...
# with the common name "user input common name"

# This part merge the "Mergeall.fas" file with "userinput.fas" file to... 
# create the merging.fas file
def mergeDirectory():
        file5 = ["ebola_edited.fas","mergeall.fas"]           # ".fas" is user input file
        with open ("virus.fas", "w") as outfile:    
                for f in file5: 
                        with open(f, "rU") as infile:
                                seq_records = SeqIO.parse(infile,"fasta")
                                SeqIO.write(seq_records, outfile,"fasta")


# the final file ouput called "virus.fas" as the input for the ClustalW alignment
# Part 1 Step 2 : Andres Wong
# This will take the user input merged with deposited database genomes
# As a fasta file
# It will perform a multiple sequencfe alignment with ClustalW
# Then look at conserved regions
# It will scan a window moving over 1 base every time 
# And provide conserved regions accounting for half coverage
# Or accounting for regions that have over n% identity
# As a FASTA format
# Then will be fed by yun's step
# And sequences will be removed as needed
# Then the last file will be sent to Group 2
# Written by Andres W
def reformat():
        shoe = open('virus.fas','r')
        zapato = open('ebola.fas', 'w')
        header = 0
        for line in shoe:
	        if line[0] == '>':
        	        header += 1
                        zapato.write('>virus_%s\n' % (header))
                else:
                        zapato.write(line)
        zapato.close()

def alineacion():
	MERGED = "ebola.fas"
	ALIGNED = "ebola.aln"
	INPUT = ClustalwCommandline("clustalw", infile = MERGED)
	INPUT()

def asterix():
	zapato = open('ebola.aln', 'r')
	blank     = 1
	howmuch   = 0
	WHERE_IS  = []
	location  = 0
	for line in zapato:
		if line[0] == ' ' and line[1] == ' ' and line[2] == ' ':
			blank += 17
			for item in line:
				howmuch += 1
				if item == '*':
					location = howmuch - blank
					WHERE_IS.append(location)
	zapato.close()
	return WHERE_IS

def stringed_genome():
	zapato = open('ebola.aln', 'r')
	WALDO = ''
	valid   = ['C','T','G','A','-']
	for line in zapato:
		if 'virus_1' in line:
			for item in line:
				if item in valid:
					WALDO += str(item)
	zapato.close()
	return WALDO 

def slip_n_slide():
	WHERE_IS = asterix()	
	lower = WHERE_IS
	WALDO    = stringed_genome()
	dresifus = ''
	conde  = 0.0 
	cuenta = 0.0
	GANAMOS = []
	campeones  = []
	grones = [] 
	lima = 0.0
	sunny   = 1 
	window  = 0 
	alianza = 100.0 #$$$
	hold_em = alianza #$$$
	min_length = 19.0 #$$$
	while hold_em >= min_length: 
                for nucleotide in WALDO:
			cuenta += 1
                        dresifus += nucleotide
                        if nucleotide == '-' and cuenta <= alianza:
                                dresifus = dresifus[:-1]
                                recuenta = cuenta - 1.0
                                if len(dresifus) >= min_length:
                                        campeones.append(dresifus)
                                        for star in lower:
						asterisco = star + conde
                                                if asterisco >= sunny and asterisco < (window + cuenta):
                                                        lima += 1
                                                	abajo = lower.index(star)
                                                if asterisco >= (window + cuenta):
                                                        break
                                        homology = lima/recuenta
                                        grones.append(homology)
                                        lower = lower[abajo:]
                                sunny += cuenta
                                window += cuenta
                                lima = 0
                                dresifus = ''
                                cuenta = 0
			elif cuenta == alianza:
				campeones.append(dresifus)
				for star in lower:
					asterisco = star + conde
					if asterisco >= sunny and asterisco <= (window + alianza):
						lima += 1
						abajo = lower.index(star)
					if asterisco > (window + alianza):
						break
				homology = lima/alianza
				grones.append(homology)
				lower = lower[abajo:]
				sunny += alianza
				window += alianza
				lima = 0
				dresifus = ''
				cuenta = 0
		lower = WHERE_IS
		cuenta  = 0.0
		sunny   = 1.0
		window  = 0.0
		dresifus = ''
		WALDO = WALDO[1:]
		hold_em = len(WALDO)
		conde += 1.0
		if hold_em < min_length:  
			criteria = numpy.median(grones)			
			win = -1			
			for candidate in grones:
				win += 1
				if candidate >= .6: #$$$
					i_did_it = campeones[win]
					GANAMOS.append(i_did_it)
	return GANAMOS 
def find_unix():
	GANAMOS = slip_n_slide()
	longway = []
	for item in GANAMOS:
		if item not in longway:
			longway.append(item)
	GANAMOS = longway
	return GANAMOS
def have_a_blast():
	reformat()
	alineacion()
	GANAMOS = find_unix()
	INCA = open('andres.fas','w')
	fasta_id = -1 	
	for GOL in GANAMOS:
		fasta_id += 1
		INCA.write('>weapon_%s\n' % (fasta_id))
		INCA.write('%s\n' % (GOL))	
	INCA.close()
def group_one():
	fasta()
	fileDirectory()
	mergeDirectory()
	have_a_blast()
# DONT FORGET TO RUN THE PROGRAM!
import re
import os, glob
from Bio import SeqIO
import fileinput
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import numpy
group_one()	
raw_input('Your file has been saved as andres.fas\nPress ENTER if you are done')
# END GROUP 1 CODE




