# This will take the user input merged with deposited database genomes
# As a fasta file
# It will perform a multiple sequencfe alignment with ClustalW
# Then look at conserved regions
# From the user input, it will make a list of conserved seqences
# As long as they meet the minimum requirements
# Written by Andres W

# BEGIN CODE

from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

MERGED = "Ebola.fas"
ALIGNED = "Ebola.aln"
DND    = "Ebola.dnd"

INPUT     = ClustalwCommandline("clustalw", infile = MERGED)
INPUT()

ALIGNED = AlignIO.read(ALIGNED, "clustal")
# print ALIGNED

# !/usr/bin/python
 
# This program will scan through the ClustalW output
# It will tell you the length of the longest continuous sequence

 
zapato = open('Ebola.aln', 'r')
 
def uninterrupted(zapato):

	longest   = 0
	das_list  = []
	blank     = 0
	howmuch   = 0
	hone_it   = []
	max_me    = 0
	min_me    = 0
	range_it  = ()
	map_it    = []
	location  = 0

	for line in zapato:
		if line[0] == ' ':
			blank += 17
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
#		print map_it
#		print howmuch
#		print location
	return map_it



# The minimum length criteria can vary, for now, 17 nucleotide minimum shRNA
def place(map_it):

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
#		print WALDO

	return WALDO


# For ceviche to work the program has to assign a name to user input variable
# It can be called EBO for now
# It cannot have any aronym of nucleotide bases
# Just change the Tai in code below to name given to user sequence

def ceviche(zapato):

	zapato = open('Ebola.aln', 'r')

	picante = ''
	valid   = ['C','T','G','A','-']

	for line in zapato:
		if line[0] == 'B' and line[1] == 'u' and line[2] == 'n':
			for item in line:
				if item in valid:
					picante += str(item)
	return picante


def candidates():

	compare  = ceviche(zapato)
	contrast = place(uninterrupted(zapato))
	ready    = 0
	go       = 0
	CONSERVED = ''
	OUTPUT   = []

	for item in contrast:
		ready      = item[0]
		go         = item[1] + 1
		CONSERVED  = compare[ready:go]
#		print ready
#		print go
#		print CONSERVED
		OUTPUT.append(CONSERVED)
	return OUTPUT


#def off_to_blast():

#	try_me = candidates() 
#	count = 0
#	secuence = ''
#	last_list = []
#	for item in try_me:
#		print item
#		for value in item:
#			count += 1
#			secuence += item
#			if count >= 17:
#				last_list.append(secuence)	
#			   	count = 0
#				secuence = ''
					
#		print last_list		
		
	


#off_to_blast()


print candidates()

# DONT FORGET TO RUN THE PROGRAM!

raw_input('Press ENTER if you are done')
print "bye dre!"

