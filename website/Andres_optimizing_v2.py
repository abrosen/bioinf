
# !/usr/bin/python

	

# This will take the user input merged with deposited database genomes
# As a fasta file
# It will perform a multiple sequencfe alignment with ClustalW
# Then look at conserved regions
# It will scan a window moving over 1 base every time 
# And provide conserved regions accounting for half coverage
# As a FASTA format to blastn on third step
# Written by Andres W

# BEGIN CODE

def reformat(filename):

	shoe = open(filename,'r')
	zapato = open('viral.fas', 'a')
	header = 0

	for line in shoe:
		if line[0] == '>':
			header += 1 
			zapato.write('>virus_%s\n' % (header))	
		else:
			zapato.write(line)
	
	zapato.close()

# Part below modified from Andrew R code

def alineacion():

	from Bio import SeqIO
	from Bio.Align.Applications import ClustalwCommandline
	from Bio import AlignIO

# This is the file from previous step, named Ebola.fas
	MERGED = "viral.fas"
	ALIGNED = "viral.aln"

	INPUT = ClustalwCommandline("clustalw", infile = MERGED)
	INPUT()

#	ALIGNED = AlignIO.read(ALIGNED, "clustal")
#	print ALIGNED
# End of Andrew R modified code
# This program will scan through the ClustalW output
# It will tell you the length of the longest continuous sequence
			  

def asterix():
	
	zapato = open('viral.aln', 'r')

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
#		print map_it
#		print howmuch
#		print location
	return WHERE_IS



# The search window will be 50 bases, for now

	
# It cannot have any aronym of nucleotide bases
# Just change the Tai in code below to name given to user sequence

def stringed_genome():

	zapato = open('viral.aln', 'r')

	WALDO = ''
	valid   = ['C','T','G','A','-']

	for line in zapato:
		if 'virus_1' in line:
			for item in line:
				if item in valid:
					WALDO += str(item)
#	print len(WALDO)

	zapato.close()
# WALDO is just the user viral genome converted from FASTA to string
	return WALDO 


def slip_n_slide():

	import numpy

	WHERE_IS = asterix()	
	WALDO    = stringed_genome()

	dresifus = ''
	cuenta   = 0
	GANAMOS = []
	campeones  = []
	grones = [] 

# Only need the three below if want to get average of homologies
#	todos = 0
#	promedio = 0
#	vamos = 0

	lima = 0
	
	sunny   = 1 

# $$$
	alianza = 50.0

# $$$
	window  = 0 
# $$$
	hold_em = alianza 
# $$$
	min_length = 18

		 
	dash = 0

# Need to figure out if complete coverege of last segment smaller than 50 is needed
# If so, it is easy coding
# Just remove commented out line below
	
	while hold_em >= alianza: 
		for nucleotide in WALDO:
			cuenta += 1
			dresifus += nucleotide 
			if nucleotide == '-' and cuenta <= alianza :
				dash += 1
				dresifus = dresifus[:-1]
				recuenta = cuenta - 1.0
				if len(dresifus) >= min_length:
					campeones.append(dresifus)
					for asterisco in WHERE_IS:
						if asterisco >= sunny and asterisco < (window + cuenta):
							lima += 1
					homology = lima/recuenta
					grones.append(homology)
				sunny += cuenta
				window += cuenta 
				lima = 0
				dresifus = ''
				cuenta = 0
				
			elif cuenta == alianza:
				campeones.append(dresifus)	
				for asterisco in WHERE_IS:
					if asterisco >= sunny and asterisco <= (window + alianza):
						lima += 1
				homology = lima/alianza
				grones.append(homology)
				sunny += alianza
				window += alianza
				lima = 0
				dresifus = ''
				cuenta = 0
					
		cuenta   = 0
		sunny   = 1 
		window  = 0 
		dresifus = ''
		WALDO = WALDO[1:]
		hold_em = len(WALDO)
	
		if hold_em < alianza:  
# criteria is the fifty percent that is most conserved 
# That means, the regions with homologies over the median
# But I will also include the median
			
			print grones	
			print len(grones)
			print len(stringed_genome())

			criteria = numpy.median(grones)			

			print criteria
			print dash
			
			win = -1			
#			over = 0
#			under = 0 
#			median = 0

			for candidate in grones:
				win += 1
				if candidate >= criteria:
# STAY AWAY FROM INDEX FUNCTION DUE TO IDNETICAL VALUES
					i_did_it = campeones[win]
					GANAMOS.append(i_did_it)
#				if item > criteria:
#					over += 1
#				if item < criteria:
#					under += 1
#				if item == criteria:
#					median += 1
#			median_test = over + median + under
# if criteria is have the the max value
#			half_max = max(grones)/2
# if criteria is the mean value
#			promedio = todos/vamos
			
#	print WHERE_IS
#	print campeones			
#	print cuenta
#	print grones 
#	print WALDO
#	print hold_em
# the below lengths for WHERE_IS, grones, and median_test should always match
#	print len(WHERE_IS)
#	print len(grones)
#	print median_test
# vamos and len(average) should always match
#	print vamos
#	print len(average)
# goles should be the total length of the sequence minus 1 minus the size of window
#	print goles	
#	print promedio
#       print half_max
#	print criteria
#	print win
#	print i_did_it
#	print slip_n_slide
	return GANAMOS 

# Below will return in FASTA format
# Gives the conserved sequences to BLAST against human genome

def have_a_blast():

	#reformat(filename)
	#alineacion()

	GANAMOS = slip_n_slide()
	print len(GANAMOS)
	INCA = open("andres.fas",'w')

	fasta_id = -1 	

	for GOL in GANAMOS:
		fasta_id += 1
		INCA.write('>pre_candidate_%s\n' % (fasta_id))
		INCA.write('%s\n' % (GOL))	
		
	INCA.close()
	
# Read the file andres.fas for the next step for BLAST

# DONT FORGET TO RUN THE PROGRAM!

# Final program to run is to_blast(), which is commented out below

if __name__ == '__main__':
	have_a_blast()

