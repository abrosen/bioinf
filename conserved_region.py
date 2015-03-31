# !/usr/bin/python

# This program will scan through the ClustalW output
# It will tell you the length of the longest continuous sequence
# Written by Andres W

import sys

target = sys.argv[1]
print target

das_boot = open(target, 'r')
# das_boot = open('Clustal', 'r')
# Above commented out to test program on smalled sequence

count    = 0
longest  = 0 
das_list = []

for line in das_boot:
	if line[0] == ' ':
		for item in line:
			if item == '*':
				count += 1
				longest += 1
			elif item != '*' and longest > 0:
				das_list.append(longest)
				longest = 0
				
print 'The sequences have %s conserved bases' % (count)
print 'Below is a list of lengths for each continous conserved sequence:'
print das_list

das_max = max(das_list)
print 'The longest continous conserved sequence is %s bp long' % (das_max)
raw_input('Press ENTER if you are done')
