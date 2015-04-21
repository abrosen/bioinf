#emulate on Ruffus
 
import os
 
def shRNA_check (self,sequence ,output):
    current_region =0
    self.sequence= sequence
    for list in open (self.sequence):
        if list[0] is not '':
            current_region += 1
            current_file = open ("%d.%s" % (current_region,output), "w")
        current_file.write(line)
 
def Blast_check (sequence, output):
    os.system ("blastn -query %s -db NT_033778.fna -out %s.output -outfmt 0" % (sequence,output))

print shRNA_check ('seq', 'sequence_segment')

for segment in range(len('sequence_segment')):
    Blast_check ('sequence_segment', 'sequence_seg_checked')
