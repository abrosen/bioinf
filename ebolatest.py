from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio import Phylo

TARGET  = "Ebola.fas"
DEST  = "Ebola.aln"
DND =  "Ebola.dnd"

"""
command  = ClustalwCommandline("clustalw", infile=TARGET)
command()
print(command)
"""



align = AlignIO.read(DEST, "clustal")
print align

tree = Phylo.read(DND, "newick")
Phylo.draw_ascii(tree)



"""
records = []
for record in SeqIO.parse(TARGET,'fasta'):
    records.append(record)

for record in records:
    print record.id, record.seq[:10]
"""


