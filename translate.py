#!/usr/bin/python
from bio2 import Sequence, loadFasta
class Residue(object):
    def __init__(self):
        self.seq = []

    def __str__(self):
        return str(self.seq)

    def __eq__(self, x):
        return  self.seq == x.seq


standard = {'PHE': 'F', 'LEU': 'L', 'ILE':'I', 'MET':'M', 'MSE':'M', 'VAL':'V', 
            'SER':'S', 'PRO':'P','THR':'T','ALA':'A', 'TYR':'Y', 'HIS':'H', 
            'GLN':'Q', 'ASN':'N', 'LYS':'K', 'ASP':'D', 'GLU':'E', 'CYS':'C',
            'TRP':'W', 'ARG':'R', 'SER':'S','GLY':'G'}

inv =  {'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY', 'F': 'PHE',
        'I': 'ILE', 'H': 'HIS', 'K': 'LYS', 'M': 'MET', 'L': 'LEU', 'N': 'ASN',
        'Q': 'GLN', 'P': 'PRO', 'S': 'SER', 'R': 'ARG', 'T': 'THR', 'W': 'TRP',
        'V': 'VAL', 'Y': 'TYR'}



def getSEQRES(target):
    data = open(target, 'r')
    line  = data.readline()
    residues = []
    while line != "":
        if line[:6] == 'SEQRES':
            res =  Residue()
            length = int(line[13:17])
            processed = 0 
            line = line[19:]
            line = line.split()
            while True:
                for r in line:
                    res.seq.append(standard[r])
                    processed +=1
                if processed < length:
                    line = data.readline()[19:].split()
                else:
                    residues.append(res)
                    print res
                    break
        line = data.readline()
    return residues


def doesBContainA(a,b):
    for res in a:
        print hasResidue(res,b)


def hasResidue(target, residues):
    for res in residues:
        if target ==  res:
            return True
    return False

if __name__ == '__main__':
    target = "4HKD.pdb"
    residues = getSEQRES(target)
    fasta_seqs = loadFasta("pdbaanr")
    doesBContainA(residues,fasta_seqs)