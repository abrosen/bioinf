from Bio.PDB import *
import matplotlib.pyplot as plt
import numpy as np

def getStruct(target):
    parser = PDBParser(PERMISSIVE=0)
    structure  = parser.get_structure('LLL', target)
    return structure
    
def getChain(structure, chain):
    model =  structure[0]
    return model[chain]

def getDists(chain):
    residues =  []
    for r in chain.get_residues():
        #dunno why water is a residue, unneeded here since I'm comparing CA atoms
        if r.get_resname() == 'HOH':
            continue
        residues.append(r)

    distancesMatrix = [] 
    for r1 in residues:
        distances = []
        for r2 in residues:
            distances.append(dist(r1,r2))
        distancesMatrix.append(distances)
    return distancesMatrix



def dist(r1, r2):
    """gets the distance between two residue objects
    defined by distance between 
    """
    #print r1, r2
    return r1['CA'] - r2['CA']


s = getStruct("4HKD.pdb")
c = getChain(s,'A')
d = getDists(c)
data =  np.array(d)


#for i in range(len(data)):
#    data[10][i] =5



# data.T for transpose
# pcolormesh faster than pcolor, and I'm demoing on a laptop
plt.pcolor(data,norm=None,cmap ="RdBu", edgecolors='k')   #(data[:,::-1])
plt.colorbar()



# display the graph
plt.show()  #savefig to save it instead
plt.close()