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

                    #print ("test while T ", counts_hits, "  start: ", query_start, "  end: ", query_end)

                elif "Score" in hsp_line:
                    parse_flag = True
                    break

                elif "Lambda" in hsp_line:
                    parse_flag = True
                    break
                
            print ("test while T ", counts_hits, "  start: ", query_start, "  end: ", query_end)
            #print ("\n")

 
        """
        record the info extracted from each hsp entry to results
        
        """
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



"""
main program
"""
handle = open ("Ebola_blast_output.txt")
#handle = open ("blastp.txt")
test= scan_blastn(handle, identity_cutoff = 90)
fasta_output (test.counts_hits, 'Ebola_hits')
 
"""
check results:
test.counts_hits
test.results[0].query_sequence



"""
