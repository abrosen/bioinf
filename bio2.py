import re
class Sequence(object):
    def __init__(self, comment,contents):
        self.header = comment
        self.data = contents
    def __str__(self):
        return self.header + '\n' + self.data

def numUniqueChars(s):
    return len(set(s))



def numUniqueWords(sentence):
    sentence = sentence.split()
    return len(set(sentence))


def numUniqueWordsInSentences(sentences):
    return numUniqueWords(" ".join(sentences))


def loadFasta(filename):
    fasta =  open(filename, 'r')
    line =  fasta.readline()
    sequences = []
    header  = ""
    data = []
    i =  0
    while line != "":
        line=line.strip()
        if line == "" and header != "":
            data = "".join(data)
            sequences.append(Sequence(header,data))
            header = ""
            data = []
        elif line == "":
            pass
        elif line[0] == ">":
            header = line
        elif line != '':
            data.append(line)
        line = fasta.readline()
        
    print "Done.  Loaded", len(sequences), "sequences."
    return sequences

def findHISTags(sequences):
    tags = []
    pattern = '^[^H]{2,3}H{5,8}'
    p = re.compile(pattern)
    for seq in sequences:
        m = p.match(seq.data)
        if m is not None:
            tag =  m.group()
            if True or tag not in tags:
                tags.append(tag)
    return tags

def main():
    print numUniqueChars("hello")
    print numUniqueWords("Hello to all my friends out out there.")
    print numUniqueWordsInSentences(["hello","from","mars"])
    s= loadFasta("pdbaanr")
    tags = findHISTags(s)
    print len(tags)

if __name__ == "__main__":
    main()
