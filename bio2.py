class Sequence(object):
    def __init__(self, comment,contents):
        self.comment = comment
        self.data = contents


def numUniqueChars(s):
    return len(set(s))



def numUniqueWords(sentence):
    sentence = sentence.split()
    return len(set(sentence))


def numUniqueWordsInSentences(sentences):
    return numUniqueWords(" ".join(sentences))


def loadFasta(filename):
    fasta =  open(filename, 'r')
    

def main():
    print numUniqueChars("hello")
    print numUniqueWords("Hello to all my friends out out there.")
    print numUniqueWordsInSentences(["hello","from","mars"])
    loadFasta("pdbaanr")

if __name__ == "__main__":
    main()
