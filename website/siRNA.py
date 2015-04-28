'''
Team Memebrs: Zhibo Wang
              Li Xu
              Chao Zhao
              Yuchen Zhang
'''
import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC
import subprocess
import random


class siRNA():
    # This funciotn is used to extract the sequence from Group 1's result
    def read_output(self, file):
        seqs = []
        temp = []
        lines = [line.strip() for line in open(file) if not line.startswith('>')]
        for x in range(len(lines)):
            for y in range(len(lines[x]) - 21):
                str = lines[x][y:y + 21]
                seqs.append(str)

        # if add length 17 18 19 20?

        # check any duplicate sequence
        for i in seqs:
            if i not in temp:
                temp.append(i)
        return temp

    def score_sys(self, seq_list, utr):
        dict_score = {}
        keys_to_remove = []
        self.seq_list = seq_list
        for x in self.seq_list:
            score = 0
            if x[0:2] == 'aa':
                score = score + 15
            elif x[1] == 'a':
                score = score + 8
            else:
                score = score

            if GC(x) >= 35 and GC(x) <= 45:
                score = score + 5
            elif GC(x) > 45 and GC(x) <= 55:
                score = score + 7
            elif GC(x) > 55 and GC(x) <= 65:
                score = score + 10
            else:
                score = score
            if 'aaaa' and 'cccc' and 'gggg' and 'tttt' not in x:
                score = score
            else:
                score = score - 1000
            if x not in utr:
                score = score
            else:
                score = score - 1000

            dict_score[x] = score
            #print dict_score

        keys_to_remove = [key for key, value in dict_score.iteritems() if value <= 0]

        for key in keys_to_remove:
            del dict_score[key]

        my_dict = sorted(dict_score, key=dict_score.get, reverse=True)[:10]
        '''
        for x in my_dict:
            print dict_score[x]
        '''
        print my_dict
        return my_dict

    # b)loop (TCAAGAG) for psiRNA vectors
    def design_shrna(self, seqs):
        self.seqs = seqs
        print self.seqs
        shrnas = []
        for seq_temp in self.seqs:
            seq =Seq(seq_temp)
            complement_seq = seq.complement()
            reversed_complement_seq = complement_seq[::-1]
            shrna = seq + 'tcaagag' + reversed_complement_seq
            shrna = shrna.transcribe()
            shrna = shrna.upper()
            shrnas.append(shrna)
        return shrnas

    # call RNAfold
    def check_2nd_structure(self, shrna_list):
        save_path = '/Users/zhibo/PycharmProjects/siRNA/'
        complete_name = os.path.join(save_path, 'inputfile.txt')
        with open(complete_name, 'w') as inputfile:
            i = 0
            for x in shrna_list:
                inputfile.write('>seq' + str(i) + '\n' + str(x) + '\n')
                i += 1
        subprocess.call('RNAfold < inputfile.txt > output.txt', shell=True)
        '''
        for x in range(len(shrna_list)):
            subprocess.call('open seq' + str(x) + '_ss.ps', shell=True)
        '''
    def determine_2nd_structure(self):
        with open ('output.txt', 'r') as f:
            name_list = []
            seq_list = []
            final_keys = []
            temp_keys = []
            temp_value = []
            final_seqs = []
            for line in f:
                if line.startswith('>'):
                    name_list.append(line[1:].replace('\n', ''))
                if line.startswith('(') or line.startswith('.'):
                    seq_list.append(line.split(' ')[0])

            seq_dict = dict(zip(name_list,seq_list))
            for key, value in seq_dict.items():
                if len(seq_dict[key].split('(.'))==2:
                    final_keys.append(key)
            with open ('inputfile.txt','r') as f:
                for line in f:
                    if line.startswith('>'):
                        temp_keys.append(line[1:].replace('\n',''))
                    else:
                        temp_value.append(line.replace('\n',''))
            seq_dict = dict(zip(temp_keys,temp_value))
            for x in final_keys:
                final_seqs.append(seq_dict[x])

            return final_seqs

    #--------Yuchen------------------------------------
    def count(self, a):
        a = a.lower()
        global seq
        seq = {'a': 0.0, 't': 0.0, 'c': 0.0, 'g': 0.0}
        seqp = dict()
        global count1
        count1 = 0
        # Count the number of each base in the sequence
        for letters in a:
            seq[letters] = seq[letters] + 1
            count1 = count1 + 1
            # Get the percentage of each base in the sequence
        for keys in seq:
            seqp[keys] = float(float(seq[keys]) / count1)
            # print seqp
        return seqp

    def matrix(self):
    # The original matrix. Mutation rate is 10^-4 per base per generation. GC has a twice rate than normal.
        ma = [
            [1 - 0.0004, 0.0001, 0.0001, 0.0002],
            [0.0001, 1 - 0.0004, 0.0002, 0.0001],
            [0.0001, 0.0002, 1 - 0.0004, 0.0001],
            [0.0002, 0.0001, 0.0001, 1 - 0.0004],
        ]
        return ma

    def randma(self, ma):
    # make a random matrix, which has a randomly 1-4 times change from the original one.
        ranma = [[0 for x in range(4)] for x in range(4)]
        for m in range(0, 4):
            for n in range(0, 4):
                if m == n:
                    ranma[m][n] = 1
                else:
                    ranma[m][n] = ma[m][n] * random.uniform(0.99, 1.01)
                    # To generate numbers for [1][1],[2][2],[3][3],[4][4]
        for i in range(0, 4):
            for j in range(0, 4):
                if i != j:
                    ranma[i][i] = ranma[i][i] - ranma[j][i]
        return ranma

    def muti(self, count, ma):
    # matrix mutiply
        global new_seq
        global judge
        judge = True
        new_seq = dict()
        new_seqp = {'a': 0.0, 't': 0.0, 'c': 0.0, 'g': 0.0}
        new_seqp['a'] = count['a'] * ma[0][0] + count['t'] * ma[0][1] + count['c'] * ma[0][2] + count['g'] * ma[0][3]
        new_seqp['t'] = count['a'] * ma[1][0] + count['t'] * ma[1][1] + count['c'] * ma[1][2] + count['g'] * ma[1][3]
        new_seqp['c'] = count['a'] * ma[2][0] + count['t'] * ma[2][1] + count['c'] * ma[2][2] + count['g'] * ma[2][3]
        new_seqp['g'] = count['a'] * ma[3][0] + count['t'] * ma[3][1] + count['c'] * ma[3][2] + count['g'] * ma[3][3]

        for keys in new_seqp.keys():
            new_seq[keys] = new_seqp[keys] * count1

            if seq[keys] - new_seq[keys] > 1 or seq[keys] - new_seq[keys] < -1:
                judge = False

        #print new_seq
        return new_seqp

    def mutation(self, list):
        dict = {}
        for x in list:
            i = 1
            count = self.count(x)
            ranma = self.matrix()
            new_count = self.muti(count, ranma)
            while judge == True:
                #print 'The sequence at ' +str(i+1)+"'s generation is:"
                ranma = self.randma(ranma)
                new_count = self.muti(new_count, ranma)
                i = i + 1
            dict[x] = i
        return dict

    def my_print(self, dict):
        for key, value in dict.items():
            print 'The sequence' + ' ' + key + ' ' + 'will change at' + ' ' + str(value) + '\'s generation!'

    #-------------------------------------------------------
    '''
    def ebola_conservative_regions (infile, outfile):
        current_region =0
        for line in open (infile):
            if line[0] == '>':
                current_region += 1
                current_file = open ("%d.%s" % (current_region, outfile), "w")
                current_file.write(line)

    def runBlast (infile, outfile):
        os.system ("blastn -query %s -db NT_033778.fna -out %s.output -outfmt 0" % (infile,outfile))

    #ebola_conservative_regions ('teflon_segs.txt', 'teflon_segment')

        no_segments =11
        for segment in range(no_segments):
        runBlast ('%d.teflon_segment' % segment, '%d.teflon_segment' % segment)
    '''
def main():
    s = siRNA()
    utr = Seq('AATGATTTGAAAGATCCGGGCCTAAAACCTTCATCA')
    seq_candidate = s.read_output('test2.fa')
    qualified_seqs = s.score_sys(seq_candidate, utr)
    print qualified_seqs
    shrnas = s.design_shrna(qualified_seqs)
    # check if the candidate sequence has the 2nd sequence
    s.check_2nd_structure(shrnas)
    # final designed shRNA
    final_seqs = s.determine_2nd_structure()
    # estimate the mutation generations
    mutation_v = s.mutation(qualified_seqs)
    s.my_print(mutation_v)


if __name__ == '__main__':
    main()
