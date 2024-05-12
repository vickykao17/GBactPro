import random
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Bio.GenBank.Record
from Bio import SeqIO
from Bio import motifs

from BPM.BPM import *


# class Fasta:

#     def __init__(self, record: Bio.SeqRecord.SeqRecord):
#         self.record = record
#         self.length = len(record.seq)

#     def __repr__(self):
#         return str((self.record.id, self.record.seq))

#     def __str__(self):
#         return str((self.record.id, self.record.seq))

#     def sequence(self, l, r):
        
#         if r > l and l >= 0 and r <= self.length:
#             ret = self.record.seq[l:r]
#         else:
#             raise Exception("Invalid index range for linear sequence")

#         return ret

#     def background(self):
#         ## Too slow
#         a = self.record.seq.count("A")/self.length
#         t = self.record.seq.count("T")/self.length
#         c = self.record.seq.count("C")/self.length
#         g = self.record.seq.count("G")/self.length

#         '''
#         ## Probabilistic
#         n = 10000
#         s = random.choices(str(self.record.seq), k=n)
#         a = s.count("A")
#         t = s.count("T")
#         c = s.count("C")
#         g = s.count("G")

#         total = a + t + c + g
#         at = (a + t) / (2 * total)
#         gc = (g + c) / (2 * total)
#         '''

#         return {'A': a, 'T': t, 'C': c, 'G': g}

'''
class Locus:

    def __init__(self, l: int, r: int, f: Fasta):
        self.fasta = f
        self.l = l
        self.r = r
        self._check_idx()

    def __repr__(self):
        return str((self.l, self.r, self.getSequence()))

    def __str__(self):
        return str((self.l, self.r, self.getSequence()))

    def __add__(self, offset: int):
        self.l += offset
        self.r += offset
        self._check_idx()
        return self

    def __radd__(self, offset: int):
        return self.__add__(offset)

    def _check_idx(self):
        if self.r <= self.l or self.l < 0 or self.r > self.fasta.length:
            raise Exception("Invalid index range for linear sequence")

    def getLoc(self):
        return (self.l, self.r)

    def setLoc(self, loc):
        self.l, self.r = loc
        self._check_idx()

    def getSequence(self):
        return self.fasta.sequence(self.l, self.r)

    def length(self):
        return len(self.getSequence())

    def toTuple(self):
        seq = str(self.getSequence())
        return (self.l, self.r, seq[:6], seq[6:-6], seq[-6:])
'''


class Genome:

    def __init__(self, record: Bio.GenBank.Record.Record):
        self.record = record
        self.length = len(record.seq)

        t = record.annotations["topology"]

        if t in ("circular", "linear"):
            self.type = t
        else:
            raise Exception(f"Invalid sequence topology {t}")

    def sequence(self, l, r):
        l = l % self.length
        r = r % self.length
        ret = None

        if r > l:
            ret = self.record.seq[l : r] 
        elif r <= l:
            ret = self.record.seq[l : ] + self.record.seq[ : r]

            if r != 0 and self.type == "linear":
                print(self.type)
                raise Exception("Invalid index range for linear sequence")

        return ret

    def background(self):
        # Too slow
        #a = self.record.seq.count("A")
        #t = self.record.seq.count("T")
        #c = self.record.seq.count("C")
        #g = self.record.seq.count("G")
        
        # Probabilistic
        n = 10000
        s = random.choices(str(self.record.seq), k=n)
        a = s.count("A")
        t = s.count("T")
        c = s.count("C")
        g = s.count("G")

        total = a + t + c + g
        at = (a + t) / (2 * total)
        gc = (g + c) / (2 * total)
        return {'A': at, 'T': at, 'C': gc, 'G': gc}

    def __repr__(self):
        return self.record.id + "\t" + self.record.description

    def __str__(self):
        return self.record.id + "\t" + self.record.description
        
    def __copy__(self):
        return Genome(self.record)


class Locus:

    def __init__(self, l: int, r: int, d: str, g: Genome):
        self.genome = g
        self.d = d
        self.l = l
        self.r = r
        self._stdidx()

    def __repr__(self):
        return str((self.l, self.r, self.d, self.getSequence()))

    def __str__(self):
        return str((self.l, self.r, self.d, self.getSequence()))

    def __add__(self, offset):
        self.l += offset
        self.r += offset
        self._stdidx()
        return self

    def __radd__(self, offset):
        return self.__add__(offset)

    ## Standardize index to the range of genome
    def _stdidx(self):

        while self.r < 0:
            self.r += self.genome.length

        while self.l < 0:
            self.l += self.genome.length

        self.l = self.l % self.genome.length
        self.r = self.r % self.genome.length

        if self.r != 0 and self.r < self.l and self.genome.type == "linear":
            raise Exception("Invalid index range for linear sequence")

    def getLoc(self):
        return (self.l, self.r)

    def getLoc_m10_end(self):
        if self.d == '+':
            return self.r
        else:
            return self.l
        
    def getLoc_m35_start(self):
        if self.d == '+':
            return self.l
        else:
            return self.r

    def getDir(self):
        return self.d

    def setLoc(self, loc):
        self.l, self.r = loc
        self._stdidx()

    def getSequence(self):
        seq = self.genome.sequence(self.l, self.r)
        if self.d == "-":
            seq = seq.reverse_complement()
        return seq

    def length(self):
        return len(self.getSequence())

    def toTuple(self):
        return (self.l, self.r, self.d)
    
#     def subLocus(self, l, r):
#         if self.d == "+":
#             return Locus(self.l + l, self.l + r - l, self.d, self.genome)
#         else:
#             return Locus(self.r - r - 1, self.r - l, self.d, self.genome)
            
#     def subLocusIterator(self, length):
#         return LocusIterator(self.l, self.r, length, self.genome, direction=self.d)




## Conversion from 1-based coordinate(NCBI) to 0-based coordinate(throughout the code)
def b1to0(l, r):
    return l - 1, r

def b0to1(l, r):
    return l + 1, r

def shuffle_sequence(x):
    x = list(x)
    random.shuffle(x)
    return "".join(x)


"""
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Promoter annotator v0.0")
    parser.add_argument("file", help="fasta file")
    
    args = parser.parse_args()
    fasta = Fasta(next(SeqIO.parse(open(args.file, "r"), "fasta")))
    fasta_name = args.file.split('/')[-1].split('.')[0]
    print(fasta_name)

    promoters = []
    expressions = [] # log10-scale

    for i in range(fasta.length-30): # 30 = 12 + 19 - 1
        ps = [Locus(0+i+j, 31+i, fasta) for j in range(5)]
        ps = [(p, predict(p.getSequence())) for p in ps]
        p = max(ps, key=lambda x: x[1])
        promoters.append(p[0])
        expressions.append(p[1])
        
    
    ## Plot ====
    exp = np.concatenate([np.repeat(np.nan, 25), 10**np.array(expressions), np.repeat(np.nan, 5)]) # 25 = 6 + 19
    num_subplots = fasta.length//500 + 1 # Set the number of subplots based on the length of fasta sequence
    fig, axes = plt.subplots(num_subplots, 1, figsize=(40, 2*num_subplots), dpi=100) # Create subplots

    for i, ax in enumerate(axes): # Plot each subplot
        start_index = i * 500
        end_index = min((i + 1) * 500, fasta.length)

        ax.plot(np.arange(start_index, end_index), exp[start_index:end_index])

        ax.set_xlim(start_index - 0.5, start_index + 500.5)
        ax.set_ylim(0 - np.nanmax(exp)/20, np.nanmax(exp)*1.05)

        ax.set_xticks(np.arange(start_index, end_index))
        ax.set_xticklabels(list(fasta.record.seq[start_index:end_index]), fontsize=8)
        ax.set_ylabel('Promoter strength (a.u.)')

    plt.tight_layout() # Adjust layout to prevent clipping of xlabel
    plt.savefig(f'output/{fasta_name}.png') 
    plt.show()
    
    ## Output promoter information
    idx = np.argsort(expressions)[::-1][:3] # top 3
    df = pd.DataFrame([promoters[i].toTuple() for i in idx], columns=['Start', 'End', 'Minus35', 'Spacer', 'Minus10'])
    df['Promoter strength'] = [10**expressions[i] for i in idx]
    df.to_csv(f'output/{fasta_name}.tsv', sep='\t', index=False)
"""
