import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from time import time
from tqdm import tqdm
from copy import deepcopy
from Bio import SeqIO
from Bio.Seq import MutableSeq, Seq

from common import *


### Functions ====
def getPromoters(tss_loc, direction, genome, tss_m10_dist, promoter_lens=[12+15, 12+16, 12+17, 12+18, 12+19]):
    promoters = []

    if direction == "+":
        p = tss_loc - tss_m10_dist - 1
        for pl in promoter_lens:
            l, r = b1to0(p - pl + 1, p) 
            promoters.append(Locus(l, r, direction, genome))

    else:
        p = tss_loc + tss_m10_dist + 1
        for pl in promoter_lens:
            l, r = b1to0(p, p + pl - 1)
            promoters.append(Locus(l+1, r+1, direction, genome))

    return promoters


def scan_promoter(tss_loc, direction, genome, tss_m10_range):

    ## enumerate all possible promoter sites
    ps = []
    for tss_m10_dist in range(tss_m10_range[0], tss_m10_range[1]+1):
        ps += getPromoters(tss_loc, direction, genome, tss_m10_dist)

    ## score these all possible promoters
    ps = [(p, score_promoter(p.getSequence())) for p in ps]

    ## find the highest RNAP binding site
    promoter = min(ps, key=lambda x: x[1])

    return (promoter[0], (tss_loc, direction))


def scan_promoter_shuffle(tss_loc, direction, genome_copied, tss_m10_range):

    ## Determine the promoter region
    if direction == "+":
        r_lim = tss_loc - tss_m10_range[0]
        l_lim = tss_loc - (12+19) - tss_m10_range[1]
    else:
        r_lim = tss_loc + (12+19) + tss_m10_range[1]
        l_lim = tss_loc + tss_m10_range[0] + 1

    l_lim, r_lim = b1to0(l_lim, r_lim)
    p_region = Locus(l_lim, r_lim, direction, genome_copied)
    
    ## Shuffle the sequence of the promoter region
    s = MutableSeq(genome_copied.record.seq)
    s[p_region.l:p_region.r] = shuffle_sequence(str(genome_copied.record.seq[p_region.l:p_region.r]))
    genome_copied.record.seq = Seq(s)
    
    return scan_promoter(tss_loc, direction, genome_copied, tss_m10_range)

def get_element_UP(locus, length=100):
    if locus.d == '+':
        up = locus.genome.sequence(locus.l-length, locus.l)
    else:
        up = locus.genome.sequence(locus.r, locus.r+length).reverse_complement()
    return up

def get_element_DOWN(locus, length=60):
    if locus.d == '+':
        down = locus.genome.sequence(locus.r, locus.r+length)
    else:
        down = locus.genome.sequence(locus.l-length, locus.l).reverse_complement()
    return down


'''
### Main ====
if __name__ == '__main__':

    t_start = time()

    tss_m10_range = [-10, 20]
    print(tss_m10_range)
    
    TSS = pd.read_csv("../data/tss_20231129.tsv", sep="\t")

    for _, df_tss in TSS.groupby("name"):
        
        ## import genome ====
        g = df_tss["genome"].iloc[0]
        name = df_tss["name"].iloc[0]
        file_g = "../data/genome/" + g + ".gb"
        genome = Genome(next(SeqIO.parse(open(file_g, "r"), "genbank")))
        genome_copied = deepcopy(genome) # for tss_shuffle
        tax = genome.record.annotations["taxonomy"]
        
        print(name)
        
        ## -35/-10 scanning ====
        tss = []
        tss_random = []
        tss_shuffle = []

        for idx, row in df_tss.iterrows():
            direction = row["dir"]
            tss_loc = int(row["loc"])

            direction_random = random.choice(['+', '-'])
            tss_loc_random = random.choice(range(0, genome.length)) + 1

            tss.append(scan_promoter(tss_loc, direction, genome, tss_m10_range))
            tss_random.append(scan_promoter(tss_loc_random, direction_random, genome, tss_m10_range))
            tss_shuffle.append(scan_promoter_shuffle(tss_loc, direction, genome_copied, tss_m10_range))

        ## output tables ====
        df_tss["UP"] = [str(get_element_UP(t[0])) for t in tss]
        df_tss["Minus35"] = [str(t[0].getSequence()[:6]) for t in tss]
        df_tss["Spacer"] = [str(t[0].getSequence()[6:-6]) for t in tss]
        df_tss["Minus10"] = [str(t[0].getSequence()[-6:]) for t in tss]
        df_tss["DOWN"] = [str(get_element_DOWN(t[0])) for t in tss]
        df_tss.to_csv(f'../table/{name}.tsv', sep='\t', index=False)

        df_tss_random = pd.DataFrame({
            "loc": [t[1][0] for t in tss_random], 
            "dir": [t[1][1] for t in tss_random], 
            "UP": [str(get_element_UP(t[0])) for t in tss_random], 
            "Minus35": [str(t[0].getSequence()[:6]) for t in tss_random], 
            "Spacer": [str(t[0].getSequence()[6:-6]) for t in tss_random], 
            "Minus10": [str(t[0].getSequence()[-6:]) for t in tss_random], 
            "DOWN": [str(get_element_DOWN(t[0])) for t in tss_random], 
            })
        df_tss_random.to_csv(f'../table/{name}_random.tsv', sep='\t', index=False)

        df_tss_shuffle = deepcopy(df_tss)
        df_tss_shuffle["UP"] = [str(get_element_UP(t[0])) for t in tss_shuffle]
        df_tss_shuffle["Minus35"] = [str(t[0].getSequence()[:6]) for t in tss_shuffle]
        df_tss_shuffle["Spacer"] = [str(t[0].getSequence()[6:-6]) for t in tss_shuffle]
        df_tss_shuffle["Minus10"] = [str(t[0].getSequence()[-6:]) for t in tss_shuffle]
        df_tss_shuffle["DOWN"] = [str(get_element_DOWN(t[0])) for t in tss_shuffle]
        df_tss_shuffle.to_csv(f'../table/{name}_shuffle.tsv', sep='\t', index=False)

    t_end = time()
    print(f"time usage: {int(t_end-t_start)} seconds")
'''