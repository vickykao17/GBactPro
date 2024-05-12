import pandas as pd
import random
import os
from Bio import SeqIO
import math
import argparse
from collections import Counter
from common import *

# name = "Escherichia coli MG1655"
# python generate_data/main.py -b Escherichia coli MG1655
# python main.py -b Escherichia coli MG1655
# python main.py -b Campylobacter jejuni NCTC11168


parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bacteria', nargs="+")
args = parser.parse_args()
name = ' '.join(args.bacteria)

regions = ["UP", "Minus35", "Spacer", "Minus10", "DOWN"]
# regions = ["DOWN"]
df_tss = pd.read_csv("./scan_result/" + name + ".tsv", sep="\t")

out_dir = "../data/datasets/" + name + "/"
if not os.path.exists(os.path.dirname(out_dir)):
    os.makedirs(os.path.dirname(out_dir))


if name == 'Escherichia coli MG1655':

    for r in regions:
        a = []
        for i in range(len(df_tss)):
            if i >= math.floor(len(df_tss) * 0.9):
                a.append('test')
            else:
                a.append('train')
        random.shuffle(a) 

        train_path = "../data/datasets/" + name + "/train/" + r + "/"
        test_path = "../data/datasets/" + name + "/test/" + r + "/"
        if not os.path.exists(os.path.dirname(train_path)):
            os.makedirs(os.path.dirname(train_path))
        if not os.path.exists(os.path.dirname(test_path)):
            os.makedirs(os.path.dirname(test_path))
        train_p_file = open(train_path + "positive.fasta", 'w')
        train_n_file = open(train_path + "negative_shuffle.fasta", 'w')
        test_p_file = open(test_path + "positive.fasta", 'w')
        test_n_file = open(test_path + "negative_shuffle.fasta", 'w')
        train_all = open(train_path + "all.fasta", 'w')
        test_all = open(test_path + "all.fasta", 'w')

        for i in range(len(df_tss)):
            if df_tss["dir"][i] == '+':
                loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m35_start"][i])+"-"+str(df_tss["loc_m10_end"][i])
            else:
                loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m10_end"][i])+"-"+str(df_tss["loc_m35_start"][i])
            
            # print(df_tss["Minus35"][i]+df_tss["Spacer"][i]+df_tss["Minus10"][i], file=p_file)
                
            if a[i] == 'train':
                print(loc + "\n" + df_tss[r][i], file=train_p_file)
                print(loc + "\n" + df_tss[r][i], file=train_all)
            else:
                print(loc + "\n" + df_tss[r][i], file=test_p_file)
                print(loc + "\n" + df_tss[r][i], file=test_all)
        for i in range(len(df_tss)):
            if df_tss["dir"][i] == '+':
                loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m35_start"][i])+"-"+str(df_tss["loc_m10_end"][i])
            else:
                loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m10_end"][i])+"-"+str(df_tss["loc_m35_start"][i])
            
            promoter = list(df_tss[r][i])
            for j in range(10):
                while random.shuffle(promoter) == list(df_tss[r][i]):
                    random.shuffle(promoter)
                if a[i] == 'train':
                    print(loc + "\n" + ''.join(promoter), file=train_n_file)
                    print(loc + "\n" + ''.join(promoter), file=train_all)
                else:
                    print(loc + "\n" + ''.join(promoter), file=test_n_file)
                    print(loc + "\n" + ''.join(promoter), file=test_all)


        train_p_file.close()
        train_n_file.close()
        test_p_file.close()
        test_n_file.close()
        train_all.close()
        test_all.close()

    train_path = "../data/datasets/" + name + "/train/35+s+10/"
    test_path = "../data/datasets/" + name + "/test/35+s+10/"
    if not os.path.exists(os.path.dirname(train_path)):
        os.makedirs(os.path.dirname(train_path))
    if not os.path.exists(os.path.dirname(test_path)):
        os.makedirs(os.path.dirname(test_path))
    train_p_file = open(train_path + "positive.fasta", 'w')
    train_n_file = open(train_path + "negative_shuffle.fasta", 'w')
    train_n_s_file = open(train_path + "negative_separate_shuffle.fasta", 'w')
    test_p_file = open(test_path + "positive.fasta", 'w')
    test_n_file = open(test_path + "negative_shuffle.fasta", 'w')
    test_n_s_file = open(test_path + "negative_separate_shuffle.fasta", 'w')
    train_all = open(train_path + "all.fasta", 'w')
    train_all2 = open(train_path + "all_s.fasta", 'w')
    test_all = open(test_path + "all.fasta", 'w')
    test_all2 = open(test_path + "all_s.fasta", 'w')
    # positive
    for i in range(len(df_tss)):
        if df_tss["dir"][i] == '+':
            loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m35_start"][i])+"-"+str(df_tss["loc_m10_end"][i])
        else:
            loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m10_end"][i])+"-"+str(df_tss["loc_m35_start"][i])

        if a[i] == 'train':
            print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Spacer"][i]+df_tss["Minus10"][i], file=train_p_file)
            print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Spacer"][i]+df_tss["Minus10"][i], file=train_all)
            print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Spacer"][i]+df_tss["Minus10"][i], file=train_all2)
        else:
            print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Spacer"][i]+df_tss["Minus10"][i], file=test_p_file)
            print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Spacer"][i]+df_tss["Minus10"][i], file=test_all)
            print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Spacer"][i]+df_tss["Minus10"][i], file=test_all2)

    # negative
    for i in range(len(df_tss)):
        if df_tss["dir"][i] == '+':
            loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m35_start"][i])+"-"+str(df_tss["loc_m10_end"][i])
        else:
            loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m10_end"][i])+"-"+str(df_tss["loc_m35_start"][i])
        # all shuffle
        promoter = list(df_tss["Minus35"][i]+df_tss["Spacer"][i]+df_tss["Minus10"][i])
        for j in range(10):
            while random.shuffle(promoter) == list(df_tss[r][i]):
                random.shuffle(promoter)
            if a[i] == 'train':
                print(loc + "\n" + ''.join(promoter), file=train_n_file)
                print(loc + "\n" + ''.join(promoter), file=train_all)
            else:
                print(loc + "\n" + ''.join(promoter), file=test_n_file)
                print(loc + "\n" + ''.join(promoter), file=test_all)
        # separate shuffle
        _10 = list(df_tss["Minus10"][i])
        spacer = list(df_tss["Spacer"][i])
        _35 = list(df_tss["Minus35"][i])
        for j in range(10):
            while random.shuffle(_10) == list(df_tss["Minus10"][i]):
                random.shuffle(_10)
            while random.shuffle(spacer) == list(df_tss["Spacer"][i]):
                random.shuffle(spacer)
            while random.shuffle(_35) == list(df_tss["Minus35"][i]):
                random.shuffle(_35)


            promoter = ''.join(_35) + ''.join(spacer) + ''.join(_10)
            if a[i] == 'train':
                print(loc + "\n" + promoter, file=train_n_s_file)
                print(loc + "\n" + promoter, file=train_all2)
            else:
                print(loc + "\n" + promoter, file=test_n_s_file)
                print(loc + "\n" + promoter, file=test_all2)

    train_p_file.close()
    train_n_file.close()
    train_n_s_file.close()
    test_p_file.close()
    test_n_file.close()
    test_n_s_file.close()
    train_all.close()
    train_all2.close()
    test_all.close()
    test_all2.close()
    
    train_path = "../data/datasets/" + name + "/train/35+10/"
    test_path = "../data/datasets/" + name + "/test/35+10/"
    if not os.path.exists(os.path.dirname(train_path)):
        os.makedirs(os.path.dirname(train_path))
    if not os.path.exists(os.path.dirname(test_path)):
        os.makedirs(os.path.dirname(test_path))
    train_p_file = open(train_path + "positive.fasta", 'w')
    train_n_file = open(train_path + "negative_shuffle.fasta", 'w')
    train_n_s_file = open(train_path + "negative_separate_shuffle.fasta", 'w')
    test_p_file = open(test_path + "positive.fasta", 'w')
    test_n_file = open(test_path + "negative_shuffle.fasta", 'w')
    test_n_s_file = open(test_path + "negative_separate_shuffle.fasta", 'w')
    train_all = open(train_path + "all.fasta", 'w')
    train_all2 = open(train_path + "all_s.fasta", 'w')
    test_all = open(test_path + "all.fasta", 'w')
    test_all2 = open(test_path + "all_s.fasta", 'w')
    # positive
    for i in range(len(df_tss)):
        if df_tss["dir"][i] == '+':
            loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m35_start"][i])+"-"+str(df_tss["loc_m10_end"][i])
        else:
            loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m10_end"][i])+"-"+str(df_tss["loc_m35_start"][i])

        if a[i] == 'train':
            print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Minus10"][i], file=train_p_file)
            print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Minus10"][i], file=train_all)
            print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Minus10"][i], file=train_all2)
        else:
            print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Minus10"][i], file=test_p_file)
            print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Minus10"][i], file=test_all)
            print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Minus10"][i], file=test_all2)

    # negative
    for i in range(len(df_tss)):
        if df_tss["dir"][i] == '+':
            loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m35_start"][i])+"-"+str(df_tss["loc_m10_end"][i])
        else:
            loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m10_end"][i])+"-"+str(df_tss["loc_m35_start"][i])
        # all shuffle
        promoter = list(df_tss["Minus35"][i]+df_tss["Minus10"][i])
        for j in range(10):
            while random.shuffle(promoter) == list(df_tss[r][i]):
                random.shuffle(promoter)
            if a[i] == 'train':
                print(loc + "\n" + ''.join(promoter), file=train_n_file)
                print(loc + "\n" + ''.join(promoter), file=train_all)
            else:
                print(loc + "\n" + ''.join(promoter), file=test_n_file)
                print(loc + "\n" + ''.join(promoter), file=test_all)
        # separate shuffle
        _10 = list(df_tss["Minus10"][i])
        _35 = list(df_tss["Minus35"][i])
        for j in range(10):
            while random.shuffle(_10) == list(df_tss["Minus10"][i]):
                random.shuffle(_10)
            while random.shuffle(_35) == list(df_tss["Minus35"][i]):
                random.shuffle(_35)

            promoter = ''.join(_35) + ''.join(_10)
            if a[i] == 'train':
                print(loc + "\n" + promoter, file=train_n_s_file)
                print(loc + "\n" + promoter, file=train_all2)
            else:
                print(loc + "\n" + promoter, file=test_n_s_file)
                print(loc + "\n" + promoter, file=test_all2)
                
    train_p_file.close()
    train_n_file.close()
    train_n_s_file.close()
    test_p_file.close()
    test_n_file.close()
    test_n_s_file.close()
    train_all.close()
    train_all2.close()
    test_all.close()
    test_all2.close()
else:
    for r in regions:
        
        path = "./dataset/" + name + "/" + r + "/"
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))
        p_file = open(path + "positive.fasta", 'w')
        g = df_tss["genome"].iloc[0]
        arr = []
        for i in range(len(df_tss)):
            if df_tss["dir"][i] == '+':
                loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m35_start"][i])+"-"+str(df_tss["loc_m10_end"][i])
            else:
                loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m10_end"][i])+"-"+str(df_tss["loc_m35_start"][i])
            
            print(loc + "\n" + df_tss[r][i], file=p_file)

        g = df_tss["genome"].iloc[0]
        file_g = "./genome/" + g + ".gb"
        genome = Genome(next(SeqIO.parse(open(file_g, "r"), "genbank")))
        print(genome.sequence)
            
        #     arr.append([g, int(df_tss["loc"].iloc[i])-1, df_tss["loc"].iloc[i]-1+len(df_tss[r][i]), df_tss["gene"].iloc[i], 0, df_tss["dir"].iloc[i]])
        
        # df = pd.DataFrame(arr)
        # df.to_csv("scan_data/"+name+"/"+r+"/positive.bed", sep = "\t", header=0, index=0)
        p_file.close()


    path = "./dataset/" + name + "/35+s+10/"
    if not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))
    p_file = open(path + "positive.fasta", 'w')


    # positive
    for i in range(len(df_tss)):
        # loc = ">"+df_tss["genome"][i]
        if df_tss["dir"][i] == '+':
            loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m35_start"][i])+"-"+str(df_tss["loc_m10_end"][i])
        else:
            loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m10_end"][i])+"-"+str(df_tss["loc_m35_start"][i])
        print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Spacer"][i]+df_tss["Minus10"][i], file=p_file)

        arr.append([g, int(df_tss["loc_m35_start"][i])-1, df_tss["loc"].iloc[i]-1+len(df_tss[r][i]), df_tss["gene"].iloc[i], 0, df_tss["dir"].iloc[i]])
        
    df = pd.DataFrame(arr)
    df.to_csv("scan_data/"+name+"/35+s+10/positive.bed", sep = "\t", header=0, index=0)

        

    path = "./dataset/" + name + "/35+10/"
    if not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))
    p_file = open(path + "positive.fasta", 'w')

    # positive
    for i in range(len(df_tss)):
        # loc = ">"+df_tss["genome"][i]
        if df_tss["dir"][i] == '+':
            loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m35_start"][i])+"-"+str(df_tss["loc_m10_end"][i])
        else:
            loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m10_end"][i])+"-"+str(df_tss["loc_m35_start"][i])
        print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Minus10"][i], file=p_file)

        arr.append([g, int(df_tss["loc"].iloc[i])-1, df_tss["loc"].iloc[i]-1+len(df_tss[r][i]), df_tss["gene"].iloc[i], 0, df_tss["dir"].iloc[i]])
        
    df = pd.DataFrame(arr)
    df.to_csv("scan_data/"+name+"/35+10/positive.bed", sep = "\t", header=0, index=0)


