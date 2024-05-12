import pandas as pd
import random
import os
from Bio import SeqIO
import math
import argparse
from collections import Counter

# name = "Escherichia coli MG1655"
# python generate_data/main.py -b Escherichia coli MG1655
# python main.py -b Escherichia coli MG1655
# python main.py -b Campylobacter jejuni NCTC11168



# parser = argparse.ArgumentParser()
# parser.add_argument('-b', '--bacteria', nargs="+")
# args = parser.parse_args()
# name = ' '.join(args.bacteria)

bacteria = {
            # 'Synechocystis sp. PCC 6803': 170,
            # 'Pseudomonas aeruginosa UCBPP-PA14': 2117,
            # 'Nanosynbacter lyticus TM7x': 188,
            # 'Mycolicibacterium smegmatis MC2 155': 3043, 
            # 'Chlamydia pneumoniae CWL029': 382,
            # 'Chlamydia trachomatis L2b': 349,
            # 'Thermotoga maritima MSB8': 550,
            # 'Burkholderia cenocepacia J2315': 1086,
            # 'Clostridioides difficile 630': 1288,
            # 'Mycoplasma pneumoniae M129': 718,
            # 'Neisseria gonorrhoeae MS11': 919,
            # 'Bordetella pertussis Tohama I': 593,
            # 'Caulobacter crescentus NA1000': 1443,
            # 'Leptospira interrogans L495': 2864,
            # 'Lachnoclostridium phytofermentans ISDg': 2057,
            # 'Streptomyces coelicolor A3(2)': 1016,
            # 'Zymomonas mobilis ZM4': 3205,
            # 'Anabaena sp. PCC7120': 3955,
            # 'Corynebacterium glutamicum ATCC 13032': 2147,
            # 'Corynebacterium diphtheriae NCTC 13129': 1202,
            # 'Mycobacterium tuberculosis H37Rv': 1778,
            'Escherichia coli MG1655': 1865,
            # 'Bacillus subtilis 168': 600,
            # 'Bacteroides thetaiotaomicron VPI-5482': 1616,
            # 'Salmonella typhimurium SL1344': 1120,
            # 'Synechococcus elongatus UTEX 2973': 2429,
            # 'Xanthomonas campestris B100': 1545,
            # 'Helicobacter pylori 26695': 716,
            # 'Staphylococcus aureus MW2': 1800,
            # 'Acinetobacter baumannii ATCC 17978':864,
            # 'Staphylococcus epidermidis ATCC 12228': 1206,
            # 'Klebsiella aerogenes KCTC 2190': 468,
            # 'Paenibacillus riograndensis SBR5': 1269,
            # 'Methylorubrum extorquens DM4': 2030,
            # 'Shewanella oneidensis MR-1': 1802,
            # 'Agrobacterium tumefaciens C58': 468,
            # 'Campylobacter jejuni NCTC11168': 675,
            # 'Bradyrhizobium japonicum USDA 110': 4451,
            # 'Listeria monocytogenes': 1079,
            # 'Listeria innocua': 911,
            # 'Streptococcus pyrogenes S119': 683,
            # 'Streptococcus agalactiae NEM316': 885,
            # 'Vibrio cholerae N16961': 1283,
            # 'Shigella flexneri M90T': 2051,
            # 'Pseudomonas putida KT2440': 1142,
            # 'Phytoplasma sp. OY-M': 82,
            # 'Geobacter sulfurreducens PCA': 875,
            # 'Novosphingobium aromaticivorans DSM 12444': 944,
            # 'Rhodobacter sphaeroides 2.4.1': 1087
            }

regions = ["UP", "Minus35", "Spacer", "Minus10", "DOWN"]
# regions = ["DOWN"]
for b in bacteria.keys():  
    name = b  
    df_tss = pd.read_csv("../data/scan_result/" + name + ".tsv", sep="\t")

    out_dir = "../data/datasets/" + name + "/"
    if not os.path.exists(os.path.dirname(out_dir)):
        os.makedirs(os.path.dirname(out_dir))

    
    if name != 'Escherichia coli MG1655':

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
            path = "../data/datasets/" + name + "/" + r + "/random/"
            if not os.path.exists(os.path.dirname(path)):
                os.makedirs(os.path.dirname(path))
            p_file = open(path + "positive.fasta", 'w')
            # n_file = open(path + "negative.fasta", 'w')
            # all = open(path + "all.fasta", 'w')
            for i in range(len(df_tss)):
                if 'N' not in df_tss["UP"][i] and 'N' not in df_tss["DOWN"][i]:
                    if df_tss["dir"][i] == '+':
                        loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m35_start"][i])+"-"+str(df_tss["loc_m10_end"][i])
                    else:
                        loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m10_end"][i])+"-"+str(df_tss["loc_m35_start"][i])
                    
                    # print(df_tss["Minus35"][i]+df_tss["Spacer"][i]+df_tss["Minus10"][i], file=p_file)
                    print(loc + "\n" + df_tss[r][i], file=p_file)
                    # print(loc + "\n" + df_tss[r][i], file=all)
            # for i in range(len(df_tss)):
            #     if 'N' not in df_tss["UP"][i] and 'N' not in df_tss["DOWN"][i]:
            #         if df_tss["dir"][i] == '+':
            #             loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m35_start"][i])+"-"+str(df_tss["loc_m10_end"][i])
            #         else:
            #             loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m10_end"][i])+"-"+str(df_tss["loc_m35_start"][i])
                    
            #         promoter = list(df_tss[r][i])
            #         for j in range(10):
            #             while random.shuffle(promoter) == list(df_tss[r][i]):
            #                 random.shuffle(promoter)
            #             print(loc + "\n" + ''.join(promoter), file = n_file)
            #             print(loc + "\n" + ''.join(promoter), file = all)

            p_file.close()
            # n_file.close()
            # all.close()

        path = "../data/datasets/" + name + "/35+s+10/"
        # if not os.path.exists(os.path.dirname(path+"all_shuffle/")):
        #     os.makedirs(os.path.dirname(path+"all_shuffle/"))
        # if not os.path.exists(os.path.dirname(path+"separate_shuffle/")):
        #     os.makedirs(os.path.dirname(path+"separate_shuffle/"))

        if not os.path.exists(os.path.dirname(path+"random/")):
            os.makedirs(os.path.dirname(path+"random/"))
        p_file = open(path + "random/positive.fasta", 'w')
        # p_file = open(path + "all_shuffle/positive.fasta", 'w')
        # n_file = open(path + "all_shuffle/negative.fasta", 'w')
        # p_s_file = open(path + "separate_shuffle/positive.fasta", 'w')
        # n_s_file = open(path + "separate_shuffle/negative.fasta", 'w')
        # all = open(path + "all_shuffle/all.fasta", 'w')
        # all2 = open(path + "separate_shuffle/all.fasta", 'w')

        # positive
        for i in range(len(df_tss)):
            if 'N' not in df_tss["UP"][i] and 'N' not in df_tss["DOWN"][i]:
                if df_tss["dir"][i] == '+':
                    loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m35_start"][i])+"-"+str(df_tss["loc_m10_end"][i])
                else:
                    loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m10_end"][i])+"-"+str(df_tss["loc_m35_start"][i])
                print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Spacer"][i]+df_tss["Minus10"][i], file=p_file)
                # print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Spacer"][i]+df_tss["Minus10"][i], file=p_s_file)
                # print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Spacer"][i]+df_tss["Minus10"][i], file=all)
                # print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Spacer"][i]+df_tss["Minus10"][i], file=all2)
        # negative
        # for i in range(len(df_tss)):
        #     if 'N' not in df_tss["UP"][i] and 'N' not in df_tss["DOWN"][i]:
        #         if df_tss["dir"][i] == '+':
        #             loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m35_start"][i])+"-"+str(df_tss["loc_m10_end"][i])
        #         else:
        #             loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m10_end"][i])+"-"+str(df_tss["loc_m35_start"][i])
        #         promoter = list(df_tss["Minus35"][i]+df_tss["Spacer"][i]+df_tss["Minus10"][i])
        #         for j in range(10):
        #             while random.shuffle(promoter) == list(df_tss[r][i]):
        #                 random.shuffle(promoter)
        #             print(loc + "\n" + ''.join(promoter), file = n_file)
        #             print(loc + "\n" + ''.join(promoter), file = all)
                
        #         _10 = list(df_tss["Minus10"][i])
        #         spacer = list(df_tss["Spacer"][i])
        #         _35 = list(df_tss["Minus35"][i])
        #         for j in range(10):
        #             while random.shuffle(_10) == list(df_tss["Minus10"][i]):
        #                 random.shuffle(_10)
        #             while random.shuffle(spacer) == list(df_tss["Spacer"][i]):
        #                 random.shuffle(spacer)
        #             while random.shuffle(_35) == list(df_tss["Minus35"][i]):
        #                 random.shuffle(_35)

        #             promoter = ''.join(_35) + ''.join(spacer) + ''.join(_10)
        #             print(loc + "\n" + promoter, file = n_s_file)
        #             print(loc + "\n" + promoter, file = all2)

        path = "../data/datasets/" + name + "/35+10/"
        # if not os.path.exists(os.path.dirname(path+"all_shuffle/")):
        #     os.makedirs(os.path.dirname(path+"all_shuffle/"))
        # if not os.path.exists(os.path.dirname(path+"separate_shuffle/")):
        #     os.makedirs(os.path.dirname(path+"separate_shuffle/"))

        if not os.path.exists(os.path.dirname(path+"random/")):
            os.makedirs(os.path.dirname(path+"random/"))
        p_file = open(path + "random/positive.fasta", 'w')
        # p_file = open(path + "all_shuffle/positive.fasta", 'w')
        # n_file = open(path + "all_shuffle/negative.fasta", 'w')
        # p_s_file = open(path + "separate_shuffle/positive.fasta", 'w')
        # n_s_file = open(path + "separate_shuffle/negative.fasta", 'w')
        # all = open(path + "all_shuffle/all.fasta", 'w')
        # all2 = open(path + "separate_shuffle/all.fasta", 'w')
        # positive
        for i in range(len(df_tss)):
            if 'N' not in df_tss["UP"][i] and 'N' not in df_tss["DOWN"][i]:
                if df_tss["dir"][i] == '+':
                    loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m35_start"][i])+"-"+str(df_tss["loc_m10_end"][i])
                else:
                    loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m10_end"][i])+"-"+str(df_tss["loc_m35_start"][i])
                print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Minus10"][i], file=p_file)
                # print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Minus10"][i], file=p_s_file)
                # print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Minus10"][i], file=all)
                # print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Minus10"][i], file=all2)
        # negative
        # for i in range(len(df_tss)):
        #     if 'N' not in df_tss["UP"][i] and 'N' not in df_tss["DOWN"][i]:
        #         if df_tss["dir"][i] == '+':
        #             loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m35_start"][i])+"-"+str(df_tss["loc_m10_end"][i])
        #         else:
        #             loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m10_end"][i])+"-"+str(df_tss["loc_m35_start"][i])
        #         promoter = list(df_tss["Minus35"][i]+df_tss["Minus10"][i])
        #         for j in range(10):
        #             while random.shuffle(promoter) == list(df_tss["Minus35"][i]+df_tss["Minus10"][i]):
        #                 random.shuffle(promoter)
        #             print(loc + "\n" + ''.join(promoter), file = n_file)
        #             print(loc + "\n" + ''.join(promoter), file = all)
        #         _10 = list(df_tss["Minus10"][i])
        #         _35 = list(df_tss["Minus35"][i])
        #         for j in range(10):
        #             while random.shuffle(_10) == list(df_tss["Minus10"][i]):
        #                 random.shuffle(_10)
        #             while random.shuffle(_35) == list(df_tss["Minus35"][i]):
        #                 random.shuffle(_35)

        #             promoter = ''.join(_35) + ''.join(_10)
        #             print(loc + "\n" + promoter, file = n_s_file)
        #             print(loc + "\n" + promoter, file = all2)
