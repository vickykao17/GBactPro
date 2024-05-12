from Bio import SeqIO
from collections import Counter
import random

# len_list = []
# input_file = "./scan_data_new_1/Escherichia_coli/negative_10.fasta"
# fasta_sequences = SeqIO.parse(open(input_file),'fasta')

# for fasta in fasta_sequences:
#     name, sequence = fasta.id, str(fasta.seq)
#     # print(len(sequence))
#     len_list.append(len(sequence))

# print(Counter(len_list))

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
            # 'Streptomyces coelicolor A3': 1016,
            # 'Zymomonas mobilis ZM4': 3205,
            # 'Anabaena sp. PCC7120': 3955,
            # 'Corynebacterium glutamicum ATCC 13032': 2147,
            # 'Corynebacterium diphtheriae NCTC 13129': 1202,
            # 'Mycobacterium tuberculosis H37Rv': 1778,
            # 'Escherichia coli MG1655': 1865,
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
            'Paenibacillus riograndensis SBR5': 1269,
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

regions = ["Minus10", "UP", "Minus35", "Spacer","DOWN", "35+s+10", "35+10"]
#            , "35+10/all_shuffle", "35+10/separate_shuffle", "35+s+10/all_shuffle", "35+s+10/separate_shuffle"

for b in bacteria.keys():
    for r in regions:
        print(b, r)
        len_list = []
        i = 0
        c = 0
        p_file = "../data/datasets/"+b+"/"+r+"/random/positive.fasta"
        all_path = "../data/datasets/"+b+"/"+r+"/random/all.fasta"
        all_file = open(all_path, 'w')
        fasta_sequences = SeqIO.parse(open(p_file),'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            # print(len(sequence))
            len_list.append(len(sequence))
            print('>'+name+'\n'+sequence, file=all_file)

        file_path = "../data/datasets/"+b+"/"+r+"/random/negative.fasta"
        n_path = "../data/datasets/"+b+"/"+r+"/random/negative_new.fasta"
        o_file = open(n_path, 'w')
        fasta_sequences = SeqIO.parse(open(file_path),'fasta')
        for fasta in fasta_sequences:
            if c == 10:
                c = 0
                i += 1
            length = len_list[i]
            name, sequence = fasta.id, str(fasta.seq)
            tmp = ''
            for s in sequence:
                if s not in ['A', 'T', "C", "G"]:
                    tmp += random.choice(["A", "T", "C", "G"])
                else:
                    tmp += s
            sequence = tmp
            print('>'+name+'\n'+sequence[:length], file=o_file)
            print('>'+name+'\n'+sequence[:length], file=all_file)
            c += 1
        o_file.close()
        all_file.close()
        

        


        