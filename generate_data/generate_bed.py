import pandas as pd
from Bio import SeqIO
from common import *

bacteria = {
            'Synechocystis sp. PCC 6803': 170,
            'Pseudomonas aeruginosa UCBPP-PA14': 2117,
            'Nanosynbacter lyticus TM7x': 188,
            'Mycolicibacterium smegmatis MC2 155': 3043, 
            'Chlamydia pneumoniae CWL029': 382,
            'Chlamydia trachomatis L2b': 349,
            'Thermotoga maritima MSB8': 550,
            'Burkholderia cenocepacia J2315': 1086,
            'Clostridioides difficile 630': 1288,
            'Mycoplasma pneumoniae M129': 718,
            'Neisseria gonorrhoeae MS11': 919,
            'Bordetella pertussis Tohama I': 593,
            'Caulobacter crescentus NA1000': 1443,
            'Leptospira interrogans L495': 2864,
            'Lachnoclostridium phytofermentans ISDg': 2057,
            'Streptomyces coelicolor A3(2)': 1016,
            'Zymomonas mobilis ZM4': 3205,
            'Anabaena sp. PCC7120': 3955,
            'Corynebacterium glutamicum ATCC 13032': 2147,
            'Corynebacterium diphtheriae NCTC 13129': 1202,
            'Mycobacterium tuberculosis H37Rv': 1778,
            'Escherichia coli MG1655': 1865,
            'Bacillus subtilis 168': 600,
            'Bacteroides thetaiotaomicron VPI-5482': 1616,
            'Salmonella typhimurium SL1344': 1120,
            'Synechococcus elongatus UTEX 2973': 2429,
            'Xanthomonas campestris B100': 1545,
            'Helicobacter pylori 26695': 716,
            'Staphylococcus aureus MW2': 1800,
            'Acinetobacter baumannii ATCC 17978':864,
            'Staphylococcus epidermidis ATCC 12228': 1206,
            'Klebsiella aerogenes KCTC 2190': 468,
            'Paenibacillus riograndensis SBR5': 1269,
            'Methylorubrum extorquens DM4': 2030,
            'Shewanella oneidensis MR-1': 1802,
            'Agrobacterium tumefaciens C58': 468,
            'Campylobacter jejuni NCTC11168': 675,
            'Bradyrhizobium japonicum USDA 110': 4451,
            'Listeria monocytogenes': 1079,
            'Listeria innocua': 911,
            'Streptococcus pyrogenes S119': 683,
            'Streptococcus agalactiae NEM316': 885,
            'Vibrio cholerae N16961': 1283,
            'Shigella flexneri M90T': 2051,
            'Pseudomonas putida KT2440': 1142,
            'Phytoplasma sp. OY-M': 82,
            'Geobacter sulfurreducens PCA': 875,
            'Novosphingobium aromaticivorans DSM 12444': 944,
            'Rhodobacter sphaeroides 2.4.1': 1087
            }

# regions = ["UP", "Minus35", "Spacer", "Minus10", "DOWN", "35+10/separate_shuffle", "35+s+10/all_shuffle"]


for b in bacteria.keys():
    name = b
    path = "../data/scan_result/" + name
    df_tss = pd.read_csv(path + ".tsv", sep = '\t')

    arr = []
    g = df_tss["genome"].iloc[0]
    # for i in range(len(df_tss)):
        
    #     # # loc = ">"+df_tss["genome"][i]
    #     # if df_tss["dir"][i] == '+':
    #     #     loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m35_start"][i])+"-"+str(df_tss["loc_m10_end"][i])
    #     # else:
    #     #     loc = ">"+df_tss["genome"][i]+":"+str(df_tss["loc_m10_end"][i])+"-"+str(df_tss["loc_m35_start"][i])
    #     # print(loc + "\n" + df_tss["Minus35"][i]+df_tss["Spacer"][i]+df_tss["Minus10"][i], file=p_file)
    #     # arr.append([g, int(df_tss["loc_m35_start"][i])-1, df_tss["loc"].iloc[i]-1+len(df_tss[r][i]), df_tss["gene"].iloc[i], 0, df_tss["dir"].iloc[i]])
        
    #     arr.append([g, int(df_tss["loc_m35_start"][i])-1, int(df_tss["loc_m10_end"][i]), df_tss["gene"].iloc[i], 0, df_tss["dir"].iloc[i]])
        
    # df = pd.DataFrame(arr)
    # df.to_csv("../data/datasets/"+name+"/positive.bed", sep = "\t", header=0, index=0)

    # print("Generate", name, "positive.bed successfully.")

    file_g = "../raw_data/genome/" + g + ".gb"
    genome = Genome(next(SeqIO.parse(open(file_g, "r"), "genbank")))
    path = "../data/datasets/"+name+"/"
    file = open(path+"complete.fasta", 'w')
    print(">"+df_tss["genome"][0]+" "+name, file=file)
    print(genome.record.seq, file=file)
    print("Generate", name, ".fasta successfully.")