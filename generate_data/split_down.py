from Bio import SeqIO
import os

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
            'Streptomyces coelicolor A3': 1016,
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


for b in bacteria.keys():
    print(b)

    p_file = "../data/datasets/"+b+"/DOWN/random/positive.fasta"
    n_file = "../data/datasets/"+b+"/DOWN/random/negative_new.fasta"

    _1_8_dir = "../data/datasets/"+b+"/DOWN_1_8/random/"
    if not os.path.exists(os.path.dirname(_1_8_dir)):
        os.makedirs(os.path.dirname(_1_8_dir))

    _9_60_dir = "../data/datasets/"+b+"/DOWN_9_60/random/"
    if not os.path.exists(os.path.dirname(_9_60_dir)):
        os.makedirs(os.path.dirname(_9_60_dir))

    _1_8_p_path = _1_8_dir + "positive.fasta"
    _1_8_n_path = _1_8_dir + "negative.fasta"
    _1_8_all_path = _1_8_dir + "all.fasta"

    _1_8_p_file = open(_1_8_p_path, 'w')
    _1_8_n_file = open(_1_8_n_path, 'w')
    _1_8_all_file = open(_1_8_all_path, 'w')

    _9_60_p_path = _9_60_dir + "positive.fasta"
    _9_60_n_path = _9_60_dir + "negative.fasta"
    _9_60_all_path = _9_60_dir + "all.fasta"

    _9_60_p_file = open(_9_60_p_path, 'w')
    _9_60_n_file = open(_9_60_n_path, 'w')
    _9_60_all_file = open(_9_60_all_path, 'w')

    
    fasta_sequences = SeqIO.parse(open(p_file),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)

        print('>'+name+'\n'+sequence[:8], file=_1_8_p_file)
        print('>'+name+'\n'+sequence[:8], file=_1_8_all_file)
        print('>'+name+'\n'+sequence[9:], file=_9_60_p_file)
        print('>'+name+'\n'+sequence[9:], file=_9_60_all_file)

    fasta_sequences = SeqIO.parse(open(n_file),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)

        print('>'+name+'\n'+sequence[:8], file=_1_8_n_file)
        print('>'+name+'\n'+sequence[:8], file=_1_8_all_file)
        print('>'+name+'\n'+sequence[9:], file=_9_60_n_file)
        print('>'+name+'\n'+sequence[9:], file=_9_60_all_file)