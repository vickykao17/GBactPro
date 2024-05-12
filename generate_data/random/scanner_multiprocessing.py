from time import time
from multiprocessing import Pool
from scanner import *
import os


# Function to process each TSS group in parallel
def process_tss_group(args):
    _, df_tss = args
    tss_m10_range = [0, 20]
    
    g = df_tss["genome"].iloc[0]
    name = df_tss["name"].iloc[0]
    file_g = "../raw_data/genome/" + g + ".gb"
    genome = Genome(next(SeqIO.parse(open(file_g, "r"), "genbank")))
    genome_copied = deepcopy(genome)  # for tss_shuffle
    tax = genome.record.annotations["taxonomy"]

    print(name)

    tss = []
    tss_random = []
    tss_shuffle = []



    for idx, row in df_tss.iterrows():
        direction = row["dir"]
        tss_loc = int(row["loc"])

        direction_random = random.choice(['+', '-'])
        tss_loc_random = random.choice(range(100, genome.length-100)) + 1

        tss.append(scan_promoter(tss_loc, direction, genome, tss_m10_range))
        tss_random.append(scan_promoter(tss_loc_random, direction_random, genome, tss_m10_range))
        tss_shuffle.append(scan_promoter_shuffle(tss_loc, direction, genome_copied, tss_m10_range))

    df_tss["loc_m10_end"] = [t[0].getLoc_m10_end() for t in tss]
    df_tss["loc_m35_start"] = [t[0].getLoc_m35_start() for t in tss]
    df_tss["UP"] = [str(get_element_UP(t[0])) for t in tss]
    df_tss["Minus35"] = [str(t[0].getSequence()[:6]) for t in tss]
    df_tss["Spacer"] = [str(t[0].getSequence()[6:-6]) for t in tss]
    df_tss["Minus10"] = [str(t[0].getSequence()[-6:]) for t in tss]
    df_tss["DOWN"] = [str(get_element_DOWN(t[0])) for t in tss]
    df_tss["label"] = ['promoter' for t in tss]
    # df_tss.to_csv(f'../table/{name}.tsv', sep='\t', index=False)
    df_tss.to_csv(f'../data/scan_result/{name}.tsv', sep='\t', index=False)

    df_tss_random = pd.DataFrame({
        "loc": [t[1][0] for t in tss_random],
        "dir": [t[1][1] for t in tss_random],
        "loc_m10_end": [t[0].getLoc_m10_end() for t in tss], 
        "UP": [str(get_element_UP(t[0])) for t in tss_random],
        "Minus35": [str(t[0].getSequence()[:6]) for t in tss_random],
        "Spacer": [str(t[0].getSequence()[6:-6]) for t in tss_random],
        "Minus10": [str(t[0].getSequence()[-6:]) for t in tss_random],
        "DOWN": [str(get_element_DOWN(t[0])) for t in tss_random],
    })
    # df_tss_random.to_csv(f'../table/{name}_random.tsv', sep='\t', index=False)
    df_tss_random.to_csv(f'../data/scan_result/{name}_random.tsv', sep='\t', index=False)

    df_tss_shuffle = deepcopy(df_tss)
    df_tss_shuffle["loc_m10_end"] = [t[0].getLoc_m10_end() for t in tss_shuffle]
    df_tss_shuffle["UP"] = [str(get_element_UP(t[0])) for t in tss_shuffle]
    df_tss_shuffle["Minus35"] = [str(t[0].getSequence()[:6]) for t in tss_shuffle]
    df_tss_shuffle["Spacer"] = [str(t[0].getSequence()[6:-6]) for t in tss_shuffle]
    df_tss_shuffle["Minus10"] = [str(t[0].getSequence()[-6:]) for t in tss_shuffle]
    df_tss_shuffle["DOWN"] = [str(get_element_DOWN(t[0])) for t in tss_shuffle]
    # df_tss_shuffle.to_csv(f'../table/{name}_shuffle.tsv', sep='\t', index=False)
    df_tss_shuffle.to_csv(f'../data/scan_result/{name}_shuffle.tsv', sep='\t', index=False)


    # out_dir = "../data/positive/" + name + "/"
    # if not os.path.exists(os.path.dirname(out_dir)):
    #     os.makedirs(os.path.dirname(out_dir))

    # file_path = open(out_dir + "positive_scan.fasta", "w")
    
    # for i in range(len(tss)):
    #     if df_tss["dir"][i] == '+':
    #         print(">"+df_tss["genome"][i]+":"+str(df_tss["loc_m35_start"][i])+"-"+str(df_tss["loc_m10_end"][i]), file=file_path)
    #     else:
    #         print(">"+df_tss["genome"][i]+":"+str(df_tss["loc_m10_end"][i])+"-"+str(df_tss["loc_m35_start"][i]), file=file_path)
    #     print(df_tss["Minus35"][i]+df_tss["Spacer"][i]+df_tss["Minus10"][i], file=file_path)
    #     print(">"+tss["genome"][i]+str(tss["promoter"][i][0].l)+":"+str(tss["promoter"][i][0].r))
    #     print(tss["promoter"][i][0].getSequence())
    # tss.to_csv(out_dir + '/report.csv')


if __name__ == '__main__':

    t_start = time()

    # tss_m10_range = [-10, 20]
    # print(tss_m10_range)

    # TSS = pd.read_csv("../data/tss_20240122.tsv", sep="\t")
    TSS = pd.read_csv("../raw_data/tss.syn_e.tsv", sep="\t")

    # path = '../data/scan_result/' + name
    # if not os.path.exists(os.path.dirname(path)):
    #     os.makedirs(os.path.dirname(path))
    # Split TSS dataframe into chunks for parallel processing
    chunks = [(i, group) for i, group in TSS.groupby("name")]

    # Create a Pool of processes and map the processing function to the chunks
    with Pool(4) as pool:
        pool.map(process_tss_group, chunks)

    t_end = time()
    print(f"time usage: {int(t_end-t_start)} seconds")