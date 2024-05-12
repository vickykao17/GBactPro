import pandas as pd
import numpy  as np 
import os
import os.path
import traceback
import subprocess
import csv
from utils import readData, fastaToCharArray, fastaToHotEncoding, fastaToTetraNucletideDic
from utils import tetranucleotide_list_to_string_list, tetraToHotEncoding, promoterToTetraFreq
from utils import joinPositiveAndNegative, generateTrainAndTestSplit, h1
from utils import fasta_to_tetranucleotide_list, tetranucleotide_list_to_string_list
from keras.preprocessing.text import Tokenizer
from keras.preprocessing.sequence import pad_sequences
import pickle
import progressbar
import warnings
import datetime
import time
from Bio.Seq import Seq 
import joblib
# from pybedtools import BedTool
import pickle
warnings.filterwarnings("ignore")
from Bio import SeqIO

using_unbalanced = False

# regions = ["UP", "Minus35", "Spacer", "Minus10", "DOWN", "35+s+10", "35+10"]
regions = ["UP"]
def createIndex():
  df = pd.DataFrame([
    #TRAINING
      ["ECOLI_2"      , "NC_000913.2"  , True],
      ["ECOLI"        , "NC_000913.3"  , True],
      ["HPYLORI"      , "NC_000915.1"  , True],
      ["HPYLORI_2"    , "NC_000915.1"  , True],
      ["C_JEJUNI"     , "NC_002163.1"  , True],
      ["CJEJUNI_2"    , "NC_002163.1"  , True],
      ["CJEJUNI_3"    , "NC_003912.7"  , True],
      ["CJEJUNI_4"    , "NC_009839.1"  , True],
      ["CJEJUNI_5"    , "NC_008787.1"  , True],
      ["SPYOGENE"     , "LR031521.1"   , True],   
      ["STYPHIRMURIUM", "NC_016810.1"  , True], 
      ["CPNEUMONIAE"  , "NC_000922.1"  , True],  
      ["SONEIDENSIS"  , "NC_004347.2"  , True],
      ["LINTERROGANS" , "NZ_LT962963.1", True],  
      ["SCOELICOLOR"  , "NC_003888.3"  , True],
      #VALIDATION
      ["MYCOBACTER"   , "NC_008596"   , False],
      ["CLOSTRIDIUM"  , "NC_010001.1" , False],
      ["RHODOBACTER_1", "NC_014034.1" , False],
      ["RHODOBACTER_2", "NC_014034.1" , False],
      ["BACILLUS"     , "CP002927.1" , False],
      # ["AERUGINOSA"   , "NC_002516.2" , False]
  ], columns=["BACTERIA", "ID", "IS_TRAINING"])
  return df

def saveIndex(df, dir):
  df.to_csv(dir)

def generateData(
  bacteria_index, save_csv=False, save_data=True, 
  out_dir="./data/promoters_new/", nextflow_path='./nextflow',
  nextflow_pipeline = "pipeline_unbalanced.nf", # 'pipeline_without_docker.nf'
  manually_balance_data = False
):
  if(using_unbalanced):
    print("GENERATING UNBALANCED DATA WITH RATIO 1:10")
  else:
    print("GENERATE DATA")

  # bacteriaDir = "./bacteria"
  bacteria_report = {}
  if bacteria_index is None:
    index = createIndex()
  else: 
    index = bacteria_index

  data_root = "./data/"
  if not os.path.exists(data_root):
      os.makedirs(data_root)

  w = csv.writer(open(data_root+"report.csv", "w"))
  vocab_size = None
  tokenizer = None


  start_time = datetime.datetime.now().time().strftime('%H:%M:%S')
  bar = progressbar.ProgressBar(max_value=len(index))
  for i, row in index.iterrows() :
    for r in regions:
      bacteria_start_time = datetime.datetime.now().time().strftime('%H:%M:%S')

      # print("\n\n", 20*"*", i+1, f". {row['BACTERIA']}", 20*"*" )
      print("\n\n {} {} {} {}".format(20*"*", i+1, row['BACTERIA'], r, 20*"*"  ) )
      #nextflow run main_pipeline.nf --bacteria ecoli && rsync outDir/ outDirOriginal/ -a --copy-links -v
      print("\n\n {} {} {} {}".format(20*"*", i+1, "NEXTFLOW DATA GENERATION", 20*"*"  ) )
      # print("\n\n", 10*"*", "NEXTFLOW DATA GENERATION",10*"*" )
      
      stderr = None      
      stdout = None

      len_list = []
      fasta_sequences = SeqIO.parse(open('scan_data/'+row['BACTERIA']+'/'+r+'/positive.fasta'), 'fasta')

      for fasta in fasta_sequences:
          name, sequence = fasta.id, str(fasta.seq)
          # print(len(sequence))
          len_list.append(len(sequence))


      print("\n\nGENERATING NEXTFLOW DATA USING PIPELINE: ", nextflow_pipeline, "\n\n")
      out = subprocess.Popen([
          nextflow_path,
          'run',
          nextflow_pipeline,   #'pipeline_without_docker.nf',   #    pipeline_unbalanced_without_docker.nf   'main_pipeline.nf',  
          '--bacteria',
          str(row['BACTERIA']),
          '--region', 
          str(r)
          # len_list
      ],
      stdout=subprocess.PIPE,
      stderr=subprocess.STDOUT)
      stdout, stderr = out.communicate()
      error_msg = ""

      print("\n\nOUT: \n\n", stdout)
      print("\n\nERRORS: \n\n ",  stderr)


      
generateData(
    pd.DataFrame([

    # ["Synechococcus elongatus UTEX 2973"          , "CP006471.1"     , "Cyanobacteria"],
    # ["Synechocystis sp. PCC 6803"                    , "NC_000911.1"    , "Cyanobacteria"],
    # ["Bacteroides thetaiotaomicron VPI-5482"     , "NC_004663.1"    , "Bacteroides"],
    ["Campylobacter jejuni NCTC11168"             , "NC_002163.1"    , "Epsilonproteobacteria"],
    # ["Escherichia coli MG1655"                 , "NC_000913.2"    , "Gammaproteobacteria"],
    # ["Nanosynbacter_lyticus"            , "NZ_CP007496.1"  , "CPR"],
    # ["Synechococcus elongatus UTEX 2973"          , "CP006471.1"     , "Cyanobacteria"],
    # ["Synechocystis sp. PCC 6803"                    , "NC_000911.1"    , "Cyanobacteria"],
    # ["Anabaena"                         , "NC_003272.1"    , "Cyanobacteria"],
    # ["Thermotoga_maritima"              , "CP004077.1"     , "Thermotogae"],
    # ["Geobacter_sulfurreducens"         , "NC_002939.5"    , "Desulfuromonadia"],
    # ["Bacillus_subtilis"                , "NC_000964.3"    , "Firmicutes"],
    # ["Staphylococcus_aureus"            , "NC_003923.1"    , "Firmicutes"],
    # ["Staphylococcus_epidermidis"       , "NC_004461.1"    , "Firmicutes"],
    # ["Paenibacillus_riograndensis"      , "NZ_LN831776.1"  , "Firmicutes"],   
    # ["Listeria_monocytogenes"           , "NC_003210.1"    , "Firmicute"], 
    # ["Listeria_innocua"                 , "NC_003212.1"    , "Firmicute"],  
    # ["Streptococcus_pyrogenes"          , "LR031521.1"     , "Firmicute"],
    # ["Streptococcus_agalactiae"         , "NC_004368.1"    , "Firmicute"],  
    # ["Clostridioides_difficile"         , "NC_009089.1"    , "Clostridia; Firmicute"],
    # ["Lachnoclostridium_phytofermentans", "NC_010001.1"    , "Clostridia; Firmicute"],
    # ["Mycoplasma_pneumoniae"            , "NC_000912.1"    , "Mollicutes; Firmicute"],
    # ["Mycobacterium_tuberculosis"       , "AL123456.3"     , "Actinobacteria"],
    # ["Mycolicibacterium_smegmatis"      , "NC_008596.1"    , "Actinobacteria"],
    # ["Corynebacterium_glutamicum"       , "BX927147.1"     , "Actinobacteria"],
    # ["Corynebacterium_diphtheriae"      , "NC_002935.2"    , "Actinobacteria"],
    # ["Streptomyces_coelicolor"          , "NC_003888.3"    , "Actinobacteria"],
    # ["Chlamydia_pneumoniae"             , "NC_000922.1"    , "Chlamydiae"],
    # ["Chlamydia_trachomatis"            , "NC_010280.1"    , "Chlamydiae"],
    # ["Bacteroides thetaiotaomicron VPI-5482"     , "NC_004663.1"    , "Bacteroides"],
    # ["Leptospira_interrogans"           , "NZ_LT962963.1"  , "Spirochaetes"],
    # ["Helicobacter_pylori"              , "NC_000915.1"    , "Epsilonproteobacteria"],
    # ["Campylobacter jejuni NCTC11168"             , "NC_002163.1"    , "Epsilonproteobacteria"],
    # ["Agrobacterium_tumefaciens"        , "NC_003062.2"    , "Alphaproteobacteria"],
    # ["Methylorubrum_extorquens"         , "NC_012988.1"    , "Alphaproteobacteria"],
    # ["Caulobacter_crescentus"           , "NC_011916.1"    , "Alphaproteobacteria"],
    # ["Zymomonas_mobilis"                , "NZ_CP023715.1"  , "Alphaproteobacteria"],
    # ["Rhodobacter_sphaeroides"          , "NC_007493.2"    , "Alphaproteobacteria"],
    # ["Bradyrhizobium_japonicum"         , "NC_004463.1"    , "Alphaproteobacteria"],
    # ["Burkholderia_cenocepacia"         , "AM747720.1"     , "Betaproteobacteria"],
    # ["Neisseria_gonorrhoeae"            , "NC_022240.1"    , "Betaproteobacteria"],
    # ["Bordetella_pertussis"             , "NC_002929.2"    , "Betaproteobacteria"],
    # ["Klebsiella_aerogenes"             , "NC_015663.1"    , "Gammaproteobacteria"],
    # ["Pseudomonas_aeruginosa"           , "NC_008463.1"    , "Gammaproteobacteria"],
    # ["Acinetobacter_baumannii"          , "CP000521.1"     , "Gammaproteobacteria"],
    # ["Escherichia coli MG1655"                 , "NC_000913.2"    , "Gammaproteobacteria"],
    # ["Salmonella_typhimurium"           , "FQ312003.1"     , "Gammaproteobacteria"],
    # ["Xanthomonas_campestris"           , "AM920689.1"     , "Gammaproteobacteria"],
    # ["Shewanella_oneidensis"            , "NC_004347.2"    , "Gammaproteobacteria"],
    # ["Vibrio_cholerae"                  , "NC_002505.1"    , "Gammaproteobacteria"],
    # ["Shigella_flexneri"                , "CP037923.1"     , "Gammaproteobacteria"],      
    # ["Pseudomonas_putida"               , "NC_002947.3"    , "Gammaproteobacteria"],  
    # ["Phytoplasma"                      , "AP006628.2"     , "Mollicutes; Firmicute"],    
    ], columns=["BACTERIA", "ID", "Lineage"])
  )
