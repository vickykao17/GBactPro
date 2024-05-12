from Bio import SeqIO
from common import *


genome = Genome(next(SeqIO.parse(open(file_g, "r"), "genbank")))