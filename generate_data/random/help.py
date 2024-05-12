import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--input-file')
parser.add_argument('--output-file')
args = parser.parse_args()

f1 = open(args.input_file,"r")
f2 = open(args.output_file,"w")

seq = ''
for line in f1:
    if line.startswith('>'):
        if seq != '':
            seq += '\n'
            f2.write(seq)
            seq = ''
        f2.write(line)
    else:
        seq += line[:-1]

seq += '\n'
f2.write(seq)
    

    
