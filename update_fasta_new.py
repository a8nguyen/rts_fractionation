import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import defaultdict
import argparse
import os
import glob
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('file', nargs=2)
args = parser.parse_args()


def parameter_file_parser(param_file):
    params = {}
    with open(param_file) as p:

        for line in p:
            if line.startswith("#") is False and line != '\n':
                lista = list(map(lambda x: x.strip(), line.split("=")))
                params[lista[0]] = lista[1]
    return params


parameter = parameter_file_parser(args.file[0])
result = max(glob.iglob(os.path.join(parameter['input_csv'],'*_realtimesearch.csv')), key=os.path.getmtime)
fasta = max(glob.iglob(os.path.join(parameter['out_dir'],'fasta','*.fasta')), key = os.path.getmtime) 
out_path = parameter['out_dir']
if not os.path.exists(out_path):
    os.makedirs(out_path)

print(f"Step 1: Reading rts results and filtering out closed-out proteins...")
print("Latest rts run is %s" %result)
closedoutPeps = int(parameter['closedoutPeps'])
Xcorr_thres = float(parameter['Xcorr_thres'])
dCN = float(parameter['dCN'])
mass_tolerance = int(parameter['mass_tolerance'])
max_pep_per_prot = int(parameter['max_pep_per_prot'])

rts_df = pd.read_csv(result, delimiter= '\t')
rts_df = rts_df[rts_df['Protein ID'].str.contains('CONTAM|CLOSED_OUT|DECOY') == False]
rts_df = rts_df[(rts_df['Passed LDA'] == True) | ((rts_df['Xcorr'] >= Xcorr_thres) & (
    rts_df['dCn'] >= 0.1) & (rts_df['Precursor PPM'] <= mass_tolerance))]

#keeping track ms3 scan
peptideAccessionCount = defaultdict(int)
proteinCount = defaultdict(int)
fractionRepeat = 0
ms3notTriggered = 0

print(f"Step 2: Updating the count of peptides...")
for _, row in rts_df.iterrows():
    ids = set(r.split('|')[0] for r in row["Protein ID"].split(' : '))
    proteins = set(r.split('_',1 )[0] for r in ids)
    for p in proteins:                  # protein count
        proteinCount[p] += 1
    seq = row["Peptide"]
    with open(out_path+'/peptides_rts.txt', 'a+') as pep_rts:
        pep_rts.seek(0, os.SEEK_SET)
        any_repeat = any(seq in line for line in pep_rts)
        if not any_repeat:
            pep_rts.write(seq + '\n')
            ms3 = False
            for accession in ids:
                prot= accession.split('_', 1)[0]
                if proteinCount[prot] < max_pep_per_prot:       # if it is ms3 triggered
                    peptideAccessionCount[accession] += 1       # add all accessions associated with such protein
                    ms3 = True
                else:
                    break
            ms3notTriggered += (not ms3)                        # keeping tab for how many ms3 is not triggered
        else:
            fractionRepeat +=1                                  # keeping tab for how many fraction gets repeated
            #print(seq)
print(f"Step 2 finished.")
print(f"     %d repeated scans." %fractionRepeat )
print(f"     %d/%d scans exceed max peptide per protein." %(ms3notTriggered, rts_df.shape[0]) )


print(f"Step 3: Writing a new fasta...")
fasta_out = open(args.file[1], 'w') 
with open(fasta, 'r') as handle:
    for record in SeqIO.parse(fasta, "fasta"):
        acc = record.id.split("|")[0]
        try:
            accCount = int( record.id.split("|")[1] ) 
        except IndexError:
            accCount = 0
        accCount += peptideAccessionCount[acc]
        if accCount >= closedoutPeps and "CLOSED_OUT" not in acc:
            acc += "_CLOSED_OUT" 
        acc += '|' + str(accCount)
        seq = record.seq
        print(">" + acc, file = fasta_out)
        print(seq, file = fasta_out)
print(f"Finished updating new fasta.")

print("Moving CSV results into output directory...")
result_files = [os.path.splitext(result)[0] +ext for ext in ['.csv', '.raw'] ] 
for f in result_files:
    shutil.move(f, os.path.join(out_path, 'results'))
# filter for newly closed out proteins
# write new FASTA with newly closed protein annotations
# compress new FASTA
# overwrite FASTA in method file

