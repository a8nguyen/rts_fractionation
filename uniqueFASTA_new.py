#!/usr/bin/env python
import sys
from Bio import SeqIO
from collections import defaultdict
from pyopenms import *

input_fasta = sys.argv[1]

dig = ProteaseDigestion()
dig.setMissedCleavages(1)

pep2prots = defaultdict(set)
prot2peps = defaultdict(set)

with open(input_fasta, "r") as handle:
    for record in SeqIO.parse(input_fasta, "fasta"):
        if "|" in record.id:
            acc = record.id.split("|")[1].strip()
        else:
            acc = "CONTAM_" + record.id.strip()
        peptides = []
        seq = AASequence.fromString(str(record.seq))
        dig.digest(seq, peptides)
        for pep in peptides:
            s = pep.toString().decode("utf-8")
            if "X" not in s and len(s) >= 5 and pep.getMonoWeight() <= 8000:
                pep2prots[s].add(acc)
                prot2peps[acc].add(s)
uni_genes = []
for prot in prot2peps:
    canonical = prot.split("-")[0]
    isoform = prot if canonical != prot else prot +'-1'

    for pep in prot2peps[prot]:
        # all peptides for protein
        print(">" + isoform + "_all")
        print(pep)
        # unique peptides for gene
        is_unique_to_gene = True
        for matching_prot in pep2prots[pep]:
            if canonical not in matching_prot:
                is_unique_to_gene = False
                uni_genes.append(matching_prot)
                break
        if is_unique_to_gene:
            print(">" + isoform + "_uni_gene")
            print(pep)
        # unique peptides for protein isoform
        if len(pep2prots[pep]) == 1:
            print(">" + isoform + "_uniso")
            print(pep)
                


        
