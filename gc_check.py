from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def GCCounter(filepath):

    for record in SeqIO.parse(filepath, "fasta"):
        sequence = str(record)
        gc_content = gc_fraction(sequence)
        print(f"Sequence ID: {record.id}, GC content: {gc_content:.2f}%")

GCCounter('/Users/sainarayananbalaji/Documents/COLLEGE/CEMB/python interp/lambda_virus.fasta')