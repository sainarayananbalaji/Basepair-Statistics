from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import math

f='/Users/sainarayananbalaji/Documents/COLLEGE/CEMB/python interp/lambda_virus.fasta'

def GCATCounter(filepath):

    for record in SeqIO.parse(filepath, "fasta"):
        seq = str(record.seq).upper()
        length = len(seq)

        aCount= seq.count("A")
        cCount= seq.count("C")
        gCount= seq.count("G")
        tCount= seq.count("T")
        nCount= seq.count("N")

        gc_content = ((gCount + cCount)/length)*100
        at_content = ((aCount + tCount)/length)*100

        freqs = [seq.count(base) / length for base in "ATGC"]
        entropy = -sum(f * math.log2(f) for f in freqs if f>0)

        gc_skew = (gCount - cCount) / (gCount + cCount) if (gCount + cCount) != 0 else 0

        print(f"Sequence ID: {record.id}")
        print(f"Length: {length}")
        print(f"GC Content: {gc_content:.2f}%") 
        print(f"AT Content: {at_content:.2f}%")
        print(f"N Counts: {nCount}")
        print(f"Shannon Entropy: {entropy:.4f}")
        print(f"GC Skew: {gc_skew:.4f}")
        print()



GCATCounter(f)