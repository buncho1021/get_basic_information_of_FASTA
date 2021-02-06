from sys import argv
import numpy as np
from Bio import SeqIO

def main():
    Sseqlen, Cseqlen, N, GC, gap=[], [], 0, 0, 0
    for line in SeqIO.parse(argv[1], "fasta"):
        sequence = line.seq.lower() #小文字にする
        len_sequence = len(sequence)
        if "n" in sequence:#scaffold
            N   += sequence.count("n")
            Sseqlen += [len_sequence]
            for j in range(1, len_sequence):#gap counting
                if sequence[j-1] != "n" and sequence[j] == "n":
                    gap += 1   
        else:#contig
            Cseqlen += [len_sequence]
        GC  += sequence.count("g") + sequence.count("c")
    length = Sseqlen + Cseqlen ;length.sort(reverse=True)
    sum_length = np.sum(length)

    lcum, half = np.cumsum(length), sum_length*.5
    N50 = length[np.where(lcum == lcum[lcum >= half][0])[0][0]]
    L50 = len(lcum[lcum <= half]) + 1

    print("Total length\t", "{:,}".format(sum_length))
    print("Number of contigs\t", len(length) + gap)
    print("Number of scaffolds\t", len(Sseqlen))
    print("Number of gaps\t", gap)
    print("Number of Ns\t", N, "\n")
    if len(Cseqlen) !=0:
        print("Max contig length\t", np.max(Cseqlen))#scaffoldを構成するcontigは考慮していません。
        print("Minmum contig length\t", np.min(Cseqlen))
        print("Mean contig length\t", round(np.mean(Cseqlen), 2))
        print("Median contig length\t", np.median(Cseqlen),"\n")
    print(f'N50\t{N50}\nL50\t{L50}')
    print("GC content\t",round(GC/sum_length, 3),'\n')
    if len(Sseqlen) != 0:
        print("Max scaffold length\t", np.max(Sseqlen))
        print("Minmum scaffold length\t", np.min(Sseqlen))
        print("Mean scaffold length\t", round(np.mean(Sseqlen), 2))
        print("Median scaffold length\t", np.median(Sseqlen))

if __name__ == "__main__":
    main()