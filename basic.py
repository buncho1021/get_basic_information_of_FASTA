import sys
import numpy as np
from Bio import SeqIO

class get_basic:
    def __init__(self, file_path=None):
        self.file_path = file_path
        self.Sseqlen, self.Cseqlen = [], [],
        self.N, self.GC, self.gap = 0, 0, 0
    
    def main(self):
        
        for line in SeqIO.parse(self.file_path, "fasta"):
            sequence = line.seq
            len_sequence = len(sequence)
            if "n" in sequence.lower():#scaffold
                self.Sseqlen += [len_sequence]
                for j in range(1, len_sequence):
                    if sequence[j-1]!="N" and sequence[j]=="N":
                        self.gap += 1  
            else:#contig
                self.Cseqlen += [len_sequence]
            self.N   += sequence.lower().count("n")
            self.GC  += sequence.lower().count("g") + sequence.lower().count("c")
        self.length = self.Sseqlen + self.Cseqlen ;self.length.sort(reverse=True)
        self.sum_length = np.sum(self.length)

        lcum, half = np.cumsum(self.length), self.sum_length*.5
        self.N50 = self.length[np.where(lcum == lcum[lcum >= half][0])[0][0]]
        self.L50 = len(lcum[lcum <= half]) + 1

    def printer(self):
        print("Total length\t","{:,}".format(self.sum_length))
        print("Number of contigs\t",len(self.length) + self.gap)
        print("Number of scaffolds\t",len(self.Sseqlen))
        print("Number of gaps\t",self.gap)
        print("Number of Ns\t",self.N,"\n")
        if len(self.Cseqlen)!=0:
            print("Max contig length\t",np.max(self.Cseqlen))#scaffoldを構成するcontigは考慮していません。update予定
            print("Minmum contig length\t",np.min(self.Cseqlen))
            print("Mean contig length\t",round(np.sum(self.Cseqlen)/len(self.Cseqlen), 2))
            print("Median contig length\t",np.median(self.Cseqlen),"\n")
        print("N50\t", self.N50)
        print("L50\t", self.L50)
        print("GC content\t",round(self.GC/self.sum_length, 3),'\n')
        if len(self.Sseqlen) != 0 :
            print("Max scaffold length\t",np.max(self.Sseqlen))
            print("Minmum scaffold length\t",np.min(self.Sseqlen))
            print("Mean scaffold length\t",np.sum(self.Sseqlen)/len(self.Sseqlen))
            print("Median scaffold length\t",np.median(self.Sseqlen))

if __name__ == "__main__":
    instance = get_basic(file_path=sys.argv[1])
    instance.main() ;instance.printer()
    
    
    
    
    
    
    
