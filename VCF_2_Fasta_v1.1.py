##########################################################
### Import Necessary Modules

import argparse                       #provides options at the command line
import sys                       #take command line arguments and uses it in the script
import gzip                       #allows gzipped files to be read
import re                       #allows regular expressions to be used
import subprocess                   #allows external programs to be called

##########################################################
### Command-line Arguments
parser = argparse.ArgumentParser(description="A script to read in a vcf file and return a fasta file.  Requires samtools to be installed and in path.")
parser.add_argument("-vcf", help = "The location of the VCF file (accepts .gz format)", default=sys.stdin, required=True)
parser.add_argument("-fasta", help = "The location of the fasta file. Can be indexed (which will speed up the script).", default=sys.stdin, required=True)
parser.add_argument("-fai", help = "The location of the index of a fasta file (.fai format), optional, default = NA", default="NA")
parser.add_argument("-flanking", help = "The distance of flanking sequence around the SNP variant, if the sequence is too close to an edge it will not be returned, default = 50", default=50)
args = parser.parse_args()

#########################################################
### Open file (object-oriented programming)

class OpenFile():
    sequences = {}
    def __init__ (self, v, t):
        if re.search(".gz$", v):
            self.file = gzip.open(v, 'r')        
        else:
            self.file = open(v, 'r')
        if t == "vcf":
            sys.stderr.write("Opened vcf file\n")             
            self.readLinesVCF(self.file)
        elif t == "fai":
            sys.stderr.write("Opened index of fasta file\n")             
            self.readLinesFastaFai(self.file)
        else:
            sys.stderr.write("Opened fasta file\n")             
            self.readLinesFasta(self.file)

    def readLinesVCF(self, vc):
        self.open_vcf = vc
        self.variants_found = 0
        self.variants_used = 0
        for line in self.open_vcf:
            try:
                line = line.decode('utf-8')
            except:
                pass
            line = line.rstrip('\n')
            if not re.search("^#", line):     
                line_list = line.split("\t")
                self.chrom = line_list[0]
                self.pos = line_list[1]
                self.ref = line_list[3]
                self.alt = line_list[4]
                self.start = int(self.pos) - int(args.flanking)
                self.end = int(self.pos) + int(args.flanking)
                self.snpPosInFlanking = int(args.flanking) + 1
                self.variants_found += 1
                if int(self.start) > 0 and int(OpenFile.sequences[self.chrom]) >= int(self.end):
                    ###Affx-88959408:A:G:36
                    self.location = str(self.chrom) + ":" + str(self.start) + "-" + str(self.end)
                    self.header = str(self.chrom) + "_" + str(self.pos) + ":" + str(self.ref) + ":" + str(self.alt) + ":" + str(self.snpPosInFlanking)
                    self.sequence = "NA"
                    process = subprocess.Popen(["samtools", "faidx", args.fasta, self.location], stdout=subprocess.PIPE)
                    process.wait()
                    for line in process.stdout:
                        try:
                            line = line.decode('utf-8')
                        except:
                            pass
                        line = line.rstrip("\n")
                        if not re.search("^>", line) and re.search("\w", line):
                            if self.sequence == "NA":
                                self.sequence = line
                            else:
                                self.sequence += line
                    if self.sequence != "NA":
                        print (">{}\n{}".format(self.header, self.sequence)) 
                        self.variants_used += 1
        self.open_vcf.close()
        sys.stderr.write("\tFinished reading vcf file: Found {} variants, used {}\n".format(self.variants_found, self.variants_used))

    def readLinesFasta(self, f):
        """Measures the lengths of the scaffolds in the fasta file"""
        self.filename = f
        self.header = "NA"
        self.seq = ""
        self.number_scaffolds = 0
        self.total_size = 0
        for line in self.filename:
            try:
                line = line.decode('utf-8')
            except:
                pass
            line = line.rstrip('\n')
            if re.search("^\>", line):
                if self.header == "NA":
                    self.header = line.split(" ")[0][1:]
                else:
                    if not (int(len(self.header)) > 0 and int(len(self.seq)) > 0):
                        sys.stderr.write("\tSequence Found not conforming to fasta: {}\n".format(self.header))
                    OpenFile.sequences[self.header] = int(len(self.seq))
                    self.total_size += int(len(self.seq))
                    self.seq = ""
                    self.header = line.split(" ")[0][1:]
                self.number_scaffolds += 1
            elif re.search("\w", line):
                self.seq += line
        if not (int(len(self.header)) > 0 and int(len(self.seq)) > 0):
            sys.stderr.write("\tSequence Found not conforming to fasta: {}\n".format(self.header))
        OpenFile.sequences[self.header] = int(len(self.seq))
        self.total_size += int(len(self.seq))
        self.seq = ""
        sys.stderr.write("\tFinished reading scaffold fasta file: Found {} sequence(s)\n".format(self.number_scaffolds)) 
        sys.stderr.write("                    Total Nucleotide(s): {}\n\n".format(self.total_size))
        self.filename.close()   

    def readLinesFastaFai(self, f):
        """Measures the lengths of the scaffolds in the fai file"""
        self.filename = f
        self.number_scaffolds = 0
        self.total_size = 0
        for line in self.filename:
            try:
                line = line.decode('utf-8')
            except:
                pass
            line = line.rstrip('\n')
            (self.header, self.nucleotideCount, self.offset, self.lineBases, self.lineWidth) = line.split("\t")
            OpenFile.sequences[self.header] = int(self.nucleotideCount)
            self.total_size += int(self.nucleotideCount)
            self.number_scaffolds += 1
        sys.stderr.write("\tFinished reading scaffold fasta fai file: Found {} sequence(s)\n".format(self.number_scaffolds)) 
        sys.stderr.write("                    Total Nucleotide(s): {}\n\n".format(self.total_size))
        self.filename.close()   

if __name__ == '__main__':
    if args.fai == "NA":
        open_fasta = OpenFile(args.fasta, "fasta")
    else:
        open_fai = OpenFile(args.fai, "fai")
    open_vcf = OpenFile(args.vcf, "vcf")
