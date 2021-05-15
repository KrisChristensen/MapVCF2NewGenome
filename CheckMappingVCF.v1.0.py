##########################################################
### Import Necessary Modules

import argparse                       #provides options at the command line
import sys                       #take command line arguments and uses it in the script
import gzip                       #allows gzipped files to be read
import re                       #allows regular expressions to be used
import textwrap                       #allows the use of textwrapping for long sequences
import subprocess                   #allows external programs to be called

##########################################################
### Command-line Arguments
parser = argparse.ArgumentParser(description="A script to check that variants that were mapped with VCF_2_Fasta_v1.1.py, bwa mem, Bam2SNPPosition.v1.0.py pipeline were correct.  Input genomes must be indexed with samtools faidx tool.")
parser.add_argument("-map", help = "The location of the map file (format: contig<tab>position<tab>newContig<tab>newPosition<tab>length<tab>pid)(accepts .gz format)", default=sys.stdin, required=True)
parser.add_argument("-fasta1", help = "The location of the first fasta file (must match column one of map file). Must be indexed.", default=sys.stdin, required=True)
parser.add_argument("-fasta2", help = "The location of the second fasta file (must match column three of map file). Must be indexed.", default=sys.stdin, required=True)
parser.add_argument("-fai1", help = "The location of the index of a fasta file1 (.fai format), optional, default = NA", default="NA")
parser.add_argument("-fai2", help = "The location of the index of a fasta file2 (.fai format), optional, default = NA", default="NA")
parser.add_argument("-flanking", help = "The distance of flanking sequence around the SNP variant, if the sequence is too close to an edge it will not be checked, default = 10.  The sequences must be exact matches to be counted as accurate", default=10)
parser.add_argument("-output", help = "Allows certain mappings to be output, the default is to report stats only.  'exact' will only output mappings with the exact sequences, 'most' will output exact and mappings where the variant isn't the same, but everything else matches, and 'partial' will output exact, most, and matches from the variant and at least half of the sequence, default: NA, options: exact, partial, most", default="NA")
parser.add_argument("-verbose", help = "Allows printing of the sequences around the variant and why the category, default: NA, option: yes", default="NA")
args = parser.parse_args()

#########################################################
### Open file (object-oriented programming)
class Variables():
    sequences = {}

class OpenFile():
    def __init__ (self, v, t):
        if re.search(".gz$", v):
            self.file = gzip.open(v, 'r')        
        else:
            self.file = open(v, 'r')
        if t == "map":
            sys.stderr.write("Opened map file\n")             
            self.readLinesMap(self.file)
        elif t == "fai":
            sys.stderr.write("Opened index of fasta file\n")             
            self.readLinesFastaFai(self.file)

    def readLinesMap(self, vc):
        self.open_map = vc
        self.variants_found = 0
        self.variants_used = 0
        self.variants_unknown = 0
        self.variants_mismatch = 0
        self.variants_partial = 0
        self.variants_varMissing = 0
        for line in self.open_map:
            try:
                line = line.decode('utf-8')
            except:
                pass
            line = line.rstrip('\n')
            if not re.search("^#", line):
                self.out = line     
                line_list = line.split("\t")
                self.contig1 = line_list[0]
                self.pos1 = line_list[1]
                self.contig2 = line_list[2]
                self.pos2 = line_list[3]
                self.start1 = int(self.pos1) - int(args.flanking)
                self.end1 = int(self.pos1) + int(args.flanking)
                self.start2 = int(self.pos2) - int(args.flanking)
                self.end2 = int(self.pos2) + int(args.flanking)
                self.snpPosInFlanking = int(args.flanking) + 1
                self.variants_found += 1
                if ((int(self.start1) > 0 and int(Variables.sequences[self.contig1]) >= int(self.end1)) and
                (int(self.start2) > 0 and int(Variables.sequences[self.contig2]) >= int(self.end2))):
                    self.location1 = str(self.contig1) + ":" + str(self.start1) + "-" + str(self.end1)
                    self.location2 = str(self.contig2) + ":" + str(self.start2) + "-" + str(self.end2)
                    self.sequence1 = "NA"
                    self.sequence2 = "NA"
                    self.sequence2RC = "NA"                   
                    process = subprocess.Popen(["samtools", "faidx", args.fasta1, self.location1], stdout=subprocess.PIPE)
                    process.wait()
                    for line in process.stdout:
                        try:
                            line = line.decode('utf-8')
                        except:
                            pass
                        line = line.rstrip("\n")
                        if not re.search("^>", line) and re.search("\w", line):
                            if self.sequence1 == "NA":
                                self.sequence1 = line
                            else:
                                self.sequence1 += line
                    process = subprocess.Popen(["samtools", "faidx", args.fasta2, self.location2], stdout=subprocess.PIPE)
                    process.wait()
                    for line in process.stdout:
                        try:
                            line = line.decode('utf-8')
                        except:
                            pass
                        line = line.rstrip("\n")
                        if not re.search("^>", line) and re.search("\w", line):
                            if self.sequence2 == "NA":
                                self.sequence2 = line
                            else:
                                self.sequence2 += line
                    if self.sequence1 != "NA" and self.sequence2 != "NA":
                        self.sequence1 = self.sequence1.lower()
                        self.sequence2 = self.sequence2.lower()
                        self.sequence2RC = self.reverseComplement(self.sequence2)
                        if self.sequence1 == self.sequence2 or self.sequence1 == self.sequence2RC:
                            if args.verbose != "NA":
                                sys.stderr.write("exact-match\t{}\t{}\trc:{}\n".format(self.sequence1, self.sequence2, self.sequence2RC))
                            self.variants_used += 1
                            if args.output != "NA":
                                print("{}\t{}".format(self.out, "exact"))
                        else:
                            self.seq1a = self.sequence1[0:int(args.flanking) - 1]
                            self.seq1b = self.sequence1[int(args.flanking)]
                            self.seq1c = self.sequence1[int(args.flanking)+1:]
                            self.seq2a = self.sequence2[0:int(args.flanking) - 1]
                            self.seq2b = self.sequence2[int(args.flanking)]
                            self.seq2c = self.sequence2[int(args.flanking)+1:]
                            self.seq2ar = self.sequence2RC[0:int(args.flanking) - 1]
                            self.seq2br = self.sequence2RC[int(args.flanking)]
                            self.seq2cr = self.sequence2RC[int(args.flanking)+1:]
                            self.first = "no"
                            self.second = "no" 
                            self.third = "no"                                              
                            if self.seq1a == self.seq2a or self.seq1a == self.seq2ar:
                                self.first = "yes"
                            if self.seq1b == self.seq2b or self.seq1b == self.seq2br:
                                self.second = "yes"
                            if self.seq1c == self.seq2c or self.seq1c == self.seq2cr:
                                self.third = "yes"
                            if self.first == "yes" and self.third == "yes":
                                self.variants_varMissing += 1
                                if args.verbose != "NA":
                                    sys.stderr.write("most-match\t{}\t{}\trc:{}\n".format(self.sequence1, self.sequence2, self.sequence2RC))
                                if args.output != "NA" and args.output != "exact":
                                    print("{}\t{}".format(self.out, "most"))                                                                         
                            elif (self.first == "yes" or self.third == "yes") and self.second == "yes":
                                self.variants_partial += 1
                                if args.verbose != "NA":
                                    sys.stderr.write("partial-match\t{}\t{}\trc:{}\n".format(self.sequence1, self.sequence2, self.sequence2RC))
                                if args.output != "NA" and args.output != "exact" and args.output != "most":
                                    print("{}\t{}".format(self.out, "partial"))                                                                
                            else:
                                self.variants_mismatch += 1
                                if args.verbose != "NA":
                                    sys.stderr.write("no-match\t{}\t{}\trc:{}\n".format(self.sequence1, self.sequence2, self.sequence2RC))

                else:
                    self.variants_unknown += 1
        self.open_map.close()
        sys.stderr.write("\tFinished reading vcf file: Found {} variants, correct variants:{}, unknown variants:{}, CompleteMismatch:{}, variantOnlyMismatch:{}, partialMatch:{}\n".format(self.variants_found, self.variants_used, self.variants_unknown, self.variants_mismatch, self.variants_varMissing, self.variants_partial))
        
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
            Variables.sequences[self.header] = int(self.nucleotideCount)
            self.total_size += int(self.nucleotideCount)
            self.number_scaffolds += 1
        sys.stderr.write("\tFinished reading scaffold fasta fai file: Found {} sequence(s)\n".format(self.number_scaffolds)) 
        sys.stderr.write("                    Total Nucleotide(s): {}\n\n".format(self.total_size))
        self.filename.close()
        
    def reverseComplement (self, seq):
        """Finds the reverse complement of sequence"""
        self.complement = {'a':'t','t':'a','c':'g','g':'c','n':'n','r':'n','y':'n','s':'n','w':'n','k':'n','m':'n','b':'n','d':'n','h':'n','v':'n','.':'n'}
        self.new = []
        for self.base in seq[::-1]:
            self.base = self.complement[self.base]
            self.new.append(self.base)
        return("".join(self.new))

if __name__ == '__main__':
    Variables()
    open_fai = OpenFile(args.fai1, "fai")
    open_fai = OpenFile(args.fai2, "fai")
    open_map = OpenFile(args.map, "map")


           
