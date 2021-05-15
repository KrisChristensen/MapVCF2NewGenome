##########################################################
### Import Necessary Modules

import argparse                       #provides options at the command line
import sys                       #take command line arguments and uses it in the script
import gzip                       #allows gzipped files to be read
import re                       #allows regular expressions to be used

##########################################################
### Command-line Arguments
parser = argparse.ArgumentParser(description="Output filtered repositioned SNPs of a vcf file based on supplied positions.")
parser.add_argument("-map", help = "The location of the file output by Bam2SNPPosition.v1.0.py with the old positions and new positions", default=sys.stdin, required=True)
parser.add_argument("-vcf", help = "The location of a vcf file to replace snp locations", default=sys.stdin, required=True) 
parser.add_argument("-fai", help = "The location of a fai file (samtools faidx) for the genome mapping to (for use in generating the vcf file contig positions)", default=sys.stdin, required=True)    
args = parser.parse_args()

#########################################################
### Open file (object-oriented programming)
class Variables():
    loci = {}
    contigs = {}

class OpenFile():
    ### Opens the file and directs it to a line-reader for the specific file-type
    def __init__ (self, filename, fileType):
        if re.search(".gz$", filename):
            self.filename = gzip.open(filename, 'rb')
        else:
            self.filename = open(filename, 'r')             
        if fileType == "vcf":
            sys.stderr.write ("Opened vcf file: {}\n\n".format (filename))
            self.readLinesVcf()
        elif fileType == "map":
            sys.stderr.write ("Opened map file: {}\n\n".format (filename))
            self.readLinesMap()
        elif fileType == "fai":
            sys.stderr.write ("Opened fai file: {}\n\n".format (filename))
            self.readLinesFai()        

    ### Reads alignment files
    def readLinesMap(self):
        for line in self.filename:
            try:
                line = line.decode('utf-8')
            except:
                pass                
            line = line.rstrip('\n')
            self.original, self.originalPosition, self.new, self.newPosition, self.len, self.pid = line.split()[0:6]
            Variables.loci["{}\t{}".format(self.original, self.originalPosition)] = "{}\t{}".format(self.new, self.newPosition)
        self.filename.close()

    ### Reads alignment files
    def readLinesFai(self):
        for line in self.filename:
            try:
                line = line.decode('utf-8')
            except:
                pass                
            line = line.rstrip('\n')
            self.contig, self.length, self.other1, self.other2, self.other3 = line.split()
            Variables.contigs[self.contig] = self.length
        self.filename.close()

    ### Reads vcf file and return remapped SNPs
    def readLinesVcf(self):
        self.contigPrint = 0
        for line in self.filename:
            try:
                line = line.decode('utf-8')
            except:
                pass
            line = line.rstrip('\n')
            if not re.search("^#", line):
                line_info = line.split("\t")
                self.chr = line_info[0]
                self.pos = line_info[1]
                if "{}\t{}".format(self.chr, self.pos) in Variables.loci:
                    self.newChr, self.newPos = Variables.loci["{}\t{}".format(self.chr, self.pos)].split("\t")
                    line_info[0] = self.newChr
                    line_info[1] = self.newPos
                    print("{}".format("\t".join(line_info)))
            else:
                if re.search("^##fileformat", line) or re.search("^##ALT=", line) or re.search("^##FILTER=", line) or re.search("^##FORMAT=", line) or re.search("^##INFO=", line) or re.search("^#CHROM", line):
                    print ("{}".format(line))
                elif re.search("^##contig=", line) and self.contigPrint == 0:
                    self.contigPrint = 1
                    for self.contig in Variables.contigs:
                        print ("##contig=<ID={},length={}>".format(self.contig, Variables.contigs[self.contig]))
        self.filename.close()
       
if __name__ == '__main__':
    Variables()
    open_fai = OpenFile(args.fai, "fai")
    open_map = OpenFile(args.map, "map")
    open_vcf = OpenFile(args.vcf, "vcf")
