##########################################################
### Import Necessary Modules

import argparse                       #provides options at the command line
import sys                       #take command line arguments and uses it in the script
import gzip                       #allows gzipped files to be read
import re                       #allows regular expressions to be used

##########################################################
### Command-line Arguments
parser = argparse.ArgumentParser(description="Output filtered repositioned SNPs of a vcf file based on supplied positions.")
parser.add_argument("-map", help = "The location of the file output by Bam2SNPPosition.v1.0.py and verified with CheckMappingVCF.v1.1.py with the old positions and new positions", default=sys.stdin, required=True)
parser.add_argument("-vcf", help = "The location of a vcf file to replace snp locations", default=sys.stdin, required=True)
parser.add_argument("-gen", help = "The full path of the new genome (for use in generating the new vcf file in proper format)", default=sys.stdin, required=True)
parser.add_argument("-fai", help = "The location of a fai file (samtools faidx) for the genome mapping to (for use in generating the vcf file contig positions)", default=sys.stdin, required=True)    
args = parser.parse_args()

#########################################################
### Open file (object-oriented programming)
class Variables():
    loci = {}
    contigs = {}
    duplicates = {}

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
            #previous format
            #contig_4431	451779	CM029847.1	60662	201	99.50248756218906	most	a	c
            #new format
            #NC_042551.1_4892465:A:G:51      51      NC_056429.1     13439842        101     98.01980198019803       most
            #NC_042561.1_60157999:A:G:51     51      NC_056429.1     13691221        102     94.11764705882352       exact
            self.originalFull, self.subsetPos, self.new, self.newPosition, self.len, self.pid, self.validation = line.split()[0:7]
            self.contig1andPos = self.originalFull.split(":")[0]
            self.original = "_".join(self.contig1andPos.split("_")[:-1])
            self.originalPosition = self.contig1andPos.split("_")[-1]
            self.previousReference = self.originalFull.split(":")[1]
            Variables.loci["{}\t{}".format(self.original, self.originalPosition)] = "{}\t{}".format(self.new, self.newPosition)
            ###Check to make sure there are no duplicates and if there are removes all of them later
            if "{}\t{}".format(self.new, self.newPosition) in Variables.duplicates:
                Variables.duplicates["{}\t{}".format(self.new, self.newPosition)] += 1
                sys.stderr.write("\tWarning: tried mapping, position {} {} to {} {}, but that position has already been mapped to.  Removing all instances of this position.\n\n".format(self.original, self.originalPosition, self.new, self.newPosition))
            else:
                Variables.duplicates["{}\t{}".format(self.new, self.newPosition)] = 1
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
                self.refAllele = line_info[3]
                self.altAllele = line_info[4]
                ###Check to make sure this position had a corresponding location in the new genome
                if "{}\t{}".format(self.chr, self.pos) in Variables.loci:
                    self.newChr, self.newPos = Variables.loci["{}\t{}".format(self.chr, self.pos)].split("\t")
                    ###Check to see if new position is uniquely mapped (doesn't use otherwise)
                    if "{}\t{}".format(self.newChr, self.newPos) in Variables.duplicates and Variables.duplicates["{}\t{}".format(self.newChr, self.newPos)] == 1:
                        #self.allelesMatch, self.newRefAllele, self.newAltAllele = self.checkAlleles(self.chr, self.pos, self.refAllele, self.altAllele)
                        ###Check to verify that the alleles match the new genome
                        #if self.allelesMatch == 1:
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
                elif re.search("^##reference=", line):
                    print ("##reference=file://{}".format(args.gen))
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
    open_fai = OpenFile(args.fai, "fai")
    open_map = OpenFile(args.map, "map")
    open_vcf = OpenFile(args.vcf, "vcf")
