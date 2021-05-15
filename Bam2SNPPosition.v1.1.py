##########################################################
### Import Necessary Modules

import argparse		               #provides options at the command line
import sys		               #take command line arguments and uses it in the script
import gzip		               #allows gzipped files to be read
import re		               #allows regular expressions to be used
import pysam			       #allows sam and bam files to be read
import SNPPosBamv10                       #finds reference sequence and haplotype

##########################################################
### Command-line Arguments
parser = argparse.ArgumentParser(description="A script to read sam or bam alignment files and identifies new SNP position.  Does not output variant position if in the middle of a query insertion.")
parser.add_argument("-file", help = "The location of the sorted bam or sam file, default=<stdin>, comma-separated for multiple files", nargs='+')
parser.add_argument("-pos", help = "The position of the SNP (output, should all be the same, generated using VCF_2_Fasta_v1.1.py), default=101", default=101)
parser.add_argument("-qual", help = "The minimum quality alignment score to output the result, default=10", default=10)
parser.add_argument("-multiple", help = "The number of alternative alignments allowed before rejecting variant, (bwa outputs 1-3), default=1", default=1)
parser.add_argument("-minLen", help = "The minimum length of the alignment, default=150", default=150)
parser.add_argument("-minPid", help = "The minimum percent identity, default=97", default=97)
args = parser.parse_args()

#########################################################
### Opening files and functions
class OpenFile():
    def __init__ (self, f, inde):
        """Opens a file -- sam or bam accepted"""
        if re.search(".sam$", f):
            self.aln = pysam.AlignmentFile(f, "r")
            OpenAln(self.aln, inde)
        elif re.search(".bam$", f):
            self.aln = pysam.AlignmentFile(f, "rb")
            OpenAln(self.aln, inde)

class OpenAln():
    def __init__ (self, aln, ind):
        """Reads the alignment file"""
        #########################################
        ###For loop to iterate over reads in  ###
        #########################################
        self.inter = aln.fetch()
        for self.line in self.inter:
            self.cigar = self.line.cigarstring
            try:
                self.refSeq = self.line.get_reference_sequence()
            except:
                self.refSeq = "NA"
            self.queSeq = self.line.query_sequence
            self.reference = self.line.reference_name 
            self.query = self.line.query_name #contig_1_191:A:T:101
            #self.oldPos = self.query.split(":")[0].split("_")[-1]
            self.alternatives = 0
            self.sPos = self.line.reference_start
            self.ePos = self.line.reference_end
            self.alternatives = []
            self.md = "NA"
            self.mapping_quality = self.line.mapping_quality
            if (self.line.has_tag('XA')):
                ####XA:Z:NW_021816887.1,+4638,55S45M,2; (Chromosome,(strand) position, Cigar, NM(number of mismatches));
                self.alternatives = self.line.get_tag('XA')[:-1].split(";")
            if (self.line.has_tag('MD')):
                ####MD:Z:61T16 or MD:Z:100
                self.md = self.line.get_tag('MD')
            if int(self.mapping_quality) >= int(args.qual) and len(self.alternatives) + 1 <= int(args.multiple):
                (self.referenceSeq, self.haplotype, self.numVar, self.seqAligned, self.newSNPPos) = SNPPosBamv10.Reference(self.queSeq, self.cigar, self.md, self.sPos, int(args.pos))
                if self.newSNPPos != "NA":
                    self.pid = (1 - int(self.numVar)/float(self.seqAligned)) * 100
                    self.newSNPPos = int(self.newSNPPos) + int(self.sPos)
                    #print ("{}\t{}\t{}\t{}\t{}".format(self.referenceSeq, self.haplotype, self.numVar, self.seqAligned, self.newSNPPos))
                    if float(self.pid) >= float(args.minPid) and int(self.seqAligned) >= int(args.minLen):
                        print ("{}\t{}\t{}\t{}\t{}\t{}".format(self.query, args.pos, self.reference, self.newSNPPos, self.seqAligned, self.pid))


#########################################################
###Order of things to be called
if __name__ == '__main__':
    if len(args.file) > 0:
        for index, bamFile in enumerate(args.file):
            open_file = OpenFile(bamFile, index)
