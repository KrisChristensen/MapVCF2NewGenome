##########################################################
### Import Necessary Modules

import sys    #take command line arguments and uses it in the script
import re     #allows regular expressions to be used

"""
This script is meant to be called from another script.  It regenerates the specific position of a nucleotide from one reference to another accounting for indels.
"""

def Reference(seq, cig, md, refS, snpPos):
    mainPos = 0
    refSeq = ""
    refSeq2 = ""
    haplotype = ""
    number_variants = 0
    seq_aligned = 0
    indel_count = 0
    md = md.split(":")[-1]
    snpPosO = int(snpPos)
    ###First the reference sequence is put into the general shape it is supposed to be in
    #print ("\t{}\t{}\t{}\t{}".format(snpPos, cig, md, refS))
    if cig is None:
        return ("NA", "NA", "NA", "NA", "NA")
    for cigar in re.findall('\d+\w', cig, re.I):
        number, character = re.findall('\d+|\w+', cigar, re.I)
        ###Handles Clipping and insertions from the query
        if str(character) == "S" or str(character) == "H" or str(character) == "I":
            if str(character) == "I":
                haplotype += str(int(refS) + int(mainPos)) + "I"
                number_variants += 1
                indel_count += int(number)
                seq_aligned += 1
            #print ("\t\tRemoving\t{} to SNPpos {}, now {}".format(number, snpPos, int(snpPos)-int(number)))            
            if snpPos != "NA" and (int(snpPosO) > int(mainPos) and int(snpPosO) <= int(mainPos) + int(number)):
                snpPos = "NA"            
            elif snpPos != "NA" and int(snpPosO) > int(mainPos) + int(number):
                snpPos = int(snpPos) - int(number)
            mainPos += int(number)
        ###Handles Deletions or Unknown Sequence from the query
        elif str(character) == "D" or str(character) == "N" or str(character) == "P":
            refSeq += int(number)*"N"
            if snpPos != "NA" and int(snpPosO) > int(mainPos):
                #print ("\t\tAdding\t{} to SNPpos {}, now {}, main: {}".format(number, snpPos, int(snpPos)+int(number), mainPos))
                snpPos += int(number)
        ###Handles Matching
        elif str(character) == "M":
            refSeq += seq[int(mainPos) :int(number) + int(mainPos)]
            mainPos += int(number)
            seq_aligned += int(number)
        else:
            sys.stderr.write("\tWarning, {} found in cigar string from alignment file and is not supported\n".format(character))
    #print ("\t\t\tSNPpos now {}".format(snpPos))
    ###Next the reference sequence is modified to output the correct sequence
    mainPos = 0
    for mdfull in re.findall('\d+[a-zA-Z|^][a-zA-Z]*|\d+$', md, re.I):
        if len(re.findall('\d+|\w+|\^\w+', mdfull, re.I)) > 1:
            number, character = re.findall('\d+|\w+|\^\w+', mdfull, re.I)
        else:
            character = "NA"
            number = mdfull
        if re.search("^\^", str(character)):
            character = character[1:]
            refSeq2 += refSeq[int(mainPos) :int(number) + int(mainPos)] + str(character)
            haplotype += str(int(refS) + int(mainPos) + int(number)) + "[" + str(character) + "]"
            mainPos += int(number) + 1
            number_variants += 1
            seq_aligned += 1
        elif str(character) == "NA":
            refSeq2 += refSeq[int(mainPos):]
        else:
            refSeq2 += refSeq[int(mainPos) :int(number) + int(mainPos)] + str(character)
            haplotype += str(int(refS) + int(mainPos) + int(number)) + str(character)
            mainPos += int(number) + 1
            number_variants += 1
    #if snpPos != "NA" and (int(len(seq)) + int(refS) - int(indel_count) < int(snpPos) or int(refS) > int(snpPos)):
    #    snpPos = "NA"
    #elif(snpPos != "NA"):
    #    snpPos = int(snpPos) + int(refS) - 1
    return (refSeq2, haplotype, number_variants, seq_aligned, snpPos)

if __name__ == '__main__':
    seqEx = "ATAGACCAGGACCAAGTTACGATAGGAGGACCCAATTTTACGCAGACCATACATAGCGGACATATTGGACCGAGGTAGGGCCAGAGATATCAGAGTTTAG"
    cigEx = "35M1I25M1D39M"
    mdEx = "MD:Z:41G18^T27T11"
    refStart = "1"
    snpPos = "43"
    refSeqEx, haplotypeF, nVar, sAligned, newSNP = Reference(seqEx, cigEx, mdEx, refStart, snpPos)    
    print ("{}".format(refSeqEx))
    print ("{}".format(seqEx))
    print ("{}".format(haplotypeF))
    print ("{}".format(nVar))
    print ("{}".format(sAligned))
    print ("{}".format(newSNP))


####################################################################
####  Example

#1    97    CHR1    1    60    35M1I25M1D39M    =    111    210    ATAGACCAGGACCAAGTTACGATAGGAGGACCCAATTTTACGCAGACCATACATAGCGGACATATTGGACCGAGGTAGGGCCAGAGATATCAGAGTTTAG    JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ    NM:i:4    MD:Z:41G18^T27T11    MC:Z:27M1D36M1I36M    AS:i:75    XS:i:0

### Example Alignment ###
#         10        20        30         40	     50        60        70        80        90        100 
#12345678901234567890123456789012345 67890123456789012345678901234567890123456789012345678901234567890
#ATAGACCAGGACCAAGTTACGATAGGAGGACCCAA TTTACGGAGACCATACATAGCGGACTATATTGGACCGAGGTAGGGCCAGAGATTTCAGAGTTTAG
#||||||||||||||||||||||||||||||||||| |||||| |||||||||||||||||| ||||||||||||||||||||||||||| |||||||||||
#ATAGACCAGGACCAAGTTACGATAGGAGGACCCAATTTTACGCAGACCATACATAGCGGAC ATATTGGACCGAGGTAGGGCCAGAGATATCAGAGTTTAG
#1234567890123456789012345678901234567890123456789012345678901 234567890123456789012345678901234567890
#         10        20        30        40	    50        60         70        80        90        100

### Actual Reference sequence###
#ATAGACCAGGACCAAGTTACGATAGGAGGACCCAATTTACGGAGACCATACATAGCGGACTATATTGGACCGAGGTAGGGCCAGAGATTTCAGAGTTTAG

### SNP position = 42
