# MapVCF2NewGenome
A pipeline for mapping SNPs in a VCF file to a genome that wasn't used to call the variants.  It is expected that this strategy could be used for mapping VCF files to an updated genome assembly or for mapping to another closely related species' genome (to be able to directly compare nucleotide variants).  Right now this only works for SNPs, but could be modified for other variants.


<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#requirements">Requirements</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
  </ol>
</details>

<!-- requirements -->
## Requirements

These scripts have been tested with Python 2.7 and 3 and should work with either.
The script requires the following programs and files.

Programs:<br /><br />
&nbsp;&nbsp;&nbsp;bwa (or equivalent alignment program to produce .bam files)<br />
&nbsp;&nbsp;&nbsp;samtools<br />
&nbsp;&nbsp;&nbsp;bcftools<br />
&nbsp;&nbsp;&nbsp;python<br />
&nbsp;&nbsp;&nbsp;pysam (python)
    
Files:<br /><br />
&nbsp;&nbsp;&nbsp;vcf file that you wish to map to a new genome<br />
&nbsp;&nbsp;&nbsp;old genome file (fasta format and indexed -- samtools faidx oldGenome.fasta, bwa index oldGenome.fasta)<br />
&nbsp;&nbsp;&nbsp;new genome file (fasta format and indexed -- samtools faidx newGenome.fasta)<br />

<!-- usage -->
## Usage

1) Generate flanking sequences around the variant(s) so that they can be mapped to the new genome:<br /><br />
&nbsp;&nbsp;&nbsp;python VCF_2_Fasta_v1.1.py -vcf old.vcf.gz -fasta oldGenome.fasta -fai oldGenome.fasta.fai -flanking 100 > flankingSeqs.fasta<br /><br />
&nbsp;&nbsp;&nbsp;help (and further explanations): python VCF_2_Fasta_v1.1.py -h
    
2) Map fasta sequences to new genome:<br /><br />
&nbsp;&nbsp;&nbsp;bwa mem -M newGenome.fasta flankingSeqs.fasta | samtools sort -o newGenome.v.flankingSeqs.bam<br />
    
3) Index alignment file:<br /><br />
&nbsp;&nbsp;&nbsp;samtools index newGenome.v.flankingSeqs.bam<br />
    
4) Identify variant positions in new genome:<br /><br />
&nbsp;&nbsp;&nbsp;python Bam2SNPPosition.v1.1.py -file newGenome.v.flankingSeqs.bam -qual 60 -multiple 1 -minLen 200 -minPid 99 -pos 101 > old2new.positions.txt<br /><br />
&nbsp;&nbsp;&nbsp;help (and further explanations): python Bam2SNPPosition.v1.1.py -h<br />
&nbsp;&nbsp;&nbsp;note: requires SNPPosBamv10 in same working directory
    
5) Quality control:<br /><br />
&nbsp;&nbsp;&nbsp;python CheckMappingVCF.v1.0.py -map old2new.positions.txt -fasta1 oldGenome.fasta -fai1 oldGenome.fasta.fai -fasta2 newGenome.fasta -fai2 newGenome.fasta.fai -flanking 20 -output most > old2new.positions.filtered.txt<br /><br />
&nbsp;&nbsp;&nbsp;help (and further explanations): python CheckMappingVCF.v1.0.py -h<br />

6) Generate new vcf file:<br /><br />
&nbsp;&nbsp;&nbsp;python Assign_New_SNP_Positions.v1.1.py -map old2new.positions.filtered.txt -vcf old.vcf.gz -fai newGenome.fasta.fai > new.vcf<br /><br />
&nbsp;&nbsp;&nbsp;help (and further explanations): python Assign_New_SNP_Positions.v1.1.py -h<br />

7) Sort vcf file<br /><br />
&nbsp;&nbsp;&nbsp;bcftools sort new.vcf > new.sorted.vcf<br />
    

<!-- license -->
## License 

Distributed under the MIT License.
