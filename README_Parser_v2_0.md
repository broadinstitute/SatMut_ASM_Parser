# SatMut_ASM_Parser

**INTRODUCTION**

Saturation Mutagenesis screen (and library QC sequencing) data are complex and gigantic. A typical screen produces 3-5 TB sequencing data in fastq format. ASMv1.0, developed by Ted Sharpe of The Broad Institute, can process the fastq files to identify all (planned and not planned) 'variants' and assign the 'counts' to the detected 'variants.' ASMv1.0 produces many output files for various purposes. Each next-gen sequencing index pair (typically one experimental sample) gets one set of output files. Two ASMv1.0 output files, variantCounts and codonCounts, are subjected to this parser. The codonCounts files are a matrix of fixed dimension, where ncol = 64 (i.e. 64 codons), nrow = 'aa length of the ORF.' The variantCounts files are huge files, where each row represents a species that was detected in read-pairs. It is not uncommon to see near million of rows in a variantCounts file. 

It is challenging for non-data-scientists to lay their hands on either of the ASMv1.0 output files: 

- codonCounts files: one matrix per sample. At issue are: How to convert the matrix into variants and variant counts? How to tease apart ’planned' and 'not planned' variants? How to remove wildtype codon counts? How/Whether to scale the data?  

- variantCounts files: one giant file per sample and variants are presented in form of nucleotide variations in relation to the reference sequence. How to convert such presentations of the variants into codon changes? amino acid changes? 'planned' vs 'not planned?' among the 'not planned' if part of its changes is a 'planned' change? frame-shift or not? substitution or insertion or deletion? and so on so forth.

This parser is to fill the gap. It is written in R. It takes two kinds of output files of ASMv1.0, codonCounts and variantÇounts, and parse/aggregate/merge/annotate them into csv data files that can be digested by non-data-scientists.

**What is ASMv1.0?**

ASMv1.0, developed by Ted Sharpe of The Broad Institute (tsharpe@broadinstitute.org), takes NGS .fastq files to assign counts to the detected variants.

*READS FILTERING*: The reads from an NGS fragment are trimmed (1) by base quality (from 2 ends, finding first 15-bp that every base has above base-quality threshold), (2) by removing non-ORF transposon sequences, (3) by requiring that the outermost nucleotide variations in the reads have a number of reference bases as specified by --min-flanking-length parameter (2-5, depending on how the library is designed, e.g. if the library has all single codon changes, --min-flanking-length=2 will suffice to ensure a codon-signature of a variant is not split between two Nextera fragments), and (4) by requiring the trimmed reads be longer than the –min-length = 40, for example.

*VARIANT CALLING*: The reads that have passed the filtering are used differently by ORFcall (an earlier version software also by Ted Sharpe, as used in: Giacomelli et al (2018) Nat Genet 50: 1381–1387) and ASMv1.0:

(1) Under ORFcall, variants are called in context of a single codon space, that is, if a read, or read pair captures multiple codon variations, they are called as multiple independent molecules. This seems not fair as the read counts are assigned to each of the molecules represented by each codon change, without regard of the association of the multiple changes physically in a single molecule. Yet since a good quality of a library generally don’t have many unplanned variations, the codon-centric calls by ORFcall is in fact a good proximation. 

(2) Under the molecule-centric AnalyzeSaturationMutagenesis, the NGS reads are analyzed in its entirety, and with following new capabilities:
 
-	The trimmed reads are evaluated full-length, that is, variants are called in the context of entire read (or read pair). The programmed variant is called when 2 conditions are met: the detection of the programmed codon changes, AND the absence of any ‘additional’ nucleotide variations throughout the entire read or read pair. 

-	Indel detection: This feature was lacking in ORFcall. This new feature of ASMv1.0 allows us to process indel libraries with up to 10-codon indels.

**What is ASM_Parser?**

Written by Xiaoping Yang of The Broad Institute (xyang@broadinstitute.org), ASM_Parser produces 3 Data Files:

- ASMv1.0 reports .codonCounts as ‘back-compatible ORFcall’ results (similar to but not the same as the original ORFcall, due to difference in READS FILTERING). Same as the ORFcall’s .cdn reports, users of  these results need to (1) transform the data files from a table per sample into a column per sample (e.g. using R melt function), (2) aggregate (e.g. using R merge function) data into a single data file, a column per sample, (3) scaling to mitigate Nextera sequencing coverage biases. 

File 0: Data File of variants called by ‘ORFCall’

- The full benefit of ASMv1.0 lies in the file set .variantCounts, one file per sample. The .variantCounts lists the counts of the said variants, one variant per line. The detected variants include ‘substitution’, ‘insertion’, ‘deletion’ relative to the ‘reference file’ that is the ORF template with 2 flanks (~150bp each); the reference file represents the PCR amplicon when the ORF was extracted from genomic DNA. Users need to parse the .variantCounts files into forms that wet-lab biologists can consume. Below are 2 files from ASM_Parser data processing:

File 1: sub-dataset representing the ‘intended’ variants, with following columns:

- 11 leading columns that annotate the detected variants:

11 Lead-Column Names		|                  Example_variant|
--------------------------------------------------------|---------------------------------------|
  [1] "variant.by.nt.at.codon.level"  |     GCG\|2\|CCC |                      
  [2] "variant.by.aa"                 |            A2P   |           
  [3] "POS"                            |              2 |         
  [4] "Wt_codon"                        |         GCG   |             
  [5] "Vt_codon"                         |         CCC  |             
  [6] "Wt_aa"                             |          A   |          
  [7] "Vt_aa"                              |          P   |       
  [8] "delta_nt"                            |         2    |        
  [9] "is.the.whole.mol.planned"        |  Intended         |             
 [10] "del_ins_or_sub"                   |     sub           |        
 [11] "mutationType" 		     | missense|


- The above 11 columns are followed by sample specific 3-column sets: 
per screen sample, (1) raw count, (2) reference count, and (3) scaled abundance. 

The ‘raw count’ column, for example ‘Replicate1_Day0.ct’,  is for the count of the described variant in the said sample, Replicate1_Day0

The ‘ref count’ column, ‘Replicate1_Day0.ref_ct’, is the count sum of all molecules that span the same region that defined the said variant in the said sample.

The ‘fraction’ column,  ‘ Replicate1_Day0.ct_frctn’, if the scaled fraction abundance of the said variant in the said sample. The ‘fraction’ columns can be conveniently used to combine/compare replicates, compare screen arms, and thereby call screen hits. 

If a clonal plasmid sample is processed in same manner as screen samples, the ‘fraction’ value of the clonal sample can be used as the correction factor (e.g. PCR/NGS errors) to subtract from the ‘fraction’ of every experimental sample. Such correction has to be done with caution - over-correction can be worse than no correction.

With this File1, collaborators should be able to plot one sample against the other in forms of: 

- scatter plots: where each dot is a variant and axises are conditions, and in some fancy ways, one can color the data points, for example, by mutation type ’silent’, ‘nonsense’, or ‘missense’  so on so forth…  

- heat maps: in 2d-grids of ‘amino acid positions’ by ‘amino acids’ with ‘fold changes’ as the value.

File 2: The version of data file (with following 23 lead columns) would allow more in-depth analysis.

23 Lead-Column Names							             |                  Example_variant|
--------------------------------------------------------|---------------------------------------|
[1] "variant.description.by.unit.of.wtCdn.cdnPos" |              	ATG\|2\|GAC,CTG\|43\|GTG |
[2] "variant.description.by.unit.of.wt.nt_Ted"       |          	131:A>G, 132:T>A, 133:G>C, 254:C>G|
[3] "variant.description.by.unit.of.cdnPos"                 |   	2\|GAC,43\|GTG  |
[4] "num.whole.mol.nt.changes_Ted"                         |    	4|
[5] "num.whole.mol.cdn.changes_upto.a.stop_Ted"           |     	2|
[6] "list.whole.mol.all.cdn.changes_upto.a.stop_Ted"     |      	2:ATG>GAC, 43:CTG>GTG|
[7] "list.whole.mol.all.aa.changes_upto.a.stop_Ted"     |       	M:M>D, M:L>V|
[8] "is.the.whole.mol.planned"                        |         	Unintended|
[9] "pos.of.cdnChange.if.whole.mol.is.intended"   |             	|
[10] "del_ins_or_sub"                               |            sub|
[11] "mutationType"                                  |           missense|
[12] "is.there.frame.shift"                           |          FALSE|
[13] "subList.all.intended.cdnChanges.ifAny.detected.in.the.mol"| 2\|GAC,43\|GTG |
[14] "pick.the.highest.intendedCdnDeltaNtChange.cdnPos"|         2|
[15] "pick.the.highest.intendedCdnDeltaNtChange.wtCdn"  |        ATG|
[16] "pick.the.highest.intendedCdnDeltaNtChange.wtAA"    |       M|
[17] "pick.the.highest.intendedCdnDeltaNtChange.vtCdn"    |      GAC|
[18] "pick.the.highest.intendedCdnDeltaNtChange.vtAA"      |     D|
[19] "pick.the.highest.intendedCdnDeltaNtChange.deltaNT"    |    3|
[20] "num.insertion.nt"                                      |   0|
[21] "num.deletion.nt"                                        |  0|
[22] "the.first.mutation.event.occurs.at.ntPos"                | 4|
[23] "the.first.mutation.event.occurs.at.nthNT.in.a.cdn"    	|	1|



With File2, one can do everything that can be done with data File1, and more. Personally, I would like to assess library purity (the pristine intended molecules vs. those with unwanted changes), and the artificial variant calls from PCR/NGS errors. 

**INSTRUCTION TO RUN THE PARSER**
It saves a lot of trouble by setting up the R work directory in such that in all we will have 2 folders, and 3 supporting files including the R code itself (see the folder/file structure at the end of this document).

*TWO FOLDERS:*
1. Move all .variantCounts files into a folder, e.g. variantCounts/ 
	The individual .variantCounts flies should have been named to have ’Sample[0-9][0-9]’ in the file names. When you ran ASMv1.0, the setup of it should already be good to name the ASMv1.0 output files. If not, you can always rename these files to have  ‘Sample[0-9][0-9]’ in the name.

2. Move all .codonCounts files into a folder, e.g. codonCounts/ The individual .codonCounts flies should have been names to have ’Sample[0-9][0-9]’ in the file names.

*THREE FILES:*
1. “SampleAnnot.csv" file. This is a sample mapping file: two columns named 'Sample' and 'Experiment'. Be advised that 'Sample' column has to be the below format (case sensitive), and do not have 'spaces', or 'dashes' in 'Experiment.’ Note if you have multiple sequencing runs from the same sample ane/or you wish to combine the data from runs of sequencing data into one, you may name those runs with the same Sample number but followed by '_1', '_2' etc. (e.g., Sample08_1) including in both .codonCounts and .variantCounts data file names and the "SampleAnnot.csv" file, AND the 'Experiment' should be the same entry for all samples involved. If you would like to keep these replicates separate, as I often do, name the sample/experiment uniquely.

For example:
Sample	Experiment
Sample01	experiment1
Sample02	experiment2
Sample03	experiment3
Sample04	experiment4
Sample05	experiment5
Sample06	experiment6
Sample07	experiment7
Sample08_1	experiment8
Sample08_2	experiment8
Sample08_3	experiment8
 
2. "Intended_codon_list_1col.csv" file. This is a 1-column csv file listing the planned codon changes, with the column named 'key,’ formatted as aa_position|planned_codon.

For example
key
1|AAA
1|AAT
1|ACT

3. The file the R-code parser itself. You need to edit the top portion of the code to reflect your work. 
-  Wild-type ORF sequence 
-  Screen name
-  Gene name
-  Sample## for pDNA library
-  Sample## for clonal wild-type plasmid
-  Sample## of those un-slected, or less selected that may be used to filter out the called, but low-abundance variants. 

Then RUN the entire .r code. The data files and plots will be written into the 'outbox' folders inside of your ’variantCounts’  or ‘codonCounts’ folders. 

As a training set, in the zip file ‘WorkingDir_for_Parser_v2_0.zip,’  real data and the folder/file structures are already properly set up. You only need to update the ‘work.dir’  variable in the R code ‘SatMut_ASM_Parser_v2_0.R’ to the actual directory on your machine and then run the code.

****************

Good luck and direct your questions to xyang@broadinstitute.org



