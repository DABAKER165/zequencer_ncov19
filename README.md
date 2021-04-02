# Zequencer Update for SARS2-COVID-19
- Zequencer is a variant calling pipeline that has been optimized for Illumina's MiSeq generated reads
- It assumes: 
    - The reads are that high quality, high accuracy as expected with a properly calibrated MiSeq platform.
    - Input files are in the fastq.gz format 
    - Input files use the were generated with the artic ncov-2019 v3 primers.
    - more information is located at https://artic.network/ncov-2019 
    - File names are of the naming format SAMPLE_NAME-R1.fastq.gz and SAMPLE_NAME-R2.fastq.gz
        - Where R1 is the has one set of reads, and R2 has the matching reverse paired reads.
    - headers are in the miseq nomenclature:
        - example @M01472:366:000000000-JDJFD:1:1101:14487:2478 2:N:0:15 matches to: 
        @M01472:366:000000000-JDJFD:1:1101:14487:2478 1:N:0:15
- Additional Files needed for each run:
     - bed file that gives the start and end location of each primer agains the reference sequence
     - reference sequence fasta of the full genome (reference genome fasta that the bed file mapped)
     - snpEff configuration files that was build with snpEff build
 
 ## Considerations for using this version of Zequencer
- Zequencer has been modified over the years.
- The pipeline one uses is specific to the primers set.
- Considerations to modifying the pipeline include:
    - Are one or more of the primer multiplexed and amplify differently in the overlapping region?
        - The primers are based on the SARS-COV2 - MN908947.3 strain. Some primers have alternates that are offset with the goal of increasing the amplification in that part of the genome
        - The after trimming, consecutive amplicons overlap by 20-60 bp's (see bed file)
    - Are all primers the same length?
        - This workflow assumes the primers are not all the same length. It attempts to trim only the primer region instead of just trimming a set bp from each side.
    - Are the forward primers paired with a specific reverse primer?
        - This worflow assumes the use of forward and reverse primers.  Additionally, each amplicon should map to one location in the genome.
        - If there is recombination, this pipeline may not be able to identify the recombination. The reads are mapped and filtered using the reference genome or expected amplicon regions of the genome
        - A pipeline that specializes in recombination (and possibly a different primer system )
    - Are there alternates in the primers (offset by a few base pairs)?
        - This primer system has alternates that are offset by a few base pairs.  
        - When set to strict_adaptors, the forward and reverse primers much be a perfect match and occur on opposites ends of the reads.
            -  Alternate primers are allowed to be considered a match regaurdless if the reverse primer is alternate.
    - Do primers in the same (or about the same) locations have different sequences meant to amplify different variants?
        - The primers are not multiplexed to pick up variants.  If the downsampling is set to a positive number (200 or 2000 recommended) downsampling normalization can still be done 
    - Do primers map to multiple locations in the genome?
        - The primers do not map to multiple locations in the genome.  
        - This version merges the forward and reverse reads and requires a miniumum overlap of 12 BP to be considered a valid read.  
        - As stated before reads resulting from recombination will likely be removed and not called as variants.  
        - This merge step may make it difficult to find recombination, but it this helps remove PCR Chimera's or other PCR errors, where recombination happens during the PCR process
- Not all pipelines work for all situations.  Do not blindly apply this pipeline to another primer/adapter system. 
   
     
# Arguments
- Defaults are based on the docker deployed version.
- Bed path, ref_adaptor_path, ref_path and config_path should be changed as needed.
- Many of these settings directly change the arguments of the BBMAP suite of tools
- Setting strict_adaptors to TRUE will use a custom script to locate the primers, ensure they are a perfect match, a matching reverse paired end, occur on the ends, and trim the primers.

| ﻿argument               | single_character | type    | default                     | help                                                                                                                                                                                                                                                              |
| ----------------------- | ---------------- | ------- | --------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| sample_dir              | a                | string  |                             | sample_dir where your samples reside must be SAMPLENAME_R1.fastq.gz and SAMPLENAME_R2.fastq.gz format                                                                                                                                                             |
| ref_adaptor_path        | b                | string  | /nCoV-2019_adapters.fasta   | ref_adaptor_path filepath of your fasta adaptor/ primer file (forward only, no rev. comp's) is                                                                                                                                                                    |
| ref_path                | c                | string  | /ref/MN908947.3.fasta       | filepath of the reference sequence                                                                                                                                                                                                                                |
| config_path             | d                | string  | /ref/MN908947.3.config      | filepath of .config file to generate snpeff                                                                                                                                                                                                                       |
| ram                     | e                | int     | 4000                        | set ram explicitly for bbmap suite (in MB 4000 = 4GB)                                                                                                                                                                                                             |
| minallelefraction       | f                | float   | 0.1                         | fraction of called variant at that position (Allele Depth / Total Depth)                                                                                                                                                                                          |
| minlength               | g                | int     | 250                         | minimum merged read length. Shorter reads will be removed                                                                                                                                                                                                         |
| maxlength               | h                | int     | 450                         | maximum merged read length.  Longer reads will be removed.                                                                                                                                                                                                        |
| ncbi_accession          | i                | string  | MN908947.3                  | Ncbi_accession name, must match the snp_eff name in config                                                                                                                                                                                                        |
| minlength_diff          | j                | int     | 60                          | bp lower limit of merged read length compared to  the mapped inner primer region length,  strict_adapters must be set to True (t) to use                                                                                                                          |
| maxlength_diff          | k                | int     | 24                          | bp upper limit of merged read length compared to  the mapped inner primer region length,  strict_adapters must be set to True (t) to use                                                                                                                          |
| mincov                  | l                | int     | 2                           | minimum coverage at that position                                                                                                                                                                                                                                 |
| minreads                | m                | int     | 10                          | minimum read count of the allele/variant of interest                                                                                                                                                                                                              |
| callsub                 | n                | boolean | t                           | Set to false (f) to remove all substitutions from the variant file                                                                                                                                                                                                |
| callins                 | o                | boolean | t                           | Set to false (f) to remove all insertions from the variant file                                                                                                                                                                                                   |
| calldel                 | p                | boolean | t                           | Set to false (f) to remove all insertions from the variant file                                                                                                                                                                                                   |
| strict_adaptors         | q                | boolean | f                           | Set to false (f) to use bbduk to filter. Set to True (looks for first character = T or t) to use a custom python script that ensures merged reads are of appropriate length and the reverse paired read primer matches the forward read primer                    |
| minscore                | r                | int     | 10                          | minimum score in the read (not average)                                                                                                                                                                                                                           |
| minquality              | s                | int     | 10                          | minimum quality reads with lower value will be remboed                                                                                                                                                                                                            |
| rarity                  | t                | float   | 0.1                         | set to the same value as minallele fraction                                                                                                                                                                                                                       |
| coverage                | u                | boolean | t                           | set to true to use coverage scores.                                                                                                                                                                                                                               |
| maxindel                | v                | int     | 100                         | maximum insertion (N's) a read will be allowed to add and still be considered a valid mapping.                                                                                                                                                                    |
| k                       | w                | int     | 21                          | bbduk setting                                                                                                                                                                                                                                                     |
| ktrim                   | x                | string  | l                           | bbduk setting (r, l)                                                                                                                                                                                                                                              |
| qtrim                   | y                | string  | rl                          | bbduk setting (r, l, rl)                                                                                                                                                                                                                                          |
| trimq                   | z                | int     | 30                          | bbduk setting                                                                                                                                                                                                                                                     |
| hdist                   | A                | int     | 3                           | bbduk setting                                                                                                                                                                                                                                                     |
| rcomp                   | B                | boolean | f                           | bbduk setting                                                                                                                                                                                                                                                     |
| minlength_bbduk         | C                | int     | 75                          | bbduk setting                                                                                                                                                                                                                                                     |
| restrictleft            | D                | int     | 32                          | bbduk setting                                                                                                                                                                                                                                                     |
| minavgquality           | E                | int     | 20                          | bbduk setting                                                                                                                                                                                                                                                     |
| minbasequality          | F                | int     | 20                          | bbduk setting                                                                                                                                                                                                                                                     |
| use_diff                | G                | boolean | t                           | Use a +/- of each amplicon as set by minlength_diff and maxlength_diff instead of globally declaring minlength/maxlength regaurdless of the amplicon of interet.  Only works if strict_adaptors is set                                                            |
| bed_path                | H                | string  | /nCoV-2019/V3/nCoV-2019.bed | File path of the bed file.  Must be set to use strict_adaptors.                                                                                                                                                                                                   |
| outdir                  | I                | string  | ../out                      | if you want to manually set the outdir set it                                                                                                                                                                                                                     |
| reads_per_amplicon      | J                | int     | 0                           | used for downsampling normalization 0 will not down sample.  200 or 2000 is typical                                                                                                                                                                               |
| downsampling_fasta_path | K                | string  |                             | fasta of all the amplicons (start of forward to end of reverse primer to the reference genome given.  If your amplicons map to multiple references you cannot use the short cut.  If it is not declared it can be generated from the bedfile and reference fasta. |
| fast_normalize          | L                | boolean | t                           | uses a customized, faster normalization, but can only be used in conjunction with strict_adaptors = True, if using normalization this technique is upto 100x faster.                                                                                              | 
| rem_int_files           | M                | boolean | t                           | Removes all intermediate files created in the pipeline leaving the sorted bam, indexed bam, vcf and annotated vcf files                                                                                                                                           |

# Installation:
## Method 1 Local install:
- There is a lot of prerequisites you will need the following packages, directories and files:
This pipeline has been successfully tested on the following operating systems:
MacOSX 10.14 with some packages installed with homebrew.
Ubuntu 18.04, 64 bit
Debian (Buster 10.8, 64 bit)
It is not compatable with a native Windows operating system (see docker instructions for Windows's operating systems)


```
BBMap_38.90.tar.gz
bcftools 
samtools
htslib-1.12 (bgzip, tabix)
VarScan.v2.4.4.jar
snpEff_latest_core
seqtk-1.3
openjdk=11.0.1
python3 (3.5 to 3.7)
snakmake (for future workflow)
biopython 1.76
pyfasta
pandas
```
Additional files from this github repository (for COVID-19)
trim_strict.py
ncov-2019_adapters.fasta
zequencer.sh
nCoV-2019 -- directory and contents
ref -- directory and contents
## Method 2 docker:
This has not been tested with docker on a Windows platform.
This pipeline has been successfully tested on the following operating systems:
MacOSX 10.14 with some packages installed with homebrew.
Ubuntu 18.04, 64 bit

Download the directory zequencer_docker from github.
Decompress the directory as needed, and change to the directory you downloaded.
example if downloaded to you home directory ~/:
```bash
cd ~/zequencer_docker
# zequencer_ncov:v1 can be changed to a different tag as needed
docker build -t zequencer_ncov:v1 .
# run docker image (this uses in interactive mode, but interactive mode is not required), 
# mounted volumes with -v can be changed as needed
# cpus can be changed as needed
docker run --cpus 4 -it -v /Volumes:/Volumes -v /Users:/Users zequencer_ncov:v1

```
# Example usage:
- Volumes have been mounted using the above Docker code
- In this example, the sample_dir is "/Volumes/T7/zequencer/south_africa/samples"
- Change this path to reflect your sample directory of paired fastq.gz files
- The paired files must have nomenclature of: SAMPLE_NAME-R1.fastq.gz and SAMPLE_NAME-R2.fastq.gz
- I am using the default relative path (that is created if it does not exist of ../out )
    - Which for this case is: /Volumes/T7/zequencer/south_africa/out
```bash
# downsample to 2000 reads per amplicon, 
# when calling variants use: allele fraction of 0.01 
# and use BBDUK
 /zequencer.sh \
 --sample_dir /Volumes/T7/zequencer/south_africa/samples \
 --minallelefraction 0.01 \
 --rarity 0.01 \
 --strict_adaptors f \
 --fast_normalize f \
 --reads_per_amplicon 2000
 
 # Do not downsample (much faster),reads per amplicon, 
 # when calling variants use:allele fraction of 0.01 
 # and use BBDUK
 # you can also set reads_per_amplicon to 0 instead of blank.
 /zequencer.sh \
 --sample_dir /Volumes/T7/zequencer/south_africa/samples \
 --minallelefraction 0.01 \
 --rarity 0.01 \
 --strict_adaptors f \
 --fast_normalize f

 # Downsample to 2000 reads per amplicon(,reads per amplicon,
 # when calling variants use:  allele fraction of 0.01 
 # filter using only exact and paired matches of primers that are at the end of reads
 /zequencer.sh \
 --sample_dir /Volumes/T7/zequencer/south_africa/samples \
 --minallelefraction 0.01 \
 --rarity 0.01 \
 --strict_adaptors t \
 --fast_normalize f \
 --reads_per_amplicon 2000
 
 # Downsample to 2000 reads per amplicon(,reads per amplicon,
 # when calling variants use:  allele fraction of 0.01 
 # filter using only exact and paired matches of primers that are at the end of reads
 # use a faster upto 100x algorithm (safe for the ncov2019, V3 primer set).
 /zequencer.sh \
 --sample_dir /Volumes/T7/zequencer/south_africa/samples \
 --minallelefraction 0.01 \
 --rarity 0.01 \
 --strict_adaptors t \
 --fast_normalize t \
 --reads_per_amplicon 2000
 
 
 # when calling variants use:  allele fraction of 0.01 
 # filter using only exact and paired matches of primers that are at the end of reads
 # Fastest option, but does not downsample.
 /zequencer.sh \
 --sample_dir /Volumes/T7/zequencer/south_africa/samples \
 --minallelefraction 0.01 \
 --rarity 0.01 \
 --strict_adaptors t
```

# About fast_normalize downsampling normalization
- The existing nomalization was developed by Dr. Nick Loman, Professor of Microbial Genomics and Bioinformatics at  University of Birmingham, and adapted by Dr. David O'Connors lab (fast_normalize = false) maps each read to each primer (all 196+)
    - Then it will only retain up to the amount defined in reads_per_adaptor per each mapped primer
    - Example: 
        - if reads_per_adaptor=2000, and 6000 reads map to primer one and 1500 reads to primer 2,
        - only 2000 of the 6000 reads will be retained for primer one and all 1500 reads mapping to primer 2. 2.
    - This is particularly useful if you have variants the only map to one primer set that is overlapping with another primer set.
- The fast_normalize method also downsamples, but does only one mapping step (instead of one per primer), and uses the position the start position of the mapped read to the reference and the positions of the primers in the bed file
    - This uses a custom python script and pysam to conduct the filtering.
    - This has been shown to be upto 100x faster than the other downsampling method
    - This mimics the algorithm used in Nexstrain's Minion step for ONT downsampling.
    - The main advantage is speed, however if a future primer set is developed with primers that are specialized to find variants, this technique is not adequate for downsampling.
    - Generally downsampling is optional, as the primer sets use the same reference sequence and is more to try to get a more accurate percentage for a variant in the case of one primer in an overlapping region doing better to amplify reads than another primer.
# Trouble shooting:
- set rem_int_files to false (--rem_int_files f) to keep the intermediate files generated by the pipeline.  This will require about 5-10x more hardrive space but may help figure out issues
- During testing docker has randomly lost the ability to use relative paths because cwd or pwd failed on the  mounted Volumes.  
    - You will see this as a "pwd cannot find ." or "os.getcwd()" error. 
    - using the command "pw"d will result in an error in bash.
    - Exiting and rerunning the docker image and/or restarting docker resolved this issue  during testing
# Citations
- BBMAP
    - Bushnell, Brian. BBMap: A Fast, Accurate, Splice-Aware Aligner. United States: N. p., 2014.
- samtools, bcftools, htslib
    - Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools.
- VarScan
    - VarScan: variant detection in massively parallel sequencing of individual and pooled samples
- Daniel C. Koboldt, Ken Chen, Todd Wylie, David E. Larson, Michael D. McLellan, Elaine R. Mardis, George M. Weinstock, Richard K. Wilson, Li Ding
    - Bioinformatics. 2009 Sep 1; 25(17): 2283–2285.
- SnpEff
    - "A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.", Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (Austin). 2012 Apr-Jun;6(2):80-92. PMID: 22728672
- seqtk
    - https://github.com/lh3/seqtk 
- pyfasta
    - https://pypi.org/project/pyfasta/
- biopython
    - Peter J. A. Cock, Tiago Antao, Jeffrey T. Chang, Brad A. Chapman, Cymon J. Cox, Andrew Dalke, Iddo Friedberg, Thomas Hamelryck, Frank Kauff, Bartek Wilczynski, Michiel J. L. de Hoon, Biopython: freely available Python tools for computational molecular biology and bioinformatics, Bioinformatics, Volume 25, Issue 11, 1 June 2009, Pages 1422–1423
- pandas
    - McKinney, W., & others. (2010). Data structures for statistical computing in python. In Proceedings of the 9th Python in Science Conference (Vol. 445, pp. 51–56).
- artic network and nextstrain
	- https://artic.network/ncov-2019
	- https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V3
# How to cite
- use url https://github.com/dabaker165/zequencer_ncov19

# Acknowledgements
- David A. Baker bioinformatician in Department of Pathology and Laboratory Medicine at the University of Wisconsin lead programmer of this pipeline
- Dr. David O'Connor,  Medical Foundation Professor, Department of Pathology and Laboratory Medicine at the University of Wisconsin for writing the first Zequencer Pipeline,
- Dr. Nick Loman, Professor of Microbial Genomics and Bioinformatics at  University of Birmingham and Shelby O'Connor, Associate Professor, Department of Pathology and Laboratory Medicine at the University of Wisconsin 
- John (JJ) Baczenas: Associate Research Specialist in Department of Pathology and Laboratory Medicine at the University off Wisconsin for defining Pipeline criteria, generating test data, and QA testing.
- Katarina (Kat) Braun, MD-PhD, and Gage Moreno PhD students in Department of Pathology and Laboratory Medicine at the University of Wisconsin for defining pipeline criteria and generating test data.


# Roadmap items (future)
- Add Snakemake capability to manage threads and process multiple files concurrently
