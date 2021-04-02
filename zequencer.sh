#! /bin/bash

echo $1

echo $@

usage () { echo "Usage : $0
--sample_dir -a (string) default: <blank>   | Help sample_dir where your samples reside must be SAMPLENAME_R1.fastq.gz and SAMPLENAME_R2.fastq.gz format
--ref_adaptor_path -b (string) default: /nCoV-2019_adapters.fasta   | Help ref_adaptor_path filepath of your fasta adaptor/ primer file (forward only, no rev. comp's) is
--ref_path -c (string) default: /ref/MN908947.3.fasta   | Help filepath of the reference sequence
--config_path -d (string) default: /ref/MN908947.3.config   | Help filepath of .config file to generate snpeff
--ram -e (int) default: 4000   | Help set ram explicitly for bbmap suite (in MB 4000 = 4GB)
--minallelefraction -f (float) default: 0.1   | Help fraction of called variant at that position (Allele Depth / Total Depth)
--minlength -g (int) default: 250   | Help minimum merged read length. Shorter reads will be removed
--maxlength -h (int) default: 450   | Help maximum merged read length.  Longer reads will be removed.
--ncbi_accession -i (string) default: MN908947.3   | Help Ncbi_accession name, must match the snp_eff name in config
--minlength_diff -j (int) default: 60   | Help bp lower limit of merged read length compared to  the mapped inner primer region length,  strict_adapters must be set to True (t) to use
--maxlength_diff -k (int) default: 24   | Help bp upper limit of merged read length compared to  the mapped inner primer region length,  strict_adapters must be set to True (t) to use
--mincov -l (int) default: 2   | Help minimum coverage at that position
--minreads -m (int) default: 10   | Help minimum read count of the allele/variant of interest
--callsub -n (boolean) default: t   | Help Set to false (f) to remove all substitutions from the variant file
--callins -o (boolean) default: t   | Help Set to false (f) to remove all insertions from the variant file
--calldel -p (boolean) default: t   | Help Set to false (f) to remove all insertions from the variant file
--strict_adaptors -q (boolean) default: f   | Help Set to false (f) to use bbduk to filter. Set to TRUE to use a custom python script that ensures merged reads are of appropriate length and the reverse paired read primer matches the forward read primer
--minscore -r (int) default: 10   | Help minimum score in the read (not average)
--minquality -s (int) default: 10   | Help minimum quality reads with lower value will be remboed
--rarity -t (float) default: 0.1   | Help set to the same value as minallele fraction
--coverage -u (boolean) default: t   | Help set to true to use coverage scores.
--maxindel -v (int) default: 100   | Help maximum insertion (N's) a read will be allowed to add and still be considered a valid mapping.
--k -w (int) default: 21   | Help bbduk setting
--ktrim -x (string ) default: l   | Help bbduk setting (r, l)
--qtrim -y (string) default: rl   | Help bbduk setting (r, l, rl)
--trimq -z (int) default: 30   | Help bbduk setting
--hdist -A (int) default: 3   | Help bbduk setting
--rcomp -B (boolean) default: f   | Help bbduk setting
--minlength_bbduk -C (int) default: 75   | Help bbduk setting
--restrictleft -D (int) default: 32   | Help bbduk setting
--minavgquality -E (int) default: 20   | Help bbduk setting
--minbasequality -F (int) default: 20   | Help bbduk setting
--use_diff -G (boolean) default: t   | Help Use a +/- of each amplicon as set by minlength_diff and maxlength_diff instead of globally declaring minlength/maxlength regaurdless of the amplicon of interet.  Only works if strict_adaptors is set
--bed_path -H (string) default: /nCoV-2019/V3/nCoV-2019.bed   | Help File path of the bed file.  Must be set to use strict_adaptors.
--outdir -I (string) default: ../out   | Help if you want to manually set the outdir set it
--reads_per_amplicon -J (int) default: 0   | Help used for downsampling normalization 0 will not down sample.  200 or 2000 is typical
--downsampling_fasta_path -K (string) default: <blank>   | Help fasta of all the amplicons (start of forward to end of reverse primer to the reference genome given.  If your amplicons map to multiple references you cannot use the short cut.  If it is not declared it can be generated from the bedfile and reference fasta.
--fast_normalize -L (boolean) default: f   | Help use a custom script to downsample about 10-100x faster than the original algorithm"; }

ref_adaptor_path=/nCoV-2019_adapters.fasta
ref_path=/ref/MN908947.3.fasta
config_path=/ref/MN908947.3.config
ram=4000
minallelefraction=0.1
minlength=250
maxlength=450
ncbi_accession=MN908947.3
minlength_diff=60
maxlength_diff=24
mincov=2
minreads=10
callsub=t
callins=t
calldel=t
strict_adaptors=f
minscore=10
minquality=10
rarity=0.1
coverage=t
maxindel=100
k=21
ktrim=l
qtrim=rl
trimq=30
hdist=3
rcomp=f
minlength_bbduk=75
restrictleft=32
minavgquality=20
minbasequality=20
use_diff=t
bed_path=/nCoV-2019/V3/nCoV-2019.bed
outdir=../out
reads_per_amplicon=0
fast_normalize=f

for arg in "$@"; do
shift
case "$arg" in
	"--sample_dir") set -- "$@" "-a" ;;
	"--ref_adaptor_path") set -- "$@" "-b" ;;
	"--ref_path") set -- "$@" "-c" ;;
	"--config_path") set -- "$@" "-d" ;;
	"--ram") set -- "$@" "-e" ;;
	"--minallelefraction") set -- "$@" "-f" ;;
	"--minlength") set -- "$@" "-g" ;;
	"--maxlength") set -- "$@" "-h" ;;
	"--ncbi_accession") set -- "$@" "-i" ;;
	"--minlength_diff") set -- "$@" "-j" ;;
	"--maxlength_diff") set -- "$@" "-k" ;;
	"--mincov") set -- "$@" "-l" ;;
	"--minreads") set -- "$@" "-m" ;;
	"--callsub") set -- "$@" "-n" ;;
	"--callins") set -- "$@" "-o" ;;
	"--calldel") set -- "$@" "-p" ;;
	"--strict_adaptors") set -- "$@" "-q" ;;
	"--minscore") set -- "$@" "-r" ;;
	"--minquality") set -- "$@" "-s" ;;
	"--rarity") set -- "$@" "-t" ;;
	"--coverage") set -- "$@" "-u" ;;
	"--maxindel") set -- "$@" "-v" ;;
	"--k") set -- "$@" "-w" ;;
	"--ktrim") set -- "$@" "-x" ;;
	"--qtrim") set -- "$@" "-y" ;;
	"--trimq") set -- "$@" "-z" ;;
	"--hdist") set -- "$@" "-A" ;;
	"--rcomp") set -- "$@" "-B" ;;
	"--minlength_bbduk") set -- "$@" "-C" ;;
	"--restrictleft") set -- "$@" "-D" ;;
	"--minavgquality") set -- "$@" "-E" ;;
	"--minbasequality") set -- "$@" "-F" ;;
	"--use_diff") set -- "$@" "-G" ;;
	"--bed_path") set -- "$@" "-H" ;;
	"--outdir") set -- "$@" "-I" ;;
	"--reads_per_amplicon") set -- "$@" "-J" ;;
	"--downsampling_fasta_path") set -- "$@" "-K" ;;
	"--fast_normalize") set -- "$@" "-L" ;;
	*) set -- "$@" "$arg"
esac
done

while getopts a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:A:B:C:D:E:F:G:H:I:J:K:L: opt ; do
case $opt in
	a) sample_dir=$OPTARG ;;
	b) ref_adaptor_path=$OPTARG ;;
	c) ref_path=$OPTARG ;;
	d) config_path=$OPTARG ;;
	e) ram=$OPTARG ;;
	f) minallelefraction=$OPTARG ;;
	g) minlength=$OPTARG ;;
	h) maxlength=$OPTARG ;;
	i) ncbi_accession=$OPTARG ;;
	j) minlength_diff=$OPTARG ;;
	k) maxlength_diff=$OPTARG ;;
	l) mincov=$OPTARG ;;
	m) minreads=$OPTARG ;;
	n) callsub=$OPTARG ;;
	o) callins=$OPTARG ;;
	p) calldel=$OPTARG ;;
	q) strict_adaptors=$OPTARG ;;
	r) minscore=$OPTARG ;;
	s) minquality=$OPTARG ;;
	t) rarity=$OPTARG ;;
	u) coverage=$OPTARG ;;
	v) maxindel=$OPTARG ;;
	w) k=$OPTARG ;;
	x) ktrim=$OPTARG ;;
	y) qtrim=$OPTARG ;;
	z) trimq=$OPTARG ;;
	A) hdist=$OPTARG ;;
	B) rcomp=$OPTARG ;;
	C) minlength_bbduk=$OPTARG ;;
	D) restrictleft=$OPTARG ;;
	E) minavgquality=$OPTARG ;;
	F) minbasequality=$OPTARG ;;
	G) use_diff=$OPTARG ;;
	H) bed_path=$OPTARG ;;
	I) outdir=$OPTARG ;;
	J) reads_per_amplicon=$OPTARG ;;
	K) downsampling_fasta_path=$OPTARG ;;
	L) fast_normalize=$OPTARG ;;
	*) usage; exit 1;;
esac
done


# untar file

cd ${sample_dir}
#pwd
#
echo ${outdir}
if [[ "${outdir:0:2}" = ".." ]]
    then
    out_dir=`pwd`
    parent_dir="$( dirname "${out_dir}" )"
    outdir="${parent_dir}${outdir:2}"
fi


echo ${outdir}
echo ${sample_dir}
find * -maxdepth 1 -type f -name "*-R1.fastq.gz" > fastq_list.txt
# loop trough match file names
for fastq_file in `cat fastq_list.txt`;  do

    sample_name=${fastq_file%????????????}
    echo ${sample_name}
    r1_fastq=${sample_name}-R1.fastq.gz
    r2_fastq=${sample_name}-R2.fastq.gz

    # make an out folder in the parent directory and dont error if it already exists
    mkdir -p ${outdir}

    ######
    # Normaliztion through Downsampling
    # Set downsampling to a positive number,
    ######
    # get the first character of the string
    strict_a=${strict_adaptors:0:1}
    # make the first character uppercase if it is T as is t, T, True, true TRUE etc. it will pass the condition
    if [[ "${strict_a^^}" = "T" ]]
    then

        echo "USING STRICT ADAPTORS"
        java -ea -Xmx${ram}m -Xms${ram}m -Djava.library.path=/bin/bbmap/jni/ -cp /bin/bbmap/current/ jgi.BBMerge \
        in=${r1_fastq} \
        in2=${r2_fastq} \
        out=${outdir}/${sample_name}_custom.not_norm.merged.fastq

        if [[ "${reads_per_amplicon}" -ne "0" ]]
        then
            # get the first character of the string
            fast_n=${fast_normalize:0:1}
            # make the first character uppercase if it is T as is t, T, True, true TRUE etc. it will pass the condition
            if [[ "${fast_n^^}" = "T" ]]
            then
                python3 /trim_strict.py \
                --bed_filepath=${bed_path} \
                --trim_primer=f \
                --minlength=${minlength} \
                --maxlength=${maxlength} \
                --use_diff=${use_diff} \
                --minlength_diff=${minlength_diff} \
                --maxlength_diff=${maxlength_diff} \
                --sample_name=${sample_name} \
                --in_merged_filepath=${outdir}/${sample_name}_custom.not_norm.merged.fastq \
                --outdir=${outdir} \
                --ref_fasta_filepath=${ref_path} \
                --primer_fasta_filepath=${ref_adaptor_path}

                java -ea -Xmx${ram}m -Xms${ram}m -cp /bin/bbmap/current/ align2.BBMap build=1 in=${outdir}/${sample_name}_custom_filtered.fastq \
                outm=${outdir}/${sample_name}.premap.bam \
                ref=${ref_path} \
                overwrite=t \
                nodisk=t
                rm -f ${outdir}/${sample_name}_custom_filtered.fastq
                samtools collate ${outdir}/${sample_name}.premap.bam ${outdir}/${sample_name}.random

                python3 /amplicon_normalize_from_bam.py \
                --bam_inpath=${outdir}/${sample_name}.random.bam \
                --bam_outpath=${outdir}/${sample_name}.norm.bam \
                --bed_path=${bed_path} \
                --reads_per_amplicon=${reads_per_amplicon}
                samtools sort -o ${outdir}/${sample_name}.norm.sorted.bam --reference ${ref_path} ${outdir}/${sample_name}.norm.bam
                samtools index ${outdir}/${sample_name}.norm.sorted.bam
                java -ea -Xms${ram}m -cp /bin/bbmap/current/ jgi.ReformatReads \
                in=${outdir}/${sample_name}.norm.sorted.bam \
                out=${outdir}/${sample_name}_custom.merged.fastq

            else
                echo "python3 /amplicon_normalization.py \
    --downsampling_fasta_path=${downsampling_fasta_path} \
    --reads_per_amplicon=${reads_per_amplicon} \
    --outdir=${outdir} \
    --min_length_diff=${minlength_diff} \
    --bed_path=${bed_path} \
    --ref_path=${ref_path} \
    --merged_fasta_path=${outdir}/${sample_name}_custom.not_norm.merged.fastq \
    --output_filename=${sample_name}_custom.merged.fastq \
    --include_primer=t"
                python3 /amplicon_normalization.py \
                --downsampling_fasta_path=${downsampling_fasta_path} \
                --reads_per_amplicon=${reads_per_amplicon} \
                --outdir=${outdir} \
                --min_length_diff=${minlength_diff} \
                --bed_path=${bed_path} \
                --ref_path=${ref_path} \
                --merged_fasta_path=${outdir}/${sample_name}_custom.not_norm.merged.fastq \
                --output_filename=${sample_name}_custom.merged.fastq \
                --include_primer=t
            fi
        else
            mv ${outdir}/${sample_name}_custom.not_norm.merged.fastq ${outdir}/${sample_name}_custom.merged.fastq
        fi
        echo "python3 /trim_strict.py \
        --bed_filepath=${bed_path} \
        --trim_primer=t \
        --minlength=${minlength} \
        --maxlength=${maxlength} \
        --use_diff=${use_diff} \
        --minlength_diff=${minlength_diff} \
        --maxlength_diff=${maxlength_diff} \
        --sample_name=${sample_name} \
        --in_merged_filepath=${outdir}/${sample_name}_custom.merged.fastq \
        --outdir=${outdir} \
        --ref_fasta_filepath=${ref_path} \
        --primer_fasta_filepath=${ref_adaptor_path}"
        python3 /trim_strict.py \
        --bed_filepath=${bed_path} \
        --trim_primer=t \
        --minlength=${minlength} \
        --maxlength=${maxlength} \
        --use_diff=${use_diff} \
        --minlength_diff=${minlength_diff} \
        --maxlength_diff=${maxlength_diff} \
        --sample_name=${sample_name} \
        --in_merged_filepath=${outdir}/${sample_name}_custom.merged.fastq \
        --outdir=${outdir} \
        --ref_fasta_filepath=${ref_path} \
        --primer_fasta_filepath=${ref_adaptor_path}

        java -ea -Xmx${ram}m -Xms${ram}m -cp /bin/bbmap/current/ jgi.BBDuk \
        in=${outdir}/${sample_name}_custom_filtered.fastq \
        out=${outdir}/${sample_name}.trimmed.merged.fastq.gz \
        qtrim=${qtrim} \
        trimq=${trimq} \
        minavgquality=${minavgquality} \
        minbasequality=${minbasequality} \
        overwrite=t
    else
        # filter for the reads that start with the primer on the list
        java -ea -Xmx${ram}m -Xms${ram}m -cp /bin/bbmap/current/ jgi.BBDuk \
        in=${r1_fastq} \
        in2=${r2_fastq} \
        outm=${outdir}/${sample_name}_R1_filtered.fastq.gz \
        outm2=${outdir}/${sample_name}_R2_filtered.fastq.gz \
        ref=${ref_adaptor_path} \
        k=${k} \
        hdist=${hdist} \
        rcomp=${rcomp} \
        minlength=${minlength_bbduk} \
        restrictleft=${restrictleft} \
        overwrite=t

        # with the filtered reads, trim the reads and low quality ends
        java -ea -Xmx${ram}m -Xms${ram}m -cp /bin/bbmap/current/ jgi.BBDuk \
        in=${outdir}/${sample_name}_R1_filtered.fastq.gz \
        in2=${outdir}/${sample_name}_R2_filtered.fastq.gz \
        out=${outdir}/${sample_name}_R1_cleaned.fastq.gz \
        out2=${outdir}/${sample_name}_R2_cleaned.fastq.gz \
        ref=${ref_adaptor_path} \
        ktrim=${ktrim} \
        k=${k} \
        qtrim=${qtrim} \
        trimq=${trimq} \
        hdist=${hdist} \
        rcomp=${rcomp} \
        minlength=${minlength_bbduk} \
        restrictleft=${restrictleft} \
        overwrite=t

        # Merge the reads to ensure a overlap of somesort
        # default minoverlap=12
        java -ea -Xmx${ram}m -Xms${ram}m -Djava.library.path=/bin/bbmap/jni/ -cp /bin/bbmap/current/ jgi.BBMerge \
        in=${outdir}/${sample_name}_R1_cleaned.fastq.gz \
        in2=${outdir}/${sample_name}_R2_cleaned.fastq.gz \
        out=${outdir}/${sample_name}.not_norm.merged.fastq.gz \
        overwrite=t

        if [[ "${reads_per_amplicon}" -ne "0" ]]
        then
        echo "python3 /amplicon_normalization.py \
--downsampling_fasta_path=${downsampling_fasta_path} \
--reads_per_amplicon=${reads_per_amplicon} \
--outdir=${outdir} \
--min_length_diff=${minlength_diff} \
--bed_path=${bed_path} \
--ref_path=${ref_path} \
--merged_fasta_path=${outdir}/${sample_name}.not_norm.merged.fastq.gz \
--output_filename=${sample_name}_custom.merged.fastq.gz \
--include_primer=f \
--ram=${ram}"

            python3 /amplicon_normalization.py \
            --downsampling_fasta_path=${downsampling_fasta_path} \
            --reads_per_amplicon=${reads_per_amplicon} \
            --outdir=${outdir} \
            --min_length_diff=${minlength_diff} \
            --bed_path=${bed_path} \
            --ref_path=${ref_path} \
            --merged_fasta_path=${outdir}/${sample_name}.not_norm.merged.fastq.gz \
            --output_filename=${sample_name}.merged.fastq.gz \
            --include_primer=f \
            --ram=${ram}
        else
            mv ${outdir}/${sample_name}.not_norm.merged.fastq.gz ${outdir}/${sample_name}.merged.fastq.gz
        fi

        # Trim reads that are too short or too long (inner amplicon range is 323 to 380)
        java -ea -Xms${ram}m -cp /bin/bbmap/current/ jgi.ReformatReads \
        in=${outdir}/${sample_name}.merged.fastq.gz \
        maxlength=${maxlength} \
        minlength=${minlength} \
        out=${outdir}/${sample_name}.trimmed.merged.fastq.gz \
        overwrite=t

    fi
    # the reference and index files that are generated in previous steps do not work  as it maps incorrectly
    # delete the genome and index folders
    rm -rf ref/genome
    rm -rf ref/index
    echo "java -ea -Xmx${ram}m -Xms${ram}m -cp /bin/bbmap/current/ align2.BBMap build=1 in=${outdir}/${sample_name}.trimmed.merged.fastq.gz \
    outm=${outdir}/${sample_name}.bam \
    ref=${ref_path} \
    maxindel=${maxindel} \
    overwrite=t"

    # map using bbmap.  Do not allow for more than 100 indels
    java -ea -Xmx${ram}m -Xms${ram}m -cp /bin/bbmap/current/ align2.BBMap build=1 in=${outdir}/${sample_name}.trimmed.merged.fastq.gz \
    outm=${outdir}/${sample_name}.bam \
    ref=${ref_path} \
    maxindel=${maxindel} \
    overwrite=t \
    nodisk=t

    # sort and index the files (using samtools)
    samtools sort -o ${outdir}/${sample_name}.sorted.bam --reference ${ref_path} ${outdir}/${sample_name}.bam
    samtools index ${outdir}/${sample_name}.sorted.bam

    # call the variants
    java -ea -Xmx${ram}m -Xms${ram}m -cp /bin/bbmap/current/ var2.CallVariants in=${outdir}/${sample_name}.sorted.bam \
    ref=${ref_path} \
    minallelefraction=${minallelefraction} \
    rarity=${rarity} \
    coverage=${coverage} \
    calldel=${calldel} \
    callins=${callins} \
    callsub=${callsub} \
    mincov=${mincov} \
    minreads=${minreads} \
    vcf=${outdir}/${sample_name}.bbmap.vcf \
    minscore=${minscore} \
    minquality=${minquality} \
    overwrite=t

    snpEff -c ${config_path} -ud -onlyProtein ${ncbi_accession} ${outdir}/${sample_name}.bbmap.vcf > ${outdir}/${sample_name}.bbmap.ann.vcf

done
