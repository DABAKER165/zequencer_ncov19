import sys
import os
import shutil
import pyfasta
import subprocess
import gzip
import pandas as pd
from Bio import SeqIO

from argparse import ArgumentParser, ArgumentTypeError


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Boolean value expected.')


def parse_process_arrays_args(parser: ArgumentParser):
    """Parses the python script arguments from bash and makes sure files/inputs are valid"""
    parser.add_argument('--downsampling_fasta_path',
                        type=str,
                        help='filepath to fasta of all the amplicons (start of forward to end of reverse primer to the reference genome given.  If your amplicons map to multiple references you cannot use the short cut.  If it is not declared it can be generated from the bedfile and reference fasta.',
                        required=False)
    parser.add_argument('--reads_per_amplicon',
                        type=int,
                        default=2000,
                        help='up to howmany reads are randomly selected per fasta file',
                        required=True)
    parser.add_argument('--outdir',
                        type=str,
                        help='filepath where the output ends up.',
                        required=True)
    parser.add_argument('--min_length_diff',
                        type=int,
                        default=60,
                        help='number of base pairs less than the reference area that is acceptable. (typically not changed)',
                        required=False)
    parser.add_argument('--bed_path',
                        type=str,
                        help='File path of the bed file.  Must be set in conjunction with ref_path if the downsampling_fasta_path is not declared (or does not exist)',
                        required=False)
    parser.add_argument('--ref_path',
                        type=str,
                        help='filepath of the reference sequence. Must be set in conjunction with ref_path if the downsampling_fasta_path is not declared (or does not exist)',
                        required=False)
    parser.add_argument('--merged_fasta_path',
                        type=str,
                        help='input merged fasta file that is mapped',
                        required=True)
    parser.add_argument('--output_filename',
                        type=str,
                        default='downsampled.fastq',
                        help='filename of downsampled fastq output',
                        required=False)
    parser.add_argument('--include_primer',
                        type=str2bool,
                        help='if the mapped reads of the downsampling files include the primer.',
                        required=True)
    parser.add_argument('--ram',
                        type=int,
                        help='java ram limits.',
                        default=4000,
                        required=False)


def get_process_arrays_args():
    """	Inputs arguments from bash
    Gets the arguments, checks requirements, returns a dictionary of arguments
    Return: args - Arguments as a dictionary
    """
    parser = ArgumentParser()
    parse_process_arrays_args(parser)
    args = parser.parse_args()
    return args


args = get_process_arrays_args()
# same arguments to a local variable by same name as the argument
print(sys.argv)
print(args)
downsampling_fasta_path = args.downsampling_fasta_path
reads_per_amplicon = args.reads_per_amplicon
outdir = args.outdir
min_length_diff = args.min_length_diff
bed_path = args.bed_path
ref_path = args.ref_path
merged_fasta_path = args.merged_fasta_path
output_filename = args.output_filename
include_primer = args.include_primer
ram = args.ram
tmp_dir = outdir
if not (output_filename.endswith('.fq.gz') or output_filename.endswith('.fastq.gz')):
    if output_filename.endswith('.fq') or output_filename.endswith('.fastq'):
        output_filename = '{0}.gz'.format(output_filename)
    else:
        '{0}.fastq.gz'.format(output_filename)
output_filename_unzip = output_filename[:-3]

print('output_filename: {0}'.format(output_filename))
print(downsampling_fasta_path)
if (downsampling_fasta_path is None) or (downsampling_fasta_path == '') or (downsampling_fasta_path == "NONE"):
    # zero is the first position
    print("creating Downsampling fasta")
    # create a file path for the downsampling fasta
    downsampling_fasta_path = os.path.join(outdir, 'downsampling.fasta')
    # open the reference fasta file and extract the first sequence.
    fasta_sequences = SeqIO.parse(open(ref_path), 'fasta')
    # only one sequence so
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        break
    # open the bed file and parse the amplicon id's based on the nomenclature
    df_bed = pd.read_csv(bed_path, sep='\t', header=None,
                         names=['SEQ_ID', 'START', 'END', 'AMPLICON_ID', 'AMPLICON_GROUP', 'SYMBOL'])
    df_bed['END'] = df_bed['END']
    df_bed['SIDE'] = [x.split('_')[2] for x in df_bed['AMPLICON_ID']]
    df_bed['AMPLICON'] = ['_'.join(x.split('_')[:2]) for x in df_bed['AMPLICON_ID']]

    group_by_list = ['AMPLICON', 'SIDE']
    # because there are alternate primers, must use the max and min to encompass both sides of the primer
    df_bed = df_bed.groupby(group_by_list).agg({'START': 'min',
                                                'END': 'max'}).reset_index()
    # pivot so the start and ends of the forward and reverse primers are in the same row,
    df_bed_pivot = df_bed.pivot_table(index='AMPLICON', columns='SIDE', values=['START', 'END']).reset_index()
    # rename columns ass needed
    df_bed_pivot.columns = df_bed_pivot.columns.to_series().str.join('_')
    df_bed_pivot.rename(columns={'AMPLICON_': 'AMPLICON'}, inplace=True)
    # extract the sequences from the reference sequence using the start of left and end of right sequences positions
    if include_primer:
        df_bed_pivot['SEQUENCE'] = [sequence[x:y] for x, y in
                                    zip(df_bed_pivot['START_LEFT'], df_bed_pivot['END_RIGHT'])]
    else:
        df_bed_pivot['SEQUENCE'] = [sequence[x:y] for x, y in
                                    zip(df_bed_pivot['END_LEFT'], df_bed_pivot['START_RIGHT'])]
    # out put the fasta file using the amplicon column as the header.
    with open(downsampling_fasta_path, 'w') as out_file:
        for index, row in df_bed_pivot.iterrows():
            out_file.write('>{0}\n'.format(row['AMPLICON']))
            out_file.write('{0}\n'.format(row['SEQUENCE']))

# open the fasta file as using pyfasta.

f = pyfasta.Fasta(downsampling_fasta_path)
SPLIT_REF_DIR = os.path.join(tmp_dir,'split_reference_fasta')
mapped_reads_dir = os.path.join(tmp_dir,'mapped_reads')
filtered_reads_dir = os.path.join(tmp_dir, 'filtered_reads')

os.makedirs(SPLIT_REF_DIR, exist_ok=True)
os.makedirs(mapped_reads_dir, exist_ok=True)
os.makedirs(filtered_reads_dir, exist_ok=True)
# iterate over each entry in REF_FASTA and create new FASTA file with single sequence
for i in f:
    header = str(i)
    sequence = str(f[i])

    # create FASTA file
    file = open(os.path.join(SPLIT_REF_DIR, header + '.fasta'), "w")
    file.write('>' + header + '\n')
    file.write(sequence + '\n')
    file.close()

    # map meregd reads to each individual amplicon FASTA

# path to folder of individual reference files


# run BBMAP on each reference sequence individually
# save mapped reads for each file separately

# initialize counter
current_ct = 0
amplicon_ct = len(os.listdir(SPLIT_REF_DIR))
# remove output file if it already exists
if os.path.exists(os.path.join(outdir, output_filename_unzip)):
    os.remove(os.path.join(outdir, output_filename_unzip))
if os.path.exists(os.path.join(outdir, output_filename)):
    os.remove(os.path.join(outdir, output_filename))
for fn in os.listdir(SPLIT_REF_DIR):
    # print update on data processing
    current_ct += 1
    print('--Read normalizing is ' + str("%.1f" % (current_ct / amplicon_ct * 90)) + '% complete.--')

    # run bbmap
    # default sensitivty of minid=0.76 sufficient to capture only merged reads that map to a single
    # store reference in memory with nodisk
    fn_path = os.path.join(SPLIT_REF_DIR, fn)
    subprocess.call([' '.join(['java -ea -Xmx{0}m -Xms{0}m -cp /bin/bbmap/current/ align2.BBMap build=1'.format(ram),
                               'in=' + merged_fasta_path,
                               'nodisk=t',
                               'ref=' + fn_path,
                               'outm={0}'.format(os.path.join(mapped_reads_dir, '{0}.bam'.format(fn)))])],
                    shell=True)

    # run bbmap reformat.sh to filter bam file to include only reads spanning entire amplicon
    # Shelby's amplicon file includes 44bp of primer that has been trimmed during merging
    # so set "full length" to size of amplicon minus 60bp to allow some flexibility

    # get length of reference sequence
    g = pyfasta.Fasta(fn_path)

    for i in g:
        REF_LENGTH = len(str(g[i]))

    MIN_LENGTH = str(REF_LENGTH - min_length_diff)

    subprocess.call([' '.join(['java -ea -Xms{0}m -cp /bin/bbmap/current/ jgi.ReformatReads'.format(ram),
                               'in={0}'.format(os.path.join(mapped_reads_dir, '{0}.bam'.format(fn))),
                               'out={0}'.format(os.path.join(filtered_reads_dir, '{0}.filtered.fastq.gz'.format(fn))),
                               'minlength=' + MIN_LENGTH])], shell=True)

    # downsample filtered reads so each amplicon has same number of reads
    # use seqtk to filter
    # save stdout from seqtk sample to file

    with open(os.path.join(outdir, output_filename_unzip), 'a') as h:
        subprocess.call([' '.join(['seqtk',
                                   'sample',
                                   os.path.join(filtered_reads_dir, '{0}.filtered.fastq.gz'.format(fn)),
                                   str(reads_per_amplicon)])], stdout=h, shell=True)

# gzip compress downsampled FASTQ
subprocess.call(['gzip {0}'.format(os.path.join(outdir, output_filename_unzip))], shell=True)
# f_in = open(os.path.join(outdir, output_filename_unzip), 'rb')
# f_out = gzip.open(output_filename, 'wb')
# f_out.writelines(f_in)
# f_out.close()
# f_in.close()
print('--Read normalizing is 100% complete.--')
