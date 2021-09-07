#!/usr/bin/env python

from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
import os

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
    parser.add_argument('--bed_filepath',
                        type=str,
                        help='filepath to the bed file. Must be set if use_diff is set to true or it will error.',
                        required=False)
    parser.add_argument('--trim_primer',
                        type=str2bool,
                        default=True,
                        help='will trim the primer in the output fasta if set to true.  If set to false, will only filter if the paired end primers are not an exact match to the paired end primer set.',
                        required=False)
    parser.add_argument('--minlength',
                        type=int,
                        default=250,
                        help='minimum merged read length. Shorter reads will be removed',
                        required=False)
    parser.add_argument('--maxlength',
                        type=int,
                        default=450,
                        help='maximum merged read length.  Longer reads will be removed.',
                        required=False)
    parser.add_argument('--use_diff',
                        type=str2bool,
                        help='use the diff instead of the max and min lengths',
                        default=True,
                        required=False)
    parser.add_argument('--minlength_diff',
                        type=int,
                        default=24,
                        help='bp lower limit of merged read length compared to  the mapped inner primer region length,  strict_adapters must be set to True (t) to use',
                        required=False)
    parser.add_argument('--maxlength_diff',
                        type=int,
                        default=24,
                        help='bp upper limit of merged read length compared to  the mapped inner primer region length,  strict_adapters must be set to True (t) to use',
                        required=False)
    parser.add_argument('--sample_name',
                        type=str,
                        help='prefix for the filename',
                        required=True)
    parser.add_argument('--in_merged_filepath',
                        type=str,
                        help='merged input fasta file.  Use bbmerge to merge the forward and reverse files.',
                        required=True)
    parser.add_argument('--outdir',
                        type=str,
                        help='directory where your files are outputted to',
                        required=True)
    parser.add_argument('--ref_fasta_filepath',
                        type=str,
                        help='filepath to the reference fasta.  Either the (ref_fasta_filepath AND bed_filepath) must be set or the primer_fasta_filepath must be set.',
                        required=False)
    parser.add_argument('--primer_fasta_filepath',
                        type=str,
                        help='filepath to the primer fasta for all of the primer set.  Either the (ref_fasta_filepath AND bed_filepath) must be set or the primer_fasta_filepath must be set. ',
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
bed_filepath = args.bed_filepath
trim_primer = args.trim_primer
minlength = args.minlength
maxlength = args.maxlength
use_diff = args.use_diff
minlength_diff = args.minlength_diff
maxlength_diff = args.maxlength_diff
sample_name = args.sample_name
in_merged_filepath = args.in_merged_filepath
outdir = args.outdir
ref_fasta_filepath = args.ref_fasta_filepath
primer_fasta_filepath = args.primer_fasta_filepath


def test_primer_pairs(record, amplicon_pairings,
                      trim_primer=True,
                      maxlength=420,
                      minlength=280,
                      use_diff=False,
                      maxlength_diff=24,
                      minlength_diff=24):
    seq_len = len(record['sequence'])
    if (not use_diff) and ((seq_len > maxlength) or (seq_len < minlength)):
        return [False, '', '']
    for amplicon_pair in amplicon_pairings:
        for amplicon_left in amplicon_pair[0]:
            if record['sequence'].startswith(amplicon_left):
                for amplicon_right in amplicon_pair[1]:
                    if (record['sequence'].endswith(amplicon_right)):
                        if use_diff:
                            max_amp = df_amp_len_dict[amplicon_pair[2]] + maxlength_diff
                            min_amp = df_amp_len_dict[amplicon_pair[2]] - minlength_diff
                            if (seq_len > max_amp) or (seq_len < min_amp):
                                return [False, '', '']
                        if trim_primer:
                            return [True,
                                    record['sequence'][len(amplicon_left): - len(amplicon_right)],
                                    record['quality'][len(amplicon_left):- len(amplicon_right)]]
                        return [True,
                                record['sequence'],
                                record['quality']]

    return [False, '', '']


def process(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: str(v) for k, v in zip(ks, lines)}


df_amp_len_dict = {}

# get the sequence

if primer_fasta_filepath is None and ((ref_fasta_filepath is None) or (bed_filepath is None)):
    exit(1)

if (ref_fasta_filepath is not None) and (bed_filepath is not None):
    print('Using Bed file and Ref File')
    # zero is the first position
    df_bed = pd.read_csv(bed_filepath, sep='\t', header=None,
                         names=['SEQ_ID', 'START', 'END', 'AMPLICON_ID', 'AMPLICON_GROUP', 'SYMBOL'])


    fasta_sequences = SeqIO.parse(open(ref_fasta_filepath), 'fasta')
    # only one sequence so
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        break
    # MN908947.3	30	54	nCoV-2019_1_LEFT	nCoV-2019_1	+
    # zero is the first position
    df_bed = pd.read_csv(bed_filepath, sep='\t', header=None,
                         names=['SEQ_ID', 'START', 'END', 'AMPLICON_ID', 'AMPLICON_GROUP', 'SYMBOL'])

    df_bed['SIDE'] = [x.split('_')[2] for x in df_bed['AMPLICON_ID']]
    df_bed['AMPLICON'] = ['_'.join(x.split('_')[:2]) for x in df_bed['AMPLICON_ID']]

    if (primer_fasta_filepath is not None):
        print('Using primer_fasta')
        fasta_sequences = SeqIO.parse(open(primer_fasta_filepath), 'fasta')
        # only one sequence so
        df_primer = pd.DataFrame({})
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            new_row = {'AMPLICON_ID': name, 'PRIMER_SEQUENCE': sequence, 'AMPLICON': '_'.join(name.split('_')[:2]),
                       'SIDE': name.split('_')[2]}
            df_primer = df_primer.append(new_row, ignore_index=True)
        df_bed = df_bed.merge(df_primer, on=['AMPLICON_ID', 'AMPLICON', 'SIDE'], how='inner')
    else:

        primer_fasta_filepath = os.path.join(outdir, 'primers.fasta')
        df_bed['PRIMER_SEQUENCE'] = [sequence[int(x):int(y)] for x, y in zip(df_bed['START'], df_bed['END'])]
        df_bed['PRIMER_SEQUENCE'] = [x if y.upper() == 'LEFT' else str(Seq(x).reverse_complement()) for x, y in zip(df_bed['PRIMER_SEQUENCE'], df_bed['SIDE'])]

    df_adapt = df_bed
    group_by_list = ['AMPLICON', 'SIDE']
    df_bed_grouped = df_bed.groupby(group_by_list).agg({'START': 'min',
                                                        'END': 'max'}).reset_index()

    df_bed_pivot = df_bed_grouped.pivot_table(index='AMPLICON', columns='SIDE', values=['START', 'END']).reset_index()

    df_bed_pivot.columns = df_bed_pivot.columns.to_series().str.join('_')
    df_bed_pivot.rename(columns={'AMPLICON_': 'AMPLICON'}, inplace=True)

    for x, y, z in zip(df_bed_pivot['AMPLICON'], df_bed_pivot['END_RIGHT'], df_bed_pivot['START_LEFT']):
        df_amp_len_dict[x] = y - z
    print(df_amp_len_dict)
else:
    print('Using primer_fasta')
    fasta_sequences = SeqIO.parse(open(primer_fasta_filepath), 'fasta')
    # only one sequence so
    df_adapt = pd.DataFrame({})
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        new_row = {'AMPLICON_ID': name, 'PRIMER_SEQUENCE': sequence, 'AMPLICON': '_'.join(name.split('_')[:2]),
                   'SIDE': name.split('_')[2]}
        df_adapt = df_adapt.append(new_row, ignore_index=True)

if len(df_amp_len_dict) < 1:
    use_diff = False

amplicon_pairings = []
print(df_adapt)
side_list = [['LEFT', 'RIGHT'], ['RIGHT', 'LEFT']]
for amplicon_i in list(set(df_adapt['AMPLICON'])):
    df_adapt_i = df_adapt[df_adapt['AMPLICON'] == amplicon_i]
    for sides in side_list:
        df_adapt_side_i = df_adapt_i[df_adapt_i['SIDE'] == sides[0]]
        left_list = list(df_adapt_side_i['PRIMER_SEQUENCE'])
        df_adapt_side_i = df_adapt_i[df_adapt_i['SIDE'] == sides[1]]
        right_list = list(df_adapt_side_i['PRIMER_SEQUENCE'])
        right_list = [str(Seq(x).reverse_complement()) for x in right_list]
        amplicon_pairings.append([left_list, right_list, amplicon_i])

i = 1

with open(os.path.join(outdir, '{0}_custom_filtered.fastq'.format(sample_name)), 'w') as out_file:
    with open(in_merged_filepath, 'r') as fh:
        lines = []
        for line in fh:
            if i % 100000 == 0:
                print(i)
            i = i + 1
            lines.append(line.rstrip())
            if len(lines) == 4:
                record = process(lines)
                result = test_primer_pairs(record=record,
                                           amplicon_pairings=amplicon_pairings,
                                           trim_primer=trim_primer,
                                           use_diff=use_diff,
                                           maxlength=maxlength,
                                           minlength=minlength,
                                           maxlength_diff=maxlength_diff,
                                           minlength_diff=minlength_diff)
                if result[0]:
                    out_file.write('{0}\n'.format(record['name']))
                    out_file.write('{0}\n'.format(result[1]))
                    out_file.write('{0}\n'.format(record['optional']))
                    out_file.write('{0}\n'.format(result[2]))
                lines = []
