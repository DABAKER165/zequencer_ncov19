#!/usr/bin/env python

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
    parser.add_argument('--primer_amplicon_path',
                        type=str,
                        help='filepath to the previously created custom primer_amplicon_ref.',
                        required=True)
    parser.add_argument('--in1',
                        type=str,
                        help='input bam file that has been collated and made from merged primers.',
                        required=True)
    parser.add_argument('--outdir',
                        type=str,
                        help='directory where your files are outputted to',
                        default=None,
                        required=False)
    parser.add_argument('--outf',
                        type=str,
                        help='directory where your files are outputted to',
                        required=True)
    parser.add_argument('--sample_name',
                        type=str,
                        help='prefix for the filename',
                        required=True)
    parser.add_argument('--mask_primers',
                        type=str2bool,
                        default=True,
                        help='will trim the primer in the output fasta if set to true.  If set to false, will only filter if the paired end primers are not an exact match to the paired end primer set.',
                        required=False)

    parser.add_argument('--reads_per_amplicon',
                        type=int,
                        default=800,
                        help='maximum reads per amplicon.  Overlapping regions may get higher depth of coverage.',
                        required=False)
    parser.add_argument('--downsample',
                        type=str2bool,
                        help='downsample the reads ',
                        default=True,
                        required=False)

    parser.add_argument('--max_del',
                        type=int,
                        default=1,
                        help='maximum deletions allowed to match primer',
                        required=False)
    parser.add_argument('--max_ins',
                        type=int,
                        default=1,
                        help='maximum insertions allowed to match primer',
                        required=False)
    parser.add_argument('--max_sub',
                        type=int,
                        default=1,
                        help='maximum substitutions allowed to match primer',
                        required=False)
    # parser.add_argument('--ref_fasta_filepath',
    #                     type=str,
    #                     help='filepath to the reference fasta.  Either the (ref_fasta_filepath AND bed_filepath) must be set or the primer_fasta_filepath must be set.',
    #                     required=False)
    # parser.add_argument('--primer_fasta_filepath',
    #                     type=str,
    #                     help='filepath to the primer fasta for all of the primer set.  Either the (ref_fasta_filepath AND bed_filepath) must be set or the primer_fasta_filepath must be set. ',
    #                     required=False)


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

in1 = args.in1
outf = args.outf
primer_amplicon_path = args.primer_amplicon_path
sample_name = args.sample_name
downsample = args.downsample
reads_per_amplicon = args.reads_per_amplicon
max_del = args.max_del
max_ins = args.max_ins
max_sub = args.max_sub
mask_primers = args.mask_primers
outdir = args.outdir

import pysam
import pandas as pd
from os import path
from sys import stdout
if outdir is None:
    outdir = path.dirname(outf)


def compare_primer(primer_len=0,
                   is_left=True,
                   diff_len=0,
                   cigartuples=[],
                   max_del=0,
                   max_ins=0,
                   max_sub=0):
    consumes_reference = [0, 2, 3, 7, 8]
    # consumesQuery = [0,1,4,7,8]
    position = 0
    ins_count = 0
    del_count = 0
    sub_count = 0
    result_forward = (False, False, 0, 0)

    if diff_len < 0:
        del_count += abs(diff_len)
        if del_count > max_del:
            return (False, False, 0, 0)
    else:
        ins_count += abs(diff_len)
        if ins_count > max_ins:
            return (False, False, 0, 0)
    for cigar_i in cigartuples:
        if cigar_i[0] in consumes_reference:
            position += cigar_i[1]
            if cigar_i[0] == 2:
                del_count += cigar_i[1]
            if cigar_i[0] == 8:
                sub_count += cigar_i[1]
        else:
            if cigar_i[0] == 1:
                ins_count += cigar_i[1]
        if (ins_count > max_ins) or (del_count > max_del) or (sub_count > max_sub):
            return (False, False, ins_count, del_count)

        if position >= primer_len:
            return (True, (ins_count + del_count + sub_count) == 0, ins_count, del_count)

    return (False, False, 0, 0)


def mask_cigar(segment_pos, segment_end, cigartuples, primer_start, primer_end, left_primer):
    cigartuples_temp = copy(cigartuples)
    consume_1 = [0, 7, 8, 2, 3]
    consume_2 = [2, 3]
    consume_3 = [1, 4]
    consumesNeither = [5, 6]
    clip_length = primer_end - primer_start
    clip_length_i = clip_length
    cigar_length = 0
    extra = 0
    extra_ins = 0
    extra_del = 0
    if left_primer:

        extra = primer_start - segment_pos
    # print('pre extra: {0}'.format(extra))
    else:
        extra = segment_end - primer_end
    clip_length = primer_end - primer_start + extra
    clip_length_i = clip_length

    for cigar_i in cigartuples:
        if cigar_i[0] in consume_1:
            cigartuples_temp.pop(0)
            cigar_length = cigar_i[1] - clip_length_i
            if cigar_i[0] in consume_2:
                extra_del -= cigar_i[1]
            if cigar_i[1] < clip_length_i:
                clip_length_i = clip_length_i - cigar_i[1]
            else:
                if cigar_length > 0:
                    cigartuples_temp.insert(0, (cigar_i[0], cigar_length))  # - extra_del

                # print(clip_length, extra , extra_ins, extra_del)
                cigartuples_temp.insert(0, (4, clip_length + extra_ins + extra_del))
                extra = extra_ins + extra_del
                break

        elif cigar_i[0] in consume_3:
            cigartuples_temp.pop(0)
            extra_ins += cigar_i[1]

        elif cigar_i[0] in consumesNeither:
            cigartuples_temp.pop(0)

    return cigartuples_temp, extra


from copy import copy

perfect_both_match_count = 0
perfect_single_match_count = 0
imperfect_match_count = 0
single_bad_primer_count = 0
both_bad_primer_count = 0
match_count = 0
mismatch_count = 0

# Taken from the artic but they were missing the s
# consumesReference lookup for if a CIGAR operation consumes the reference sequence
# M alignment match (can be a sequence match or mismatch) yes yes
# I 1 insertion to the reference yes no
# D 2 deletion from the reference no yes
# N 3 skipped region from the reference no yes
# S 4 soft clipping (clipped sequences present in SEQ) yes no
# H 5 hard clipping (clipped sequences NOT present in SEQ) no no
# P 6 padding (silent deletion from padded reference) no no
# = 7 sequence match yes yes
# X 8 sequence mismatch yes yes


# you must collate the file FIRST if you want ot downsample asto shuffle the reads

df_primer_amp = pd.read_csv(primer_amplicon_path, sep=',')
amplicon_count_dict = {}
amplicon_rev_dict = {}
amplicon_fwd_dict = {}
amplicon_list = list(df_primer_amp['AMPLICON'].unique())
for item in amplicon_list:
    amplicon_count_dict[item] = 0
    amplicon_rev_dict[item] = 0
    amplicon_fwd_dict[item] = 0
left_start_list = df_primer_amp['START_LEFT']
left_end_list = df_primer_amp['END_LEFT']
right_start_list = df_primer_amp['START_RIGHT']
right_end_list = df_primer_amp['END_RIGHT']
primer_group_list = list(df_primer_amp['PRIMER_GROUP'])
primer_left_length_list = [len(x) for x in df_primer_amp['PRIMER_SEQUENCE_LEFT']]
primer_right_length_list = [len(x) for x in df_primer_amp['PRIMER_SEQUENCE_RIGHT']]
amplicon_all_list = list(df_primer_amp['AMPLICON'])
i = 0
m = 0
infile = pysam.AlignmentFile(in1, "rb")
bam_header = infile.header.copy().to_dict()
outfile = pysam.AlignmentFile(outf, "wb", header=bam_header)
primer_group_unique_list = []
outfile_groups = {}
if 'PRIMER_GROUP' in df_primer_amp.columns:
    primer_group_unique_list = list(df_primer_amp['PRIMER_GROUP'].unique())
for pg_i in primer_group_unique_list:
    outfile_groups[str(pg_i)] = pysam.AlignmentFile('{0}{1}.bam'.format(outf[:-4], pg_i), "wh", header=bam_header)
for segment in infile:
    if i % 10000 == 0:
        stdout.write('\r')

        stdout.write('{0} segments finished.'.format(i))
        stdout.flush()
    i += 1
    compare_start = segment.reference_start
    compare_end = segment.reference_end
    cigartuples = segment.cigartuples
    closest_tuple = min(enumerate(left_start_list), key=lambda x: abs(x[1] - compare_start))
    if downsample and (amplicon_count_dict[amplicon_all_list[closest_tuple[0]]] >= reads_per_amplicon):
        continue
    diff_left = closest_tuple[1] - compare_start

    diff_right = compare_end - right_end_list[closest_tuple[0]]

    result_left = compare_primer(primer_len=primer_left_length_list[closest_tuple[0]],
                                 is_left=True,
                                 diff_len=diff_left,
                                 cigartuples=cigartuples,
                                 max_del=max_del,
                                 max_ins=max_ins,
                                 max_sub=max_sub)

    cigartuples.reverse()
    result_right = compare_primer(primer_len=primer_right_length_list[closest_tuple[0]],
                                  is_left=False,
                                  diff_len=diff_right,
                                  cigartuples=cigartuples,
                                  max_del=max_del,
                                  max_ins=max_ins,
                                  max_sub=max_sub)

    if result_left[0] and result_right[0]:
        m += 1
        match_count += 1
        amplicon_count_dict[amplicon_all_list[closest_tuple[0]]] += 1
        if segment.is_reverse:
            amplicon_rev_dict[amplicon_all_list[closest_tuple[0]]] += 1
        else:
            amplicon_fwd_dict[amplicon_all_list[closest_tuple[0]]] += 1

        cigartuples.reverse()

        ref_start = segment.reference_start
        cigartuples.reverse()

        ref_ending = segment.reference_end

        cigartuples, extra = mask_cigar(segment_pos=segment.pos,
                                        segment_end=segment.reference_end,
                                        cigartuples=cigartuples,
                                        primer_start=right_start_list[closest_tuple[0]],
                                        primer_end=right_end_list[closest_tuple[0]],
                                        left_primer=False)

        cigartuples.reverse()

        segment_end1 = segment.reference_end
        cigartuples, extra = mask_cigar(segment_pos=segment.reference_start,
                                        cigartuples=cigartuples,
                                        segment_end=segment.reference_end,
                                        primer_start=left_start_list[closest_tuple[0]],
                                        primer_end=left_end_list[closest_tuple[0]],
                                        left_primer=True)

        segment.cigartuples = cigartuples
        segment.pos = segment.pos + cigartuples[0][1] - extra

        outfile_groups[str(primer_group_list[closest_tuple[0]])].write(segment)

        outfile.write(segment)
        if result_left[1] and result_right[1]:
            perfect_both_match_count += 1
        elif not result_left[1] and not result_right[1]:
            imperfect_match_count += 1
        else:
            perfect_single_match_count += 1

    else:
        mismatch_count += 1
        if not result_left[0] and not result_right[0]:
            both_bad_primer_count += 1
        else:
            single_bad_primer_count += 1

infile.close()
outfile.close()
for pg_i, outfile_g in outfile_groups.items():
    outfile_g.close()
total_count = match_count + mismatch_count
df_summary = pd.DataFrame({'MATCH_PARAMETER': ['total_count',
                                               'match_count',
                                               'mismatch_count',
                                               'perfect_both_match_count',
                                               'perfect_single_match_count',
                                               'imperfect_match_count',
                                               'single_bad_primer_count',
                                               'both_bad_primer_count'],
                           'MATCH_COUNT': [total_count,
                                           match_count,
                                           mismatch_count,
                                           perfect_both_match_count,
                                           perfect_single_match_count,
                                           imperfect_match_count,
                                           single_bad_primer_count,
                                           both_bad_primer_count],
                           'MATCH_PERCENT': [round(total_count / total_count, 4),
                                             round(match_count / total_count, 4),
                                             round(mismatch_count / total_count, 4),
                                             round(perfect_both_match_count / total_count, 4),
                                             round(perfect_single_match_count / total_count, 4),
                                             round(imperfect_match_count / total_count, 4),
                                             round(single_bad_primer_count / total_count, 4),
                                             round(both_bad_primer_count / total_count, 4)]})
df_summary.to_csv(path.join(outdir, '{0}_primer_filter_summary.csv'.format(sample_name)), index=False)
df_amplicon_count = pd.DataFrame(amplicon_count_dict.items(), columns=['AMPLICON', 'READ_COUNT'])
df_amplicon_count['DOWNSAMPLED_COUNT'] = [x if x < reads_per_amplicon else reads_per_amplicon for x in
                                          df_amplicon_count['READ_COUNT']]
df_amplicon_rev = pd.DataFrame(amplicon_rev_dict.items(), columns=['AMPLICON', 'REV_COUNT'])
df_amplicon_fwd = pd.DataFrame(amplicon_fwd_dict.items(), columns=['AMPLICON', 'FWD_COUNT'])
df_amplicon_count = df_amplicon_count.merge(df_amplicon_rev, on='AMPLICON', how='left')
df_amplicon_count = df_amplicon_count.merge(df_amplicon_fwd, on='AMPLICON', how='left')
df_amplicon_count.to_csv(path.join(outdir, '{0}_amplicon_depth_summary.csv'.format(sample_name)), index=False)
print(df_summary)
print(df_amplicon_count)
