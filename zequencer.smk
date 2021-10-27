from snakemake.utils import min_version
import os

min_version("6.0")

###########################################################################
# Set up default folders for input_dir  & output_dir if they do not exist #
###########################################################################
src_dir = "./src"
if 'src_dir' not in config:
    config['src_dir'] = src_dir
default_dict_path = os.path.join(config['src_dir'], 'variant_calling_rules.json')

def read_json(filepath):
    """
    Basic json read function, throws warnings and returns [] if file does not exist or is not of valid format
    """
    import os
    import json
    if not os.path.exists(filepath):
        print("Warning: condor_q_held.json Does not exist")
        return {}
    try:
        with open(filepath) as f_in:
            return json.load(f_in)

    except ValueError:
        print("Warning: json is blank (none held) or bad format")
        return {}
def deep_merge_dicts(original, incoming):
    """
    Deep merge two dictionaries. Modifies original.
    For key conflicts if both values are:
     a. dict: Recursively call deep_merge_dicts on both values.
     c. any other type: Value is overridden.
     d. conflicting types: Value is overridden.
    """
    for key in incoming:
        if key in original:
            if isinstance(original[key], dict) and isinstance(incoming[key], dict):
                deep_merge_dicts(original[key], incoming[key])

            else:
                original[key] = incoming[key]
        else:
            original[key] = incoming[key]
    return original


####
# read the default dictionary with the hardcoded (yuck path)
default_config_dict = read_json(default_dict_path)

config = deep_merge_dicts(default_config_dict, config)

input_dir = "./"
out_dir = "out"
ref_dir = "ref"
os.makedirs(config['ref_dir'],exist_ok=True)
os.makedirs(config['out_dir'],exist_ok=True)
primer_amplicon_path = os.path.join(config['ref_dir'],"primer_amplicon.csv")
if 'primer_amplicon_path' not in config:
    config['primer_amplicon_path'] = primer_amplicon_path

ref_adaptor_path = os.path.join(config['ref_dir'],"ref_adaptor.fasta")
if 'ref_adaptor_path' not in config:
    config['ref_adaptor_path'] = ref_adaptor_path
use_amplicon_norm_bbduk_trim = True
if 'use_amplicon_norm_bbduk_trim' in config:
    if config['use_amplicon_norm_bbduk_trim'][0].upper() == 'F':
        use_amplicon_norm_bbduk_trim = False

insert_fasta_path = os.path.join(config['ref_dir'], "insert_amplicons.fasta")
if 'insert_fasta_path' not in config:
    config['insert_fasta_path'] = insert_fasta_path

src_dir = "./src"
if 'src_dir' not in config:
    config['src_dir'] = src_dir
##################################
# Get Sample list using patterns #
##################################
# needs to end in -R1.fastq.gz and reverse file as -R2.fastq.gz
sample_list = os.listdir(config['input_dir'])
print(sample_list)
SAMPLES_PRE = [x[:-12] for x in sample_list if (len(x) > 12) and (x[-12:] == '-R1.fastq.gz' and x[:2] != '._')]
SAMPLES = []
for sample_i in SAMPLES_PRE:
    if '{0}-R2.fastq.gz'.format(sample_i) in sample_list:
        SAMPLES.append(sample_i)
# SAMPLES=['ZIKV-DAKAR-41524-06192021_Rep01']
config['SAMPLES'] = SAMPLES
print(config)

#################################
# Begin the Modulized Workflow #
################################
snake_make_rules = os.path.join(config['src_dir'], 'variant_calling_rules.smk')
module variant_calling_rules:
    # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
    snakefile: snake_make_rules
    config: config

use rule all from variant_calling_rules

use rule download_ncbi_reference from variant_calling_rules

use rule create_insert_fasta_path_from_primer_amplicon from variant_calling_rules

use rule create_primer_amp_from_ref_adaptor_fasta from variant_calling_rules

use rule create_primer_fasta_ref_bed from variant_calling_rules

if use_amplicon_norm_bbduk_trim:

    use rule bbduk_filter_primer_reads from variant_calling_rules

    use rule bbduk_trim_primer_and_quality from variant_calling_rules

    use rule bbmerge_reads_from_bbduk_filter_trim from variant_calling_rules

    use rule slow_amplicon_normalization from variant_calling_rules

    use rule bbmap_reformat_trim_read_length from variant_calling_rules

    use rule bbmap_filtered_trimmed from variant_calling_rules

else:

    use rule bbmerge_reads from variant_calling_rules

    use rule bbmap_merged_reads from variant_calling_rules

    use rule filter_downsample_mask_bam from variant_calling_rules

use rule bbmap_call_variants from variant_calling_rules


