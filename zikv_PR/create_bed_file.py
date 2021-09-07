# docker run -it -v /Users:/Users dockerreg.chtc.wisc.edu/dabaker3/zequencer_ncov19:v1
# cd /Users/dabaker3/zequencer/zequencer_ncov19/zequencer_ncov19
# python3 /Users/dabaker3/zequencer/zequencer_ncov19/zequencer_ncov19/get_config_for_ref.py \
# /Users/dabaker3/zequencer/zequencer_ncov19/zequencer_ncov19 \
# accession KU501215.1 \
# dabaker3@wisc.edu
#
#
#
# bbmap.sh in=/Users/dabaker3/zequencer/zequencer_ncov19/zequencer_ncov19/zikv_PR/ZIKV_PR_Primer_Sequences.fasta \
# ref=/Users/dabaker3/zequencer/zequencer_ncov19/zequencer_ncov19/ref/KU501215.1.fa \
# outm=/Users/dabaker3/zequencer/zequencer_ncov19/zequencer_ncov19/zikv_PR/zikv_PR_primer_map.sam
#
#
# grep -v '^@' /Users/dabaker3/zequencer/zequencer_ncov19/zequencer_ncov19/zikv_PR/zikv_PR_primer_map.sam > /Users/dabaker3/zequencer/zequencer_ncov19/zequencer_ncov19/zikv_PR/zikv_PR_primer_map.tsv
import pandas as pd

input_tsv = '/Users/dabaker3/zequencer/zequencer_ncov19/zequencer_ncov19/zikv_PR/zikv_PR_primer_map.tsv'
df = pd.read_csv(input_tsv, sep='\t',  header=None,
                         names=['NAME', 'FLAG', 'REF_NAME', 'START_POSITION', 'SCORE', 'CIGAR', 'SYMBOL', 'blank2', 'blank3', 'SEQ', 'BLANK_4','BLANK_5','BLANK_6'])

df['START_POSITION'] = df['START_POSITION'] - 1
df['SEQ_LEN'] = [len(x) for x in df['SEQ']]
df['END_POSITION'] = df['START_POSITION'] + df['SEQ_LEN']
df['PLUS'] = ['+' if x[-2:] == 'FT' else '-' for x in df['NAME']]
df['NUMBER'] = [int(x.split('_')[2]) for x in df['NAME']]
df['EVEN'] = ['ZIK_400_1' if x % 2 == 1 else 'ZIK_400_2' for x in df['NUMBER']]
df = df[['REF_NAME', 'START_POSITION', 'END_POSITION', 'NAME', 'EVEN', 'PLUS']]
df.to_csv('/Users/dabaker3/zequencer/zequencer_ncov19/zequencer_ncov19/zikv_PR/zikv_PR.bed', sep='\t', header=None, index=False)