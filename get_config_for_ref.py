import re
import sys
from Bio import SeqIO
from Bio import Entrez
import subprocess
from shutil import copyfile
import os
# pipeline_dir = '~/zequencer/zequencer_ncov19/zequencer_ncov19'
# ncbi_accession = 'KX601166.2'
# email = username@wisc.edu


pipeline_dir = sys.argv[1]
ncbi_accession = sys.argv[2]
email = sys.argv[3]

print(pipeline_dir, ncbi_accession, email)
Entrez.email = email  # Always tell NCBI who you are

# Downloading...
net_handle = Entrez.efetch(db="nucleotide", id=ncbi_accession, rettype="gb", retmode="text")
out_handle = open('ref/tmp.gbk', "w")
out_handle.write(net_handle.read())
out_handle.close()
net_handle.close()
print("Saved")

# with at least some sequences (e.g., KU501215) the NCBI sequence record contains a versioning .1 or .2suffix
# this is treated inconsistently by some tools
# to eliminate this issue, remove versioning suffix from temporary genbank file
# in testing, discovered that snpEff uses chromsome name from LOCUS field of Genbank file
# other tools use ncbi_accession field for chromosome name
# use a regular expression to change the value of the LOCUS field to match value of ncbi_accession field
# biopython expects an exact number of spaces between fields, so we need to calculate the number of spaces to add after the replaced ncbi_accession number

# there is almost certainly a more elegant way to do this with biopython

with open('ref/tmp.gbk') as infile, open('ref/tmp_cleaned.gbk', 'w') as outfile:
    for line in infile:
        # remove version info
        line = re.sub(ncbi_accession + '\.[1-9]', ncbi_accession, line)

        # overwrite locus field with ncbi_accession number

        # number of spaces after ncbi_accession
        space_ct = 23 - len(ncbi_accession)
        spacer = ' ' * space_ct

        line = re.sub('LOCUS(\s*)\S*\s*', 'LOCUS' + r'\1' + ncbi_accession + spacer, line)
        outfile.write(line)

seq_id = SeqIO.read('ref/tmp.gbk', "genbank").id

# create Genbank file
SeqIO.convert('ref/tmp.gbk', "genbank", 'ref/' + ncbi_accession + '.gbk', "genbank")

# create FASTA file
SeqIO.convert('ref/tmp.gbk', "genbank", 'ref/' + ncbi_accession + '.fa', "fasta")

subprocess.call(['samtools', 'faidx', 'ref/' + ncbi_accession + '.fa'])

# subprocess.call(['novoindex', 'ref/' + ncbi_accession + '.fa.nix', 'ref/' + ncbi_accession + '.fa'])

input_file = ['ref/' + ncbi_accession + '.gbk']
output = ['ref/' + ncbi_accession + '.config',
          'ref/' + ncbi_accession + '/snpEffectPredictor.bin',
          'ref/' + ncbi_accession + '/genes.gbk']

record = SeqIO.read(input_file[0], "genbank")
seq_id = record.id
seq_description = record.description

# write snpEff config
with open(output[0], "w") as myfile:
    myfile.write('data.dir = .\n')
    myfile.write('lof.ignoreProteinCodingAfter  : 0.95\n')
    myfile.write('lof.ignoreProteinCodingBefore : 0.05\n')
    myfile.write('lof.deleteProteinCodingBases : 0.50\n')
    myfile.write(
        'codon.Standard                                                          : TTT/F, TTC/F, TTA/L, TTG/L+, TCT/S, TCC/S, TCA/S, TCG/S, TAT/Y, TAC/Y, TAA/*, TAG/*, TGT/C, TGC/C, TGA/*, TGG/W, CTT/L, CTC/L, CTA/L, CTG/L+, CCT/P, CCC/P, CCA/P, CCG/P, CAT/H, CAC/H, CAA/Q, CAG/Q, CGT/R, CGC/R, CGA/R, CGG/R, ATT/I, ATC/I, ATA/I, ATG/M+, ACT/T, ACC/T, ACA/T, ACG/T, AAT/N, AAC/N, AAA/K, AAG/K, AGT/S, AGC/S, AGA/R, AGG/R, GTT/V, GTC/V, GTA/V, GTG/V, GCT/A, GCC/A, GCA/A, GCG/A, GAT/D, GAC/D, GAA/E, GAG/E, GGT/G, GGC/G, GGA/G, GGG/G\n')
    myfile.write('\n' + seq_id + '.genome : ' + seq_description)

# create temporary Genbank file named genes.gbk as needed by snpeff
os.makedirs(os.path.join(pipeline_dir, 'ref', ncbi_accession), exist_ok=True)
copyfile(os.path.join(pipeline_dir, input_file[0]), os.path.join(pipeline_dir, output[2]))

# create snpeff database
cmd = ['snpEff', 'build', '-c', output[0], '-genbank', seq_id]

subprocess.call(cmd)

print(' '.join(cmd))
