# ncbi_accession = 'KX601166.2'
from os import path


def create_cmd_req_list_optional_dict(main_cmd='',
                                      required_arg_list=[],
                                      optional_arg_dict={},
                                      delimiter='=') -> object:
    required_args = ' '.join(required_arg_list)
    optional_args = ''
    for key_i, value_i in optional_arg_dict.items():
        optional_args = '{0} {1}{2}{3}'.format(optional_args,key_i,delimiter,value_i)
    return ' '.join([main_cmd, required_args, optional_args])


rule all:
    input:
        # vcf=expand(path.join(config['out_dir'],"{sample}.vcf"),sample=config['SAMPLES']),
        ann_vcf=expand(path.join(config['out_dir'],"{sample}_ann.vcf"),sample=config['SAMPLES'])
    run:
        print('Finished!')

rule download_ncbi_reference:
    message:
        """download reference sequence from NCBI in Genbank format"""
    output:
        # temp('ref/tmp.gbk'),
        # temp('ref/tmp_cleaned.gbk')
        path.join(config['ref_dir'], config['ncbi_accession'], '.gbk'),
        path.join(config['ref_dir'], config['ncbi_accession'], '.fa')

    params:
        email=config['email'],
        ncbi_accession=config['ncbi_accession']
    run:
        from Bio import Entrez
        import re

        '''retrieve Genbank file from NCBI and copy to temporary file'''
        # print(pipeline_dir, ncbi_accession, email)

        Entrez.email = params.email  # Always tell NCBI who you are

        # Downloading...
        net_handle = Entrez.efetch(db="nucleotide",id=params.ncbi_accession,rettype="gb",retmode="text")
        out_handle = open('ref/tmp.gbk',"w")
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

        with open('ref/tmp.gbk') as infile, open('ref/tmp_cleaned.gbk','w') as outfile:
            for line in infile:
                # remove version info
                line = re.sub(params.ncbi_accession + '\.[1-9]',params.ncbi_accession,line)

                # overwrite locus field with ncbi_accession number

                # number of spaces after ncbi_accession
                space_ct = 23 - len(params.ncbi_accession)
                spacer = ' ' * space_ct

                line = re.sub('LOCUS(\s*)\S*\s*','LOCUS' + r'\1' + params.ncbi_accession + spacer,line)
                outfile.write(line)
        from Bio import SeqIO

        # seq_id = SeqIO.read('ref/tmp.gbk',"genbank").id
        # create Genbank file
        SeqIO.convert('ref/tmp.gbk',"genbank",path.join('ref', config['ncbi_accession'], '.gbk'),"genbank")
        # create FASTA file
        SeqIO.convert('ref/tmp.gbk',"genbank",path.join('ref', config['ncbi_accession'], '.fa'),"fasta")

rule bbmerge_reads:
    message:
        '''
        Merge files based on parameters.
        '''
    input:
        in1=path.join(config['input_dir'],"{sample}-R1.fastq.gz"),
        in2=path.join(config['input_dir'],"{sample}-R2.fastq.gz")
    output:
        out=path.join(config['out_dir'],"{sample}-merged.fastq.gz")
    params:
        optional_arg_dict={},
        ram=config['ram']
    threads:
        2
    run:
        print(config['SAMPLES'])
        print(['in={0}'.format(input.in1),
               'in2={0}'.format(input.in2),
               'out={0}'.format(output.out)])
        import subprocess

        main_cmd = 'java -ea -Xmx{0}m -Xms{0}m -Djava.library.path=/bin/bbmap/jni/ \
        -cp /bin/bbmap/current/ jgi.BBMerge'.format(params.ram)
        required_arg_list = ['in={0}'.format(input.in1),
                             'in2={0}'.format(input.in2),
                             'out={0}'.format(output.out)]
        run_cmd = create_cmd_req_list_optional_dict(main_cmd=main_cmd,
            required_arg_list=required_arg_list,
            optional_arg_dict=params.optional_arg_dict,
            delimiter='=')
        subprocess.call(run_cmd,shell=True)


rule bbmap_merged_reads:
    message:
        '''
        BBMAP the merged reads
        '''
    input:
        in1=path.join(config['out_dir'],"{sample}-merged.fastq.gz"),
        ref_path=config['ref_fasta_path']
    output:
        outm=path.join(config['out_dir'],"{sample}.bam"),
        outss=path.join(config['out_dir'],"{sample}-collated.bam")
    params:
        optional_arg_dict=config['bbmap_merged_reads'],
        ram=config['ram']
    threads:
        2
    run:
        import subprocess

        main_cmd = 'java -ea -Xmx{0}m -Xms{0}m -cp /bin/bbmap/current/ align2.BBMap build=1'.format(params.ram)
        required_arg_list = ['in={0}'.format(input.in1),
                             'ref={0}'.format(input.ref_path),
                             'outm={0}'.format(output.outm)]
        run_cmd = create_cmd_req_list_optional_dict(main_cmd=main_cmd,
            required_arg_list=required_arg_list,
            optional_arg_dict=params.optional_arg_dict,
            delimiter='=')
        subprocess.call(run_cmd,shell=True)
        # maxindel=${maxindel} \
        # overwrite=t \
        # nodisk=t
        # sort and index the files (using samtools)
        run_cmd = ' '.join(['samtools',
                            'collate',
                            '-o {0}'.format(output.outss),
                            '--reference {0}'.format(input.ref_path),
                            '{0}'.format(output.outm)])
        subprocess.call(run_cmd,shell=True)
        # run_cmd = ' '.join(['samtools',
        #                     'index',
        #                     '{0}'.format(output.outss)])
        # subprocess.call(run_cmd,shell=True)
        #
        # touch(output.outi)

rule filter_downsample_mask_bam:
    message:
        '''
        Merge files based on parameters.
        '''
    input:
        in1=path.join(config['out_dir'],"{sample}-collated.bam"),
        primer_amplicon_path=config['primer_amplicon_path'],
        ref_path=config['ref_fasta_path']
    output:
        outf=path.join(config['out_dir'],"{sample}-downsampled_trimmed.bam"),
        outss=path.join(config['out_dir'],"{sample}-downsampled_trimmed_sorted.bam"),
        outi=path.join(config['out_dir'],"{sample}-downsampled_trimmed_sorted.bam.bai")
    params:
        optional_arg_dict=config['filter_downsample_mask_bam'],
        sample_name="{sample}"
    threads:
        2
    run:
        import pandas as pd
        import subprocess
        from os import path
        script_path = path.join(config['src_dir'], 'trim_filter_bam_by_amplicon.py')
        main_cmd = 'python3 {0}'.format(script_path)
        required_arg_list = ['--in1={0}'.format(input.in1),
                             '--primer_amplicon_path={0}'.format(input.primer_amplicon_path),
                             '--outf={0}'.format(output.outf),
                             '--sample_name={0}'.format(params.sample_name)]
        run_cmd = create_cmd_req_list_optional_dict(main_cmd=main_cmd,
            required_arg_list=required_arg_list,
            optional_arg_dict=params.optional_arg_dict,
            delimiter='=')
        subprocess.call(run_cmd, shell=True)
        # downsample = args.downsample
        # reads_per_amplicon = args.reads_per_amplicon
        # max_del = args.max_del
        # max_ins = args.max_ins
        # max_sub = args.max_sub
        # mask_primers = args.mask_primers
        # out_dir = args.out_dir
        primer_groups = [1,2]

        if path.isfile(config['primer_amplicon_path']):
            df = pd.read_csv(config['primer_amplicon_path'])
            primer_groups =  list(df['PRIMER_GROUP'].unique())
        outf_files =[output.outf]
        outss_files = [output.outss]
        outi_files = [output.outi]

        for pg_i in primer_groups:
            touch('/Volumes/T7/covid_19_v4/out_amp/makenew_loop.txt')
            outf_file_i='{0}{1}.bam'.format(output.outf[:-4],pg_i)
            outf_files.insert(0,outf_file_i)
            outss_file_i='{0}_g{1}.bam'.format(output.outss[:-4],pg_i)
            outss_files.insert(0,outss_file_i)
            outi_file_i='{0}_g{1}.bam.bai'.format(output.outi[:-8],pg_i)
            outi_files.insert(0,outi_file_i)

        for outf, outss, outi in zip(outf_files, outss_files, outi_files):
            run_cmd = ' '.join(['samtools',
                                'sort',
                                '-o {0}'.format(outss),
                                '--reference {0}'.format(input.ref_path),
                                '{0}'.format(outf)])
            subprocess.call(run_cmd,shell=True)
            run_cmd = ' '.join(['samtools',
                                'index',
                                '{0}'.format(outss)])
            subprocess.call(run_cmd,shell=True)

            touch(outi)


rule bbmap_call_variants:
    message:
        '''
        BBMAP the call varinats from the mapped files
        '''
    input:
        in1=path.join(config['out_dir'],"{sample}-downsampled_trimmed_sorted.bam"),
        ref_path=config['ref_fasta_path'],
        config_path=config['config_path']
        # primer_amplicon_path=config['primer_amplicon_path']
    output:
        vcf=path.join(config['out_dir'],"{sample}.vcf"),
        ann_vcf=path.join(config['out_dir'],"{sample}_ann.vcf"),
        vcf_flag=path.join(config['out_dir'],"{sample}_vcf_primer_diff.csv")
    params:
        optional_arg_dict=config['bbmap_call_variants'],
        ram=config['ram'],
        ncbi_accession=config['ncbi_accession']
    threads:
        2
    run:
        import subprocess
        from os import path
        import pandas as pd
        primer_groups = []
        if path.isfile(config['primer_amplicon_path']):
            df = pd.read_csv(config['primer_amplicon_path'])
            primer_groups =  list(df['PRIMER_GROUP'].unique())
        ann_vcf_files =[output.ann_vcf]
        vcf_files = [output.vcf]
        bam_files = [input.in1]
        for pg_i in primer_groups:
            ann_vcf_file_i='{0}_g{1}_ann.vcf'.format(output.ann_vcf[:-8],pg_i)
            ann_vcf_files.insert(0,ann_vcf_file_i)
            vcf_file_i='{0}_g{1}.vcf'.format(output.vcf[:-4],pg_i)
            vcf_files.insert(0,vcf_file_i)
            bam_file_i='{0}_g{1}.bam'.format(input.in1[:-4],pg_i)
            bam_files.insert(0,bam_file_i)
        for vcf, ann_vcf, in1 in zip(vcf_files,ann_vcf_files,bam_files):
            if not path.isfile(in1):
                continue
            main_cmd = 'java -ea -Xmx{0}m -Xms{0}m -cp /bin/bbmap/current/ var2.CallVariants'.format(params.ram)
            required_arg_list = ['in={0}'.format(in1),
                                 'ref={0}'.format(input.ref_path),
                                 'vcf={0}'.format(vcf)]
            run_cmd = create_cmd_req_list_optional_dict(main_cmd=main_cmd,
                required_arg_list=required_arg_list,
                optional_arg_dict=params.optional_arg_dict,
                delimiter='=')
            subprocess.call(run_cmd,shell=True)

            # "minallelefraction": "0.001",
            # "rarity": "0.001",
            # "coverage": "t",
            # "calldel": "t",
            # "callins": "t",
            # "callsub": "t",
            # "mincov": "2",
            # "minreads": "7",
            # "minscore": "10",
            # "minquality": "10",
            # "overwrite": "t",
            # "minstrandratio": "0",
            # "minpairingrate": "0",
            # "usebias": "f",
            # "useedist": "f",
            # "useidentity": "f",
            # "minedist": "0",
            # "border": "0",
            # "trimq": "0",
            # "nscan": "f",
            # "minreadmapq": "0",
            # "minqualitymax": "10",
            # "covpenalty": "0"
            run_cmd = ' '.join(['snpEff',
                                '-c {0}'.format(input.config_path),
                                '-ud -onlyProtein',
                                '{0}'.format(params.ncbi_accession),
                                '{0}'.format(vcf),
                                '>',
                                '{0}'.format(ann_vcf)])
            subprocess.call(run_cmd,shell=True)
        import io
        import pandas as pd
        print(ann_vcf_files)
        if len(ann_vcf_files) >1:
            vcf1_path=ann_vcf_files[0]
            vcf2_path=ann_vcf_files[1]


            def read_vcf(path):
                with open(path, 'r') as f:
                    lines = [l for l in f if not l.startswith('##')]
                return pd.read_csv(
                    io.StringIO(''.join(lines)),
                    dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                           'QUAL': str, 'FILTER': str, 'INFO': str},
                    sep='\t'
                ).rename(columns={'#CHROM': 'CHROM'})

            df_vcf_1 = read_vcf(vcf1_path)
            df_vcf_1['GROUP'] = 1
            df_vcf_2= read_vcf(vcf2_path)
            df_vcf_2['GROUP'] = 2
            # df_pa = pd.read_csv(config['primer_amplicon_path'])
            df_pa2= df.groupby('AMPLICON')[['START_RIGHT','END_LEFT']].agg({'END_LEFT':min, 'START_RIGHT':max}).reset_index()
            df_pa2.sort_values(by = ['END_LEFT'], inplace=True)
            end_left = list(df_pa2['END_LEFT'])
            start_right = list(df_pa2['START_RIGHT'])
            end_left.pop(0)
            start_right.pop()
            range_list = [[x,y] for x,y in zip(end_left,start_right) if y - x > 0]

            col_list = ['POS', 'REF','ALT','QUAL','FILTER','INFO','GROUP']
            df_vcf = pd.concat([df_vcf_1[col_list],df_vcf_2[col_list]], ignore_index=True)
            position_list = list(set([x for x in df_vcf['POS']]))
            overlap_list = []
            for i in position_list:
                for range_i in range_list:
                    if i > range_i[0] and i <= range_i[1]:
                        overlap_list.append(i)
            df_vcf_ol = df_vcf[df_vcf['POS'].isin(overlap_list)]
            df_vcf_ol['COUNT'] = df_vcf_ol.groupby(['REF','ALT','POS'])['GROUP'].transform('count')
            df_vcf_flag = df_vcf_ol[df_vcf_ol['COUNT']<2]
            df_vcf_flag['AD'] = [float(x.split(';')[8][3:]) for x in df_vcf_flag['INFO']]
            df_vcf_flag['DP'] = [float(x.split(';')[9][3:]) for x in df_vcf_flag['INFO']]
            df_vcf_flag['AF'] = [float(x.split(';')[12][3:]) for x in df_vcf_flag['INFO']]
            df_vcf_flag.sort_values(by=['POS'], inplace=True)
            df_vcf_flag.to_csv(output.vcf_flag)
        touch(output.vcf_flag)

rule bbduk_filter_primer_reads:
    message:
        '''
        use bb duk to filter for reads with primer.  only verifies one primer if using matched pairs
        '''
    input:
        in1=path.join(config['input_dir'],"{sample}-R1.fastq.gz"),
        in2=path.join(config['input_dir'],"{sample}-R2.fastq.gz"),
        ref_adaptor_path=config['ref_adaptor_path']
    output:
        out=path.join(config['out_dir'],"{sample}-R1_trimmed.fastq.gz"),
        out2=path.join(config['out_dir'],"{sample}-R2_trimmed.fastq.gz")
    params:
        optional_arg_dict=config['bbduk_filter_primer_reads'],
        ram=config['ram'],
        ncbi_accession=config['ncbi_accession']
    threads:
        2
    run:
        import subprocess

        main_cmd = 'java -ea -Xmx{0}m -Xms{0}m -cp /bin/bbmap/current/ jgi.BBDuk'.format(params.ram)
        required_arg_list = ['in={0}'.format(input.in1),
                             'in2={0}'.format(input.in2),
                             'outm={0}'.format(output.out),
                             'outm2={0}'.format(output.out2),
                             'ref={0}'.format(input.ref_adaptor_path)]
        run_cmd = create_cmd_req_list_optional_dict(main_cmd=main_cmd,
            required_arg_list=required_arg_list,
            optional_arg_dict=params.optional_arg_dict,
            delimiter='=')
        print(run_cmd)
        subprocess.call(run_cmd,shell=True)

# k=${k} \
# hdist=${hdist} \
# rcomp=${rcomp} \
# minlength=${minlength_bbduk} \
# restrictleft=${restrictleft} \
# overwrite=t


rule bbduk_trim_primer_and_quality:
    message:
        '''
        BBMAP the call varinats from the mapped files
        '''
    input:
        in1=path.join(config['out_dir'],"{sample}-R1_trimmed.fastq.gz"),
        in2=path.join(config['out_dir'],"{sample}-R2_trimmed.fastq.gz"),
        ref_adaptor_path=config['ref_adaptor_path']
    output:
        out=path.join(config['out_dir'],"{sample}-R1_filtered_trimmed.fastq.gz"),
        out2=path.join(config['out_dir'],"{sample}-R2_filtered_trimmed.fastq.gz")
    params:
        optional_arg_dict=config['bbduk_trim_primer_and_quality'],
        ram=config['ram']
    threads:
        2
    run:
        import subprocess

        main_cmd = 'java -ea -Xmx{0}m -Xms{0}m -cp /bin/bbmap/current/ jgi.BBDuk'.format(params.ram)
        required_arg_list = ['in={0}'.format(input.in1),
                             'in2={0}'.format(input.in2),
                             'ref={0}'.format(input.ref_adaptor_path),
                             'out={0}'.format(output.out),
                             'out2={0}'.format(output.out2)]
        run_cmd = create_cmd_req_list_optional_dict(main_cmd=main_cmd,
            required_arg_list=required_arg_list,
            optional_arg_dict=params.optional_arg_dict,
            delimiter='=')
        subprocess.call(run_cmd,shell=True)

# with the filtered reads, trim the reads and low quality ends
# ktrim=${ktrim} \
# k=${k} \
# qtrim=${qtrim} \
# trimq=${trimq} \
# hdist=${hdist} \
# rcomp=${rcomp} \
# minlength=${minlength_bbduk} \
# restrictleft=${restrictleft} \
# overwrite=t
rule bbmerge_reads_from_bbduk_filter_trim:
    message:
        '''
        Merge files based on parameters.
        '''
    input:
        in1=path.join(config['out_dir'],"{sample}-R1_filtered_trimmed.fastq.gz"),
        in2=path.join(config['out_dir'],"{sample}-R2_filtered_trimmed.fastq.gz")
    output:
        out=path.join(config['out_dir'],"{sample}-merged.fastq.gz")
    params:
        optional_arg_dict={},
        ram=config['ram']
    threads:
        2
    run:
        print(config['SAMPLES'])
        print(['in={0}'.format(input.in1),
               'in2={0}'.format(input.in2),
               'out={0}'.format(output.out)])
        import subprocess

        main_cmd = 'java -ea -Xmx{0}m -Xms{0}m -Djava.library.path=/bin/bbmap/jni/ \
        -cp /bin/bbmap/current/ jgi.BBMerge'.format(params.ram)
        required_arg_list = ['in={0}'.format(input.in1),
                             'in2={0}'.format(input.in2),
                             'out={0}'.format(output.out)]
        run_cmd = create_cmd_req_list_optional_dict(main_cmd=main_cmd,
            required_arg_list=required_arg_list,
            optional_arg_dict=params.optional_arg_dict,
            delimiter='=')
        subprocess.call(run_cmd,shell=True)

rule slow_amplicon_normalization:
    message:
        '''
        Slower amplicon normaliztion use if needing to map to amplicon. 
        Might be useful if mapping a region that is similar/duplicated in the genome
        '''
    input:
        merged_fasta_path=path.join(config['out_dir'],"{sample}-merged.fastq.gz"),
        ref_path=config['ref_fasta_path'],
        bed_path=config['bed_path']
    output:
        out=path.join(config['out_dir'],"{sample}-downsampled.fastq")
    params:
        optional_arg_dict=config['slow_amplicon_normalization'],
        ram=config['ram']
    threads:
        2
    run:
        import subprocess
        script_path = path.join(config['src_dir'], 'amplicon_normalization.py')
        main_cmd = 'python3 {0}'.format(script_path)
        required_arg_list = ['--merged_fasta_path={0}'.format(input.merged_fasta_path),
                             '--ref_path={0}'.format(input.ref_path),
                             '--bed_path={0}'.format(input.bed_path),
                             '--output_filename={0}'.format(output.out),
                             '--ram={0}'.format(params.ram)]
        run_cmd = create_cmd_req_list_optional_dict(main_cmd=main_cmd,
            required_arg_list=required_arg_list,
            optional_arg_dict=params.optional_arg_dict,
            delimiter='=')
        subprocess.call(run_cmd,shell=True)

# --reads_per_amplicon=${reads_per_amplicon} \
# --outdir=${sample_dir} \
# --min_length_diff=${minlength_diff} \
# --include_primer=f \


rule bbmap_reformat_trim_read_length:
    message:
        '''
        Slower amplicon normaliztion use if needing to map to amplicon. 
        Might be useful if mapping a region that is similar/duplicated in the genome
        '''
    input:
        in1=path.join(config['out_dir'],"{sample}-downsampled.fastq"),
    output:
        out=path.join(config['out_dir'],"{sample}-downsampled_trimmed.fastq")
    params:
        optional_arg_dict=config['bbmap_reformat_trim_read_length'],
        ram=config['ram']
    threads:
        2
    run:
        import subprocess

        main_cmd = 'java -ea -Xms{0}m -cp /bin/bbmap/current/ jgi.ReformatReads'.format(params.ram)
        required_arg_list = ['in={0}'.format(input.in1),
                             'out={0}'.format(output.out)]
        run_cmd = create_cmd_req_list_optional_dict(main_cmd=main_cmd,
            required_arg_list=required_arg_list,
            optional_arg_dict=params.optional_arg_dict,
            delimiter='=')
        subprocess.call(run_cmd,shell=True)
# maxlength=${maxlength} \
# minlength=${minlength} \
# overwrite=t
rule bbmap_filtered_trimmed:
    message:
        '''
        BBMAP the filtered, trimmed reads, 
        sort and index bam with samtools
        '''
    input:
        in1=path.join(config['out_dir'],"{sample}-downsampled_trimmed.fastq"),
        ref_path=config['ref_fasta_path']
    output:
        outm=path.join(config['out_dir'],"{sample}-downsampled_trimmed.bam"),
        outss=path.join(config['out_dir'],"{sample}-downsampled_trimmed_sorted.bam"),
        outi=path.join(config['out_dir'],"{sample}-downsampled_trimmed_sorted.bam.bai")
    params:
        optional_arg_dict=config['bbmap_filtered_trimmed'],
        ram=config['ram']
    threads:
        2
    run:
        import subprocess
        main_cmd = 'java -ea -Xmx{0}m -Xms{0}m -cp /bin/bbmap/current/ align2.BBMap build=1'.format(params.ram)
        required_arg_list = ['in={0}'.format(input.in1),
                             'ref={0}'.format(input.ref_path),
                             'outm={0}'.format(output.outm)]
        run_cmd = create_cmd_req_list_optional_dict(main_cmd=main_cmd,
            required_arg_list=required_arg_list,
            optional_arg_dict=params.optional_arg_dict,
            delimiter='=')
        subprocess.call(run_cmd,shell=True)
        # maxindel=100 \
        # overwrite=t \
        # nodisk=t
        # sort and index the files (using samtools)
        run_cmd = ' '.join(['samtools',
                            'sort',
                            '-o {0}'.format(output.outss),
                            '--reference {0}'.format(input.ref_path),
                            '{0}'.format(output.outm)])
        subprocess.call(run_cmd,shell=True)
        run_cmd = ' '.join(['samtools',
                            'index',
                            '{0}'.format(output.outss)])
        subprocess.call(run_cmd,shell=True)

        touch(output.outi)

rule create_primer_fasta_ref_bed:
    message:
        '''
        Slower amplicon normaliztion use if needing to map to amplicon. 
        Might be useful if mapping a region that is similar/duplicated in the genome
        '''
    input:
        ref_fasta_path=config['ref_fasta_path'],
        bed_path=config['bed_path']
    output:
        ref_adaptor_path=config['ref_adaptor_path'],
        df_bed_temp=path.join(config['out_dir'], 'df_bed.csv')
    threads:
        2
    run:
        from Bio.Seq import Seq
        from Bio import SeqIO
        import pandas as pd
        from Bio.Seq import Seq
        import os
        from pathlib import Path
        bed_path = input.bed_path

        ref_fasta_path = input.ref_fasta_path

        ref_adaptor_path = output.ref_adaptor_path

        df_bed_temp=output.df_bed_temp
        if not os.path.isfile(bed_path):
            print("Error Bedfile does not exist")
            exit(1)
        fasta_sequences = SeqIO.parse(open(ref_fasta_path), 'fasta')

        # only one sequence so
        name = None
        sequence = None
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            break
        if sequence is None:
            print('Incoming fasta File is blank or no valid sequences')
            exit(0)

        df_bed = pd.read_csv(bed_path, sep='\t', header=None, usecols=[0, 1, 2, 3, 4],
            names=['SEQ_ID', 'START', 'END', 'AMPLICON_ID','PRIMER_GROUP'])

        df_bed['SIDE'] = [x.split('_')[2] for x in df_bed['AMPLICON_ID']]
        df_bed['AMPLICON'] = ['_'.join(x.split('_')[:2]) for x in df_bed['AMPLICON_ID']]
        df_bed['GROUP']  = [x.split('_')[1] for x in df_bed['AMPLICON_ID']]
        df_bed['PRIMER_SEQUENCE'] = [sequence[int(x):int(y)] for x, y in zip(df_bed['START'], df_bed['END'])]
        # take the reverse complement so it can match directly.
        df_bed['PRIMER_SEQUENCE'] = [x if y.upper() == 'LEFT' else str(Seq(x).reverse_complement()) for x, y in zip(df_bed['PRIMER_SEQUENCE'], df_bed['SIDE'])]
        with open(ref_adaptor_path, 'w') as out_file:
            for x,y in zip(df_bed['AMPLICON_ID'], df_bed['PRIMER_SEQUENCE']):
                out_file.write('>{0}\n'.format(x))
                out_file.write('{0}\n'.format(y))
        df_bed.to_csv(df_bed_temp, index=False)


rule create_primer_amp_from_ref_adaptor_fasta:
    message:
        '''
        Slower amplicon normaliztion use if needing to map to amplicon. 
        Might be useful if mapping a region that is similar/duplicated in the genome
        '''
    input:
        ref_fasta_path=config['ref_fasta_path'],
        df_bed_temp=path.join(config['out_dir'], 'df_bed.csv')
    output:
        primer_amplicon_path=config['primer_amplicon_path']
    threads:
        2
    run:
        import pandas as pd
        from Bio import SeqIO

        primer_amplicon_path = output.primer_amplicon_path
        ref_fasta_path = input.ref_fasta_path
        df_bed_temp=input.df_bed_temp
        fasta_sequences = SeqIO.parse(open(ref_fasta_path), 'fasta')
        name = None
        sequence = None
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            break
        df_bed_side = {}
        df_bed = pd.read_csv(df_bed_temp)
        # Cartesian join instead of pivot because some have multiple primer combinations.
        df_bed2 = df_bed.drop(columns=['AMPLICON_ID','GROUP'])
        for side_i in ['LEFT', 'RIGHT']:
            df_bed_side[side_i]  = df_bed2[df_bed2['SIDE']== side_i]
            df_bed_side[side_i].drop(columns= ['SIDE'], inplace=True)
            df_bed_side[side_i].rename(columns={'PRIMER_SEQUENCE':'PRIMER_SEQUENCE_{0}'.format(side_i),
                                                'START':'START_{0}'.format(side_i),
                                                'END':'END_{0}'.format(side_i)}, inplace=True)
        df_bed_pivot = df_bed_side['RIGHT'].merge(df_bed_side['LEFT'], on = ['SEQ_ID', 'AMPLICON','PRIMER_GROUP'])
        df_bed_pivot['AMPLICON_SEQ'] = [sequence[x:y] for x,y in zip(df_bed_pivot['END_LEFT'],df_bed_pivot['START_RIGHT'] )]
        df_bed_pivot['RANK'] = df_bed_pivot.groupby(['AMPLICON'])['AMPLICON'].cumcount().add(1)
        df_bed_pivot['AMPLICON_NAME'] = ['{0}_{1}'.format(x,y) for x,y in zip(df_bed_pivot['AMPLICON'], df_bed_pivot['RANK'])]
        # df_bed_pivot[df_bed_pivot['AMPLICON'] == 'nCoV-2019_9']


        df_bed_pivot.to_csv(primer_amplicon_path,index=False)

rule create_insert_fasta_path_from_primer_amplicon:
    message:
        '''
        Slower amplicon normaliztion use if needing to map to amplicon. 
        Might be useful if mapping a region that is similar/duplicated in the genome
        '''
    input:
        primer_amplicon_path=config['primer_amplicon_path']
    output:
        insert_fasta_path=config['insert_fasta_path']
    threads:
        2
    run:
        primer_amplicon_path=input.primer_amplicon_path
        import pandas as pd
        df_bed_pivot = pd.read_csv(primer_amplicon_path)
        with open(output.insert_fasta_path, 'w') as out_file:
            for x,y in zip(df_bed_pivot['AMPLICON_NAME'], df_bed_pivot['AMPLICON_SEQ']):
                out_file.write('>{0}\n'.format(x))
                out_file.write('{0}\n'.format(y))
