#! /usr/bin/env python3

import os
import argparse
import subprocess
import pandas as pd

# Required dependencies in $PATH include bwa, samtools, bbmap, bedtools, bcftools, usearch

"""Variables for database pathways; these need to be set until pipeline upgrade"""
adapter_database1 = '/data/opt/jarfiles/trimmomatic/v0.39/Trimmomatic-0.39/adapters/All_adapters.fa'

"""Full pathway for trimmomatic jar. Not sure how to include in sysargs argument at this point"""

trimmomatic_path1 = '/data/opt/jarfiles/trimmomatic/v0.39/Trimmomatic-0.39/trimmomatic-0.39.jar'

def create_directory(outdir):
    """Function to create the output directory"""
    if os.path.isdir(outdir):
        raise Exception("Directory already exists")
    if not os.path.isdir(outdir):
        os.system("mkdir -p %s" % outdir)
    return

def trimmomatic_command(trimmomatic_path, threads, sample_name, outdir, read1, read2, adapter_database):
    trimmomatic_cwd = ['java', '-jar', trimmomatic_path, 'PE', '-phred33',
                       '-threads', threads, read1, read2, '{0}/{1}_ILMN_R1.fastq.gz'.format(outdir, sample_name),
                       '{0}/{1}_ILMN_unpaired_R1.fastq.gz'.format(outdir, sample_name),
                       '{0}/{1}_ILMN_R2.fastq.gz'.format(outdir, sample_name),
                       '{0}/{1}_ILMN_unpaired_R2.fastq.gz'.format(outdir, sample_name),
                       'ILLUMINACLIP:{0}:2:30:10'.format(adapter_database), 'LEADING:3', 'TRAILING:15',
                       'SLIDINGWINDOW:4:15', 'MINLEN:36']
    subprocess.run(trimmomatic_cwd)
    return

def get_arguments():
    """Parse assembler arguments"""
    parser = argparse.ArgumentParser(description="Variant identification and copy number quantification", add_help=False)
    # Help arguments
    help_group = parser.add_argument_group("Help")
    help_group.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    help_group.add_argument('-V', '--version', action='version', version='%(prog)s version 0.1',
                            help="Show convict's version number")
    # input_arguments
    input_group = parser.add_argument_group("Inputs")

    input_group.add_argument('-R1', '--read1', required=True, help='forward read', type=str, default=None)
    input_group.add_argument('-R2', '--read2', required=True, help='reverse read', type=str, default=None)
    input_group.add_argument('-db', '--database', required=True, help='Database subject for query search', type=str, default=None)
    input_group.add_argument('-c', '--control', required=True, help='Control genes for gene depth normalization', type=str, default=None)

    # optional_arguments
    optional_group = parser.add_argument_group("Optional Inpugts")
    optional_group.add_argument('-s', '--sample_name', required=False, help="Sample name prefix", type=str,
                                default='SAMPLE')
    optional_group.add_argument('-t', '--threads', required=False, help='Number of threads', type=str, default='1')
    optional_group.add_argument('-p', '--ploidy', required=False, help='Ploidy of sample for variant calling. '
                                'NB: For Diploid samples, use "GRCh37"', type=str, default='1')
    optional_group.add_argument('-tr', '--trim', action='store_true')
    # output_arguments
    output_group = parser.add_argument_group("Outputs")
    output_group.add_argument('-o', '--outdir', required=True, help="Name of the output directory", type=str, default=None)

    args = parser.parse_args()
    return args

def run_conditions():
    args = get_arguments()
    database = args.database
    control = args.control
    threads = args.threads
    sample_name = args.sample_name
    if args.outdir:
        create_directory(args.outdir)
    else:
        args.outdir = os.getcwd()
    outdir = args.outdir
    tmp = '{0}/tmp'.format(outdir)
    create_directory(tmp)

    # Trim input reads if 'trim' option is provided
    if args.trim is True:
        print("\nTrimming low-quality bases and adapters from raw Illumina short-reads\n")
        trimmomatic_command(trimmomatic_path1, threads, sample_name, outdir, args.read1, args.read2, adapter_database1)
        os.system('rm %s' % ('{0}/*unpaired*'.format(outdir)))
        trim_R1 = '{0}/{1}_ILMN_R1.fastq.gz'.format(outdir, sample_name)
        trim_R2 = '{0}/{1}_ILMN_R2.fastq.gz'.format(outdir, sample_name)
    else:
        print("\nUser has provided trimmed, paired-end short-reads\n")
        trim_R1 = args.read1
        trim_R2 = args.read2
    # Perform usearch cluster for target DNA sequences and consensus fasta generation of target DB
    print("\nCluster target DNA sequences\n")
    clusters = '{0}/clusters.tsv'.format(outdir)
    consensus_db = '{0}/targetDB_nt_cons.fasta'.format(tmp)
    subprocess.Popen('usearch -cluster_fast ' + database + ' -id 0.9 -uc ' + clusters +
                     ' -consout ' + consensus_db, shell=True).wait()
    # Create bwa mem alignment to generate coverage depth and breadth of coverage stats from bbmap pileup.sh
    print("\nBWA mem alignment of trimmed short-reads to consensus fasta file from uclust results\n")
    subprocess.Popen('bwa index ' + consensus_db, shell=True).wait()
    subprocess.Popen('bwa mem -t ' + threads + ' ' + consensus_db + ' ' + trim_R1 + ' ' + trim_R2 +
                     ' | samtools sort -@ ' + threads + ' > ' + tmp + '/database_align_sort.bam ', shell=True).wait()
    target_covstats = '{0}/database_align_sort_covstats.tsv'.format(tmp)
    subprocess.Popen('pileup.sh delcoverage=f in=' + tmp + '/database_align_sort.bam ' +
                     'out=' + target_covstats, shell=True).wait()
    # Create bwa mem alignment to generate coverage depth for housekeeping genes to normalize coverage for target genes
    subprocess.Popen('bwa index ' + control, shell=True).wait()
    subprocess.Popen('bwa mem -t ' + threads + ' ' + control + ' ' + trim_R1 + ' ' + trim_R2 +
                     ' | samtools sort -@ ' + threads + ' ' + '> ' + tmp + '/control_align_sort.bam ', shell=True).wait()
    cntrl_covstats = '{0}/control_align_sort_covstats.tsv'.format(tmp)
    subprocess.Popen('pileup.sh delcoverage=f in=' + tmp + '/control_align_sort.bam ' +
                     'out=' + cntrl_covstats, shell=True).wait()
    # Convert target and cntrl covstats as well as clusters tab file into pandas dataframes
    print("\nGenerate breadth of coverage, coverage depth, and coverage depth ratios to determine copy number\n")
    header1 = ['ID', 'Avg_fold', 'Length', 'Ref_GC', 'Covered_percent', 'Covered_bases', 'Plus_reads',
                     'Minus_reads', 'Read_GC', 'Median_fold', 'Std_Dev']
    target_df = pd.read_table(target_covstats, header=0, names=header1)
    cntrl_df = pd.read_table(cntrl_covstats, header=0, names=header1)
    cluster_header = ['record_type', 'ID', 'seq_length', 'percent_id', 'strandness', 'na1', 'na2',
                      'compr_align', 'cluster_label', 'target_label']
    clusters_df = pd.read_table(clusters, names=cluster_header)
    # Subset clusters file into only cluster record types
    clusters_subset = clusters_df[clusters_df.record_type.str.contains('C', case=False)]
    # Append 'Cluster' prefix to cluster number to match target_df output
    clusters_subset.loc[clusters_subset.index, 'ID'] = 'Cluster' + clusters_subset['ID'].astype(str)
    # Append representative cluster label sequence to target df
    target_df = pd.merge(target_df, clusters_subset[['ID', 'cluster_label']], how='left', on='ID')
    # Filter cntrl dataframe values if 'Avg_fold' is less than < 10 or 'Covered_percent' < 90.0
    cntrl_filter_drop = cntrl_df[(cntrl_df['Avg_fold'] < 10) | (cntrl_df['Covered_percent'] < 90.0)]
    cntrl_filter_drop_outfile = '{0}/{1}_cntrl_filter_drop_covstats.tsv'.format(tmp, sample_name)
    cntrl_filter_drop.to_csv(cntrl_filter_drop_outfile, sep="\t", index=False)
    # Filter cntrl dataframe values if 'Avg_fold' is greater than >= 10 or 'Covered_percent' >= 90.0
    cntrl_filter_df = cntrl_df[(cntrl_df['Avg_fold'] >= 10) & (cntrl_df['Covered_percent'] >= 90.0)]
    # Calculate Avg_fold * length and create new column
    cntrl_filter_df['avg_fold_X_length'] = cntrl_filter_df['Avg_fold'] * cntrl_filter_df['Length']
    # Extract average_fold from control dataframe
    average_fold_cntrl = cntrl_filter_df['avg_fold_X_length'].sum() / cntrl_filter_df['Length'].sum()
    # Append average_fold to target dataframe
    target_df["Avg_fold_cntrl"] = average_fold_cntrl
    # Calculate [target_gene_Avg_fold / cntrl_gene_Avg_fold ] if df["Covered_percent"] > 90
    target_df.loc[target_df["Covered_percent"] >= 90, 'norm_cov_depth'] = target_df["Avg_fold"] / target_df[
        "Avg_fold_cntrl"]
    outfile = '{0}/{1}_convict_results.tsv'.format(outdir, sample_name)
    # Reformat columns for output and remove unnecessary columns
    target_df = target_df[['ID', 'cluster_label', 'Length', 'Covered_percent', 'Avg_fold', 'Median_fold', 'Std_Dev',
                           'Avg_fold_cntrl', 'norm_cov_depth', 'Covered_bases', 'Plus_reads', 'Minus_reads', 'Read_GC']]
    target_df.to_csv(outfile, sep="\t", index=False)
    # Call variants using bcftools and normalize variants for ambiguous REF/ALT sites
    print("\nVariant calling on consensus cluster fasta files to get consensus fasta file from short-read alignment\n")
    subprocess.Popen('bcftools mpileup -Ou --threads ' + threads + ' -f ' + consensus_db + ' ' +
                     tmp + '/database_align_sort.bam' + ' | bcftools call -mv -Ou --threads ' + threads + ' --ploidy' +
                     ' ' + args.ploidy + ' | bcftools norm -Oz -f ' + consensus_db + ' --threads ' + threads + ' -o ' +
                     tmp + '/norm_calls.vcf.gz', shell=True).wait()
    subprocess.Popen('bcftools index ' + tmp + '/norm_calls.vcf.gz', shell=True).wait()
    # Need to mask reference.fasta for regions with zero coverage before creating consensus fasta file
    subprocess.Popen('bedtools genomecov -bga -split -ibam ' + tmp + '/database_align_sort.bam' + ' > ' +
                     tmp + '/coverage.bed', shell=True).wait()
    subprocess.Popen('grep -w 0$ ' + tmp + '/coverage.bed' + ' > ' + tmp + '/zero-coverage.bed', shell=True).wait()
    subprocess.Popen('bedtools maskfasta -fi ' + consensus_db + ' -bed ' + tmp + '/zero-coverage.bed' + ' -fo ' +
                     tmp + '/ref_zero_masked.fasta', shell=True).wait()
    # Create consensus fasta file from alignment
    subprocess.Popen('bcftools consensus --fasta-ref ' + tmp + '/ref_zero_masked.fasta' + ' -o ' +
                     outdir + '/consensus_align.fasta' + ' ' + tmp + '/norm_calls.vcf.gz', shell=True).wait()
    blastout = '{0}/blast_results_tmp.tsv'.format(tmp)
    consensus_align_fasta = '{0}/consensus_align.fasta'.format(outdir)
    blastn_cmd = 'blastn -subject {0} -query {1} -outfmt \"6 std qcovs\" -out {2} -evalue 0.00001 -perc_identity 0.9' \
                 ' -qcov_hsp_perc 0.9 -culling_limit 1'.format(database, consensus_align_fasta, blastout)
    subprocess.Popen(blastn_cmd, shell=True).wait()
    results = pd.read_csv(blastout, sep="\t", header=None)
    headers = ['query', 'subject',
               'pc_identity', 'aln_length', 'mismatches', 'gaps_opened',
               'query_start', 'query_end', 'subject_start', 'subject_end',
               'e_value', 'bitscore', 'qcovs']
    results.columns = headers
    blast_results = '{0}/{1}_blastn_results.tsv'.format(outdir, sample_name)
    results.to_csv(blast_results, sep="\t", index=False)
    print("\nCONVICT pipeline complete!\n")

if __name__ == '__main__': run_conditions()

