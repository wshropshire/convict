#! /usr/bin/env python3

import os
import argparse
import subprocess
import pandas as pd
import numpy as np
from Bio import SeqIO

'''Required dependencies in $PATH include kma, kmerresistance, bowtie2, minimap2 (ONT), samtools, bbmap, bedtools, bcftools
This version of CONVICT will default to using kmerresistance with default species and ResFinder database
NB: If using custom databases for target/species databases, make sure they are indexed using kma index
NB: Make sure that headers of database fasta files do not include spaces as this will affect SeqIO matching of kmerresistance hits'''

"""Variables for kma database and trimmomatic database pathways"""
dir = os.path.dirname(__file__)
resfinder_database1 = os.path.join(dir, './../db/resfinder_db')
adapter_database1 = os.path.join(dir, './../db/ALL_adapters.fa')

# Bacteria database pathway has to be specified due to size constraints
species_database1 = '~/Documents/Datbases/kma_databases/bacteria.ATG'

"""Full pathway for trimmomatic jar. Must be modified for users path environment if trimming is intended to be included in pipeline.
Future versions will include this, likely in a docker container"""
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
    help_group.add_argument('-V', '--version', action='version', version='%(prog)s version 1.0',
                            help="Show convict's version number")
    # input_arguments
    input_group = parser.add_argument_group("Inputs")

    input_group.add_argument('-l', '--long_reads', required=False, help='long-read input', type=str, default=None)
    input_group.add_argument('-R1', '--read1', required=False, help='forward read', type=str, default=None)
    input_group.add_argument('-R2', '--read2', required=False, help='reverse read', type=str, default=None)
    input_group.add_argument('-c', '--control', required=True, help='Control genes for gene depth normalization', type=str, default=None)
    input_group.add_argument('-db', '--database', required=False, help='Database subject for query search', type=str,
                             default=resfinder_database1)
    input_group.add_argument('-sp', '--species', required=False, help='Species database for identification with short-reads', type=str,
                             default=species_database1)

### Need to include better error handling to specify that either R1/R2 or long reads are required. 

    # optional_arguments
    optional_group = parser.add_argument_group("Optional Inputs")
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
    control = args.control
    threads = args.threads
    database = args.database
    species = args.species
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
        print("\n# Trimming low-quality bases and adapters from raw Illumina short-reads")
        trimmomatic_command(trimmomatic_path1, threads, sample_name, outdir, args.read1, args.read2, adapter_database1)
        os.system('rm %s' % ('{0}/*unpaired*'.format(outdir)))
        trim_R1 = '{0}/{1}_ILMN_R1.fastq.gz'.format(outdir, sample_name)
        trim_R2 = '{0}/{1}_ILMN_R2.fastq.gz'.format(outdir, sample_name)
    elif args.read1 is not None:
        print("\n# User has provided trimmed, paired-end short-reads")
        trim_R1 = args.read1
        trim_R2 = args.read2
    else:
        print("\n# No short-read data is provided for trimming")

#### Adding in kmerresistance step to identify variants from short-reads which removes the user's need to specify a target database prior. 
#### i.e. can include large gene database to look for hits, similar to abricate, using kmerresistance  by first indexing, and use the results as input to look for 
#### gene amplifications or other gene variants with the samtools/bcftools variant calling
    
    kmerresistance_out = '{0}/{1}_kmerresistance'.format(tmp, sample_name)
    if args.long_reads is None:
        print("\n# Performing k-mer alignment using kmerresistance of short-reads to %s database" %(os.path.basename(database)))
        print("\n# Includes k-mer alignment using kmerresistance of short-reads to %s database for species identification\n" %(os.path.basename(species)))
        subprocess.Popen('kmerresistance ' '-i ' + trim_R1 + ' ' + trim_R2 + ' -o ' + 
            kmerresistance_out + ' -t_db ' + database + ' -s_db ' + species, shell=True).wait()
    else:
        print("\n# Performing k-mer alignment using kmerresistance of long-reads to %s database" %(os.path.basename(database)))
        print("\n# Includes k-mer alignment using kmerresistance of long-reads to %s database for species identification\n" %(os.path.basename(species)))
        subprocess.Popen('kmerresistance ' '-i ' + args.long_reads + ' -o ' + 
            kmerresistance_out + ' -t_db ' + database + ' -s_db ' + species, shell=True).wait()
    # Need to subset fsa file with hits that hit criterion for reporting on KmerRes file
    KmerRes_file = '{0}/{1}_kmerresistance.KmerRes'.format(tmp, sample_name)
    df = pd.read_table(KmerRes_file)
    list1 = df.iloc[1: , :1]
    set_id = set(list1['#Template'])
    KmerRes_fasta_out = '{0}/{1}_KmerRes.fasta'.format(tmp, sample_name)
    fasta_file = '{0}/{1}_kmerresistance.fsa'.format(tmp, sample_name)
    fasta_seq = SeqIO.parse(open(fasta_file), 'fasta')
    with open(KmerRes_fasta_out, "w") as f:
        for seq in fasta_seq:
            if seq.id in set_id:
                SeqIO.write([seq], f, "fasta")
    # Perform local alignment of kmerresistance database hits with short or long-reads using bowtie2
    if args.long_reads is None:
        print("\n# Bowtie2 alignment with paired-end short-reads of kmerresistance output\n")
        subprocess.Popen('bowtie2-build -f ' + KmerRes_fasta_out + ' target', shell=True).wait()
        subprocess.Popen('bowtie2 -p ' + threads + ' --very-sensitive-local -a -t -x target -1 ' + trim_R1 + ' -2 ' + trim_R2 +
                         ' |  samtools sort -@ ' + threads + ' > ' + tmp + '/target_align_sort.bam ', shell=True).wait()
    else:
        # Minimap2 by default already outputs 5 secondary alignments. This can be tuned using -N option.
        print("\n# Minimap2 alignment of ONT long-reads to kmerresistance output\n")
        subprocess.Popen('minimap2 -t ' + threads + ' -ax map-ont ' + KmerRes_fasta_out + ' ' + args.long_reads +
                         ' | samtools sort -@ ' + threads + ' > ' + tmp + '/target_align_sort.bam ', shell=True).wait()
    target_bam = '{0}/target_align_sort.bam'.format(tmp)
# Create pubMLST/control gene alignment to generate coverage depth of housekeeping genes to normalize coverage for target genes
# implementing same 'input all' alignments as for Kmerresistance Resfinder output alignment, however, shouldn't make significant difference; should be kept consistent
# for the purpose of keeping same parameters the same the calculate normalized coverage depth. Good to control for housekeeping genes that have variable coverage depths
    if args.long_reads is None:
        print("\n# Bowtie2 alignment with paired-end short-reads of control genes\n")
        subprocess.Popen('bowtie2-build -f ' + control + ' control', shell=True).wait()
        subprocess.Popen(
            'bowtie2 -p ' + threads + ' --very-sensitive-local -a -t -x control -1 ' + trim_R1 + ' -2 ' + trim_R2 +
            ' |  samtools sort -@ ' + threads + ' > ' + tmp + '/control_align_sort.bam ', shell=True).wait()
    else:
        print("\n# Minimap2 alignment of ONT long-reads to control gene fasta file\n")
        subprocess.Popen('minimap2 -t ' + threads + ' -ax map-ont ' + control + ' ' + args.long_reads +
                         ' | samtools sort -@ ' + threads + ' ' + '> ' + tmp + '/control_align_sort.bam ',
                         shell=True).wait()
    cntrl_bam = '{0}/control_align_sort.bam'.format(tmp)
# Create output file for binned covstats
    control_bins = '{0}/control_covbin.tsv'.format(tmp)
    control_bins_corrected = '{0}/control_covbin_corrected.tsv'.format(tmp)
    control_covstats = '{0}/control_covstats.tsv'.format(tmp)
    # Use pileup.sh to calculate coverage depth for control genes and output 'BED' like file with coverage depth/nt position
    # stdev=f to avoid mean/stdev printout on bincov
    subprocess.Popen('pileup.sh delcoverage=f binsize=100 stdev=f in=' + cntrl_bam + ' bincov=' + control_bins + ' out=' + control_covstats, shell=True).wait()
    # Convert target to pandas dataframe and add headers to binned output
    control_bins_df = pd.read_table(control_bins)
    control_bins_df.columns = ['RefName', 'Cov', 'Pos', 'RunningPos']
    # Create list of unique control genes
    gene_names = control_bins_df.RefName.unique()
    # Create dataframes
    # Might want to take control_final out
    control_final = pd.DataFrame()
    control_stats_final = pd.DataFrame()
    for x in range(int(len(gene_names))):
        # select observations for each respective control gene
        control = np.where(control_bins_df.RefName == gene_names[x])
        control = pd.DataFrame(control_bins_df.loc[control])
        # calculate average, standard deviation, and distance one standard deviation above/below average
        avg = np.mean(control["Cov"])
        stdev = np.std(control["Cov"])
        above = avg + stdev
        below = avg - stdev
        # calculate a coefficient of variance
        cv = stdev / avg
        # calculate absolute value of distance of each bin from the average of all bins
        d = np.abs(control.Cov - avg)
        control = pd.concat([control, d], axis=1)
        control.columns = ['RefName', 'Cov', 'Pos', 'RunningPos', 'distance']
        # Create a dataframe that resembles covstats
        control_stats = pd.DataFrame(
            {'Control_Gene_Name': [gene_names[x]], 'Avg_fold': [avg], 'Std_dev': [stdev],
             'CV': [cv], 'stdev_below': [below], 'stdev_above': [above]})
        # Loop through control gene stats to find
        for z in range(int(len(control.index))):
            # Append bin to final control stats dataframe if CV < 0.05, else remove windows above/below 1 stdev
            if (control_stats["CV"] < 0.05).all():
                control_stats_final = pd.concat([control_stats_final, control_stats])
            # If CV not equal to 0, and greater than or equal to 0.05 determine how many windows are above/below 1 std dev.
            # Check the all argument
            elif (control_stats["Avg_fold"] != 0).all():
                windows_below = pd.DataFrame(control[control.Cov <= below])
                windows_above = pd.DataFrame(control[control.Cov >= above])
                # Sort windows by distance
                windows = pd.concat([windows_below, windows_above]).sort_values(by="distance", ascending=False)
                # If there are windows drop the indexed window
                if not windows.empty:
                    control = control.drop([windows.index[0]])
            # Re-calculate average, standard deviation, cv
            avg = np.mean(control["Cov"])
            stdev = np.std(control["Cov"])
            above = avg + stdev
            below = avg - stdev
            cv = stdev / avg
            control_stats = pd.DataFrame(
                {'Control_Gene_Name': [gene_names[x]], 'Avg_fold': [avg], 'Std_dev': [stdev],
                 'CV': [cv], 'stdev_below': [below], 'stdev_above': [above]})
        control_final = pd.concat([control_final, control])
# Why would duplicates be present
    control_stats_final = control_stats_final.drop_duplicates(subset=['Control_Gene_Name'])
    control_final.to_csv(control_bins_corrected, sep='\t')

# Need to do same algorithm for target genes:
    # Create bin covstats for target genes
    target_bins = '{0}/target_covbin.tsv'.format(tmp)
    target_bins_corrected = '{0}/target_covbin_corrected.tsv'.format(tmp)
    target_covstats = '{0}/target_covstats.tsv'.format(tmp)
    # Use pileup.sh to calculate coverage depth for target genes and output 'BED' like file with coverage depth/nt position
    subprocess.Popen('pileup.sh delcoverage=f binsize=100 stdev=f in=' + target_bam + ' bincov=' + target_bins + ' out=' + target_covstats,
                     shell=True).wait()
    # Convert target to pandas dataframe and add headers to binned output
    target_bins_df = pd.read_table(target_bins)
    target_bins_df.columns = ['RefName', 'Cov', 'Pos', 'RunningPos']
    # Create list of unique target genes
    gene_names = target_bins_df.RefName.unique()
    # Create dataframes
    target_final = pd.DataFrame()
    target_stats_final = pd.DataFrame()
    for x in range(int(len(gene_names))):
        # select observations for each respective target gene
        target = np.where(target_bins_df.RefName == gene_names[x])
        target = pd.DataFrame(target_bins_df.loc[target])
        # calculate average, standard deviation, and distance one standard deviation above/below average
        avg = np.mean(target["Cov"])
        stdev = np.std(target["Cov"])
        above = avg + stdev
        below = avg - stdev
        # calculate a coefficient of variance
        cv = stdev / avg
        # calculate absolute value of distance of each bin from the average of all bins
        d = np.abs(target.Cov - avg)
        target = pd.concat([target, d], axis=1)
        target.columns = ['RefName', 'Cov', 'Pos', 'RunningPos', 'distance']
        target_stats = pd.DataFrame(
            {'Target_Gene_Name': [gene_names[x]], 'Avg_fold': [avg], 'Std_dev': [stdev],
             'CV': [cv], 'stdev_below': [below], 'stdev_above': [above]})
        # Loop through target gene stats
        for z in range(int(len(target.index))):
            # Append bin to final control stats dataframe if CV < 0.05, else remove windows above/below 1 stdev
            if (target_stats["CV"] < 0.05).all():
                target_stats_final = pd.concat([target_stats_final, target_stats])
            # If CV not equal to 0, and greater than or equal to 0.05 determine how many windows are above/below 1 std dev.
            elif (target_stats["Avg_fold"] != 0).all():
                windows_below = pd.DataFrame(target[target.Cov <= below])
                windows_above = pd.DataFrame(target[target.Cov >= above])
                # Sort windows by distance
                windows = pd.concat([windows_below, windows_above]).sort_values(by="distance", ascending=False)
                # If there are windows drop the indexed windows
                if not windows.empty:
                    target = target.drop([windows.index[0]])
            # Re-calculate average, standard deviation, cv
            avg = np.mean(target["Cov"])
            stdev = np.std(target["Cov"])
            above = avg + stdev
            below = avg - stdev
            cv = stdev / avg
            target_stats = pd.DataFrame(
                {'Target_Gene_Name': [gene_names[x]], 'Avg_fold': [avg], 'Std_dev': [stdev],
                 'CV': [cv], 'stdev_below': [below], 'stdev_above': [above]})
        target_final = pd.concat([target_final, target])
        # Why would duplicates be present
    target_stats_final = target_stats_final.drop_duplicates(subset=['Target_Gene_Name'])
    # Create pandas dataframes for covstats information to append Length and Covered_percent to both respective final stats databases
    print("\nGenerate breadth of coverage, coverage depth, and coverage depth ratios to determine copy number\n")
    header1 = ['ID', 'Avg_fold', 'Length', 'Ref_GC', 'Covered_percent', 'Covered_bases', 'Plus_reads',
                     'Minus_reads', 'Read_GC', 'Median_fold', 'Std_Dev']
    target_df = pd.read_table(target_covstats, header=0, names=header1)
    cntrl_df = pd.read_table(control_covstats, header=0, names=header1)
    # drop 'Avg_fold', 'Median_fold', 'Plus reads', 'Minus_reads', 'Read_GC', 'Median_fold', 'Std_Dev'
    target_df = target_df.drop(['Avg_fold', 'Plus_reads', 'Minus_reads', 'Read_GC', 'Median_fold', 'Std_Dev'], axis=1)
    cntrl_df = cntrl_df.drop(['Avg_fold', 'Plus_reads', 'Minus_reads', 'Read_GC', 'Median_fold', 'Std_Dev'], axis=1)
    # Append results
    target_df = target_df.merge(target_stats_final, left_on='ID', right_on='Target_Gene_Name')
    cntrl_df = cntrl_df.merge(control_stats_final, left_on='ID', right_on='Control_Gene_Name')
    # Filter cntrl dataframe values to drop if 'Avg_fold' is less than < 10
    cntrl_filter_drop = cntrl_df[(cntrl_df['Avg_fold'] < 10) | (cntrl_df['Covered_percent'] < 90.0)]
    cntrl_filter_drop_outfile = '{0}/{1}_cntrl_filter_drop_covstats.tsv'.format(tmp, sample_name)
    cntrl_filter_drop.to_csv(cntrl_filter_drop_outfile, sep="\t", index=False)
    # Filter cntrl dataframe values to keep if 'Avg_fold' is >= 10 or 'Covered_percent' >= 90.0
    cntrl_filter_df = cntrl_df[(cntrl_df['Avg_fold'] >= 10) & (cntrl_df['Covered_percent'] >= 90.0)]
    # Calculate Avg_fold * length and create new column
    cntrl_filter_df['avg_fold_X_length'] = cntrl_filter_df['Avg_fold'] * cntrl_filter_df['Length']
    # Extract average_fold from control dataframe
    average_fold_cntrl = cntrl_filter_df['avg_fold_X_length'].sum() / cntrl_filter_df['Length'].sum()
    # Append average_fold to target dataframe
    target_df["Avg_fold_cntrl"] = average_fold_cntrl
    # Calculate [target_gene_Avg_fold / cntrl_gene_Avg_fold ] if df["Covered_percent"] > 90
    # POSSIBLY MODIFY DUE TO ONT ISSUES
    target_df.loc[target_df["Covered_percent"] >= 90, 'norm_cov_depth'] = target_df["Avg_fold"] / target_df[
        "Avg_fold_cntrl"]
    outfile = '{0}/{1}_convict_results.tsv'.format(outdir, sample_name)
    # Reformat columns for output and remove unnecessary columns
    target_df = target_df[['ID', 'Length', 'Covered_bases', 'Covered_percent', 'Avg_fold', 'Std_dev',
                           'Avg_fold_cntrl', 'norm_cov_depth']]
    target_df.to_csv(outfile, sep="\t", index=False)
    # Potentially remove variant calling step if only gives redundant information from kmerresistance
    # Call variants using bcftools and normalize variants for ambiguous REF/ALT sites
    print("\nVariant calling on consensus target fasta files to get consensus fasta file from short-read alignment\n")
    subprocess.Popen('bcftools mpileup -Ou --threads ' + threads + ' -f ' + KmerRes_fasta_out + ' ' +
                     tmp + '/target_align_sort.bam' + ' | bcftools call -mv -Ou --threads ' + threads + ' --ploidy' +
                     ' ' + args.ploidy + ' | bcftools norm -Oz -f ' + KmerRes_fasta_out + ' --threads ' + threads + ' -o ' +
                     tmp + '/norm_calls.vcf.gz', shell=True).wait()
    subprocess.Popen('bcftools index ' + tmp + '/norm_calls.vcf.gz', shell=True).wait()
    # Need to mask reference.fasta for regions with zero coverage before creating consensus fasta file
    subprocess.Popen('bedtools genomecov -bga -split -ibam ' + tmp + '/target_align_sort.bam' + ' > ' +
                     tmp + '/coverage.bed', shell=True).wait()
    subprocess.Popen('grep -w 0$ ' + tmp + '/coverage.bed' + ' > ' + tmp + '/zero-coverage.bed', shell=True).wait()
    subprocess.Popen('bedtools maskfasta -fi ' + KmerRes_fasta_out + ' -bed ' + tmp + '/zero-coverage.bed' + ' -fo ' +
                     tmp + '/ref_zero_masked.fasta', shell=True).wait()
    # Create consensus fasta file from alignment
    subprocess.Popen('bcftools consensus --fasta-ref ' + tmp + '/ref_zero_masked.fasta' + ' -o ' +
                     outdir + '/consensus_align.fasta' + ' ' + tmp + '/norm_calls.vcf.gz', shell=True).wait()
    blastout = '{0}/blast_results_tmp.tsv'.format(tmp)
    consensus_align_fasta = '{0}/consensus_align.fasta'.format(outdir)
    blastn_cmd = 'blastn -subject {0} -query {1} -outfmt \"6 std qcovs\" -out {2} -evalue 0.00001 -perc_identity 0.9' \
                 ' -qcov_hsp_perc 0.9 -culling_limit 1'.format(KmerRes_fasta_out, consensus_align_fasta, blastout)
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

