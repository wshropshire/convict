#! /usr/bin/env python3
import argparse
import os
import glob
import pandas as pd

def create_directory(outdir):
    """Function to create the output directory"""
    if os.path.isdir(outdir):
        raise Exception("Directory already exists")
    if not os.path.isdir(outdir):
        os.system("mkdir -p %s" % outdir)
    return

def get_arguments():
    """"Parser assembly arguments"""
    parser = argparse.ArgumentParser(description="Script to combine multiple convict result output results")
    # Help Arguments
    help_group = parser.add_argument_group("Help")
#    help_group.add_argument('-h', action='help', help='Show this help message and exit')
    help_group.add_argument('-V', '--version', action='version', version='%(prog)s version 0.1',
                            help="Show convict's version number")
    # input_arguments
    input_group = parser.add_argument_group("Inputs")
    input_group.add_argument('-f', '--output_files', required=True, metavar='IN FILE',
                             type=argparse.FileType('r'), nargs='*', help='tab delimited output files to merge')
    input_group.add_argument('-p', '--prefix', required=False, type=str,
                             help="Prefix to append as suffix to output files: default=SAMPLE", default='SAMPLE')
    input_group.add_argument('-c', '--column', required=False, type=str,
                             help='Column to select for detailed results by sample: default=norm_cov_depth', 
                             default='norm_cov_depth')

    # output_arguments
    output_group = parser.add_argument_group("Outputs")
    output_group.add_argument('-o', '--outdir', required=False, help="Name of the output directory", type=str)

    args = parser.parse_args()
    return args

def run_conditions():
    args = get_arguments()
    df_list = args.output_files
    prefix = args.prefix
    if args.outdir:
        create_directory(args.outdir)
    else:
        args.outdir = os.getcwd()
    outdir = args.outdir
    data = []
    for tsv in df_list:
        frame = pd.read_csv(tsv, sep='\t', header=0)
        frame['filename'] = tsv.name
        data.append(frame)
    merge_df = pd.concat(data, ignore_index=True)
    pivot_cnv = merge_df.pivot(index='filename', columns='ID', values=args.column)
    results = '{0}/{1}_concatenated_results.tsv'.format(outdir, prefix)
    results1 = '{0}/{1}_{2}_results.tsv'.format(outdir, prefix, args.column)
    merge_df.to_csv(results, sep='\t')
    pivot_cnv.to_csv(results1, sep='\t')


if __name__ == '__main__': run_conditions()