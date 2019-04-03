#!/usr/bin/env python
import os, sys, json
import argparse
from fastq_utils import run_fq_util
if __name__ == "__main__":
    #modelling input parameters after rockhopper
    parser = argparse.ArgumentParser()
    #if you want to support multiple genomes for alignment you should make this json payload an nargs+ parameter
    parser.add_argument('--jfile',
            help='json file for job {"reference_genome_id": "1310806.3", \
                    "output_file": "rnaseq_baumanii_1505311", \
                    "recipe": ["FASTQC","TRIM","ALIGN"], "output_path": "/anwarren@patricbrc.org/home/test",\
                    "paired_end_libs": [{"read1": "/anwarren@patricbrc.org/home/rnaseq_test/MHB_R1.fq.gz",\
                    "read2": "/anwarren@patricbrc.org/home/rnaseq_test/MHB_R2.fq.gz"},\
                    {"read1": "/anwarren@patricbrc.org/home/rnaseq_test/MERO_75_R1.fq.gz",\
                    "read2": "/anwarren@patricbrc.org/home/rnaseq_test/MERO_75_R2.fq.gz"}]}', required=True)
    parser.add_argument('--sstring', help='json server string specifying api {"data_api":"url"}', required=True, default=None)
    parser.add_argument('--index', help='flag for enabling using HISAT2 indices', action='store_true', required=False)
    #parser.add_argument('-L', help='csv list of library names for comparison', required=False)
    #parser.add_argument('-C', help='csv list of comparisons. comparisons are library names separated by percent. ', required=False)
    parser.add_argument('-p', help='JSON formatted parameter list for info keyed to program', default=None, required=False)
    parser.add_argument('-o', help='output directory. defaults to current directory.', required=False, default=None)
    if len(sys.argv) ==1:
        parser.print_help()
        sys.exit(2)
    map_args = parser.parse_args()
    assert map_args.d.startswith(".") # job object folder name needs a .
        
    #create library dict
    with open(map_args.jfile, 'r') as job_handle:
        job_data = json.load(job_handle)
    server_info = json.loads(map_args.sstring)
    for k,d in server_info.iteritems():
        job_data[k]=d
    if map_args.o == None:
        output_dir="./"
    else:
        output_dir=map_args.o
    count=0
    if map_args.p != None and os.path.exists(map_args.p):
    	tool_params=json.load(open(map_args.p,'r'))
    else:
        tool_params={}

    run_fq_util(job_data,output_dir,tool_params)
