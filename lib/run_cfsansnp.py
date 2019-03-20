#!/usr/bin/env python3

# author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os
import sys
import shutil
import argparse
import re
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))
from CutShaw.core import fileparser
from CutShaw.core import calldocker
from CutShaw.lib import run_fastani


class CfsanSnp:
    # class object to contain fastq file information
    runfiles = None
    # path to fastq files
    path = None
    # output directory
    output_dir = None

    def __init__(self, runfiles=None, path=None, output_dir = None):
        reference_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))[:-4] + "/PT_genomes"
        if output_dir:
            self.output_dir = os.path.abspath(output_dir)
        else:
            self.output_dir = os.getcwd()

        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)

        if runfiles:
            self.runfiles = runfiles
        else:
            self.path = path
            self.runfiles = fileparser.RunFiles(self.path, output_dir=output_dir)

        self.db=reference_dir

        self.cfsansnp_out_dir = self.output_dir + "/cfsansnp_output/"

    def cfsansnp(self):
        # create output directory
        cfsansnp_out_dir = self.cfsansnp_out_dir
        db_name = os.path.basename(self.db)
        cfsan_read_dir = cfsansnp_out_dir + "/cfsan-reads/"

        if not os.path.isdir(cfsansnp_out_dir):
            os.makedirs(cfsansnp_out_dir)
            print("Directory for cfsansnp output made: ", cfsansnp_out_dir)
        if not os.path.isdir(cfsansnp_out_dir):
            os.makedirs(cfsansnp_out_dir)
            print("Directory for cfsansnp read dir made: ", cfsan_read_dir)

        for read in self.runfiles.reads:
            # get id
            id = self.runfiles.reads[read].id
            cfsansnp_result = "/%s/%s/snpma.fasta"%(cfsansnp_out_dir, id)
            # change self.path to local dir if path is a basemounted dir
            if os.path.isdir(self.path + "/AppResults"):
                self.path = self.output_dir

            fastani_obj = run_fastani.FastANI(path=path, output_dir=output_dir)
            fastani_reference = fastani_obj.fastani()[1][id]

            # create paths for data
            mounting = {self.path: '/datain', cfsansnp_out_dir: '/dataout', self.db: '/db'}
            out_dir = '/dataout/'
            in_dir = '/datain/'
            db = '/db/'

            fwd_read = "/%s/raw_reads/"%in_dir + os.path.basename(self.runfiles.reads[read].fwd)
            rev_read = "/%s/raw_reads/"%in_dir + os.path.basename(self.runfiles.reads[read].rev)

            if not os.path.isdir(cfsan_read_dir + id):
                os.makedirs(cfsan_read_dir + id)
            if not os.path.islink(cfsan_read_dir + id + "/" + os.path.basename(self.runfiles.reads[read].fwd)):
                os.symlink(fwd_read, cfsan_read_dir + id + "/" + os.path.basename(self.runfiles.reads[read].fwd))
            if not os.path.islink(cfsan_read_dir + id + "/" + os.path.basename(self.runfiles.reads[read].rev)):
                os.symlink(rev_read, cfsan_read_dir + id + "/" + os.path.basename(self.runfiles.reads[read].rev))

            reference_genome = "/%s"%fastani_reference
            command = "bash -c 'run_snp_pipeline.sh -m soft -o {out_dir}{id} -s {out_dir}cfsan-reads " \
                      "{db}{reference_genome}'".format(out_dir=out_dir,db=db, in_dir=in_dir, id=id,
                                                       reference_genome=reference_genome)
            # call the docker process
            if not os.path.isfile(cfsansnp_result):
                print("Generating cfsansnp output for sample " + id)
                 calldocker.call("staphb/cfsan-snp-pipeline:2.0.2",command,'/dataout',mounting)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="run_cfsansnp.py <input> [options]")
    parser.add_argument("input", type=str, nargs='?', help="path to dir containing read files")
    parser.add_argument("-o", default="", nargs='?', type=str, help="Name of output_dir")


    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    path = os.path.abspath(args.input)
    output_dir = args.o

    if not output_dir:
        output_dir = os.getcwd()

    cfsansnp_obj = CfsanSnp(path=path,output_dir=output_dir)
    cfsansnp_obj.cfsansnp()
