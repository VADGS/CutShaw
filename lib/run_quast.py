#!/usr/bin/env python3

# author: Kevin Libuit
# email: kevin.libuit@dgs.virginia.gov

import os
import sys
import argparse
import glob
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))
from CutShaw.core import fileparser
from CutShaw.core import calldocker


class Quast:
    # class object to contain fastq file information
    runfiles = None
    # path to fastq files
    path = None
    # output directory
    output_dir = None

    def __init__(self, threads=None, runfiles=None, path=None, output_dir = None, extra_params=None):
        if output_dir:
            self.output_dir = os.path.abspath(output_dir)
        else:
            self.output_dir = os.getcwd()

        if threads:
            self.threads = threads
        else:
            self.threads = 1

        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)

        if runfiles:
            self.runfiles = runfiles
        else:
            self.path = path
            self.runfiles = fileparser.RunFiles(self.path, output_dir=output_dir)

        if extra_params:
            self.extra_params = " ".join(extra_params)
        else:
            self.extra_params = ""

        self.quast_out_dir = self.output_dir + "/quast_output/"

    def quast(self):
        # create output directory
        quast_out_dir = self.quast_out_dir

        assembly_files = []

        for read in self.runfiles.reads:
            # get id
            id = self.runfiles.reads[read].id
            assembly_files.append("/data/spades_output/%s/contigs.fasta" % (id))
        # change self.path to local dir if path is a basemounted dir
        if os.path.isdir(self.path + "/AppResults"):
            self.path = self.output_dir


        assembly_files = ' '.join(assembly_files)
        # create paths for data
        mounting = {self.path: '/data'}
        out_dir = '/data'
        in_dir = '/data'

        command = "bash -c 'quast.py {assembly_files} -o {out_dir}'".format(assembly_files=assembly_files,
                                                                            out_dir=out_dir + "/quast_output/")


        # call the docker process
        #print("Generating Quast assembly metrics")
        calldocker.call("staphb/quast",command,'/dataout',mounting)

        #print("Quast assembly metrics saved to: %s"%(quast_out_dir))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="run_quast.py <input> [options]")
    parser.add_argument("input", type=str, nargs='?', help="path to dir containing read files")
    parser.add_argument("-o", default="", nargs='?', type=str, help="Name of output_dir")
    parser.add_argument("-t",default=16,type=int,help="number of threads")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    path = os.path.abspath(args.input)
    output_dir = args.o
    threads = args.t

    if not output_dir:
        output_dir = os.getcwd()

    quast_obj = Quast(path=path,threads=threads,output_dir=output_dir)
    quast_obj.quast()