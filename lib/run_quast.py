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
from CutShaw.lib import run_fastani


class Quast:
    # class object to contain fastq file information
    runfiles = None
    # path to fastq files
    path = None
    # output directory
    output_dir = None

    def __init__(self, runfiles=None, path=None, output_dir=None):
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

        self.db = reference_dir

        self.quast_out_dir = self.output_dir + "/quast_output/"

    def quast(self):
        # create output directory
        quast_out_dir = self.quast_out_dir

        if not os.path.isdir(quast_out_dir):
            os.makedirs(quast_out_dir)
            print("Directory for Quast output made: ", quast_out_dir)

        for read in self.runfiles.reads:
            # get id
            id = self.runfiles.reads[read].id
            # change self.path to local dir if path is a basemounted dir
            if os.path.isdir(self.path + "/AppResults"):
                self.path = self.output_dir
            print(self.path)
            fastani_obj = run_fastani.FastANI(path=path, output_dir=output_dir)
            fastani_reference = fastani_obj.fastani()[1][id]
            reference_genome = "/%s" % fastani_reference

            assembly = "/spades_output/%s/contigs.fasta" % id
            quast_results = "%s/%s/"%(quast_out_dir, id)
            if not os.path.isdir(quast_results):
                os.makedirs(quast_results)

            # create paths for data
            mounting = {self.path:'/datain', quast_results:'/dataout', self.db: '/db'}
            out_dir = '/dataout/'
            in_dir = '/datain/'
            db = '/db/'

            command = "bash -c 'quast.py {in_dir}{assembly} -r {db}{reference_genome} -o {out_dir}'".format(
                assembly=assembly, id=id, out_dir=out_dir, reference_genome=reference_genome, in_dir=in_dir,
                db=db)


            # call the docker process
            #print("Generating Quast assembly metrics")
            calldocker.call("staphb/quast:5.0.0",command,'/dataout/',mounting)

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

    quast_obj = Quast(path=path,output_dir=output_dir)
    quast_obj.quast()
