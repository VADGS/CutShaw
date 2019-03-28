#!/usr/bin/env python3

# author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os
import sys
import shutil
import csv
import argparse
import re
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))
from CutShaw.core import fileparser
from CutShaw.core import calldocker
from CutShaw.lib import run_spades

class FastANI:
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

        self.fastani_out_dir = self.output_dir + "/fastani_output/"

    def fastani(self):
        # create output directory
        fastani_out_dir = self.fastani_out_dir

        if not os.path.isdir(fastani_out_dir):
            os.makedirs(fastani_out_dir)
            print("Directory for fastani output made: ", fastani_out_dir)

        taxons = {}
        reference_genomes = {}

        for read in self.runfiles.reads:
            # get id
            id = self.runfiles.reads[read].id
            fastani_result = "/fastani_%s.out"%id
            # change self.path to local dir if path is a basemounted dir
            if os.path.isdir(self.path + "/AppResults"):
                self.path = self.output_dir
            assembly = "/spades_output/%s/contigs.fasta" % id

            if not os.path.isfile(assembly):
                spades_obj = run_spades.Spades(path=self.path, output_dir=self.output_dir)
                spades_obj.spades()

            # create paths for data
            mounting = {self.path:'/datain',fastani_out_dir:'/dataout', self.db: '/db'}
            ref_list = "/reference_list.txt"
            out_dir = '/dataout/'
            in_dir = '/datain/'
            db = '/db/'

            command = "bash -c 'fastANI -q {in_dir}{assembly} --rl {db}{ref_list} -o " \
                      "{out_dir}/{fastani_result}'".format(assembly=assembly,out_dir=out_dir,db=db, in_dir=in_dir,
                                                           ref_list=ref_list, fastani_result=fastani_result)

            # call the docker process
            if not os.path.isfile("%s/%s"%(fastani_out_dir, fastani_result)):
                print("Generating FastANI report for sample " + id)
                calldocker.call("staphb/fastani",command,'/dataout',mounting)

            with open("%s/%s"%(fastani_out_dir, fastani_result)) as file:
                tsv_reader = csv.reader(file, delimiter="\t", quotechar='"')
                predicted_taxon = ""
                reference_genome = str(os.path.basename(next(tsv_reader)[1]))
                if "SAP18-0432" in reference_genome:
                    predicted_taxon = "Salmonella enterica subsp. enterica serover Enteritidis"
                elif "SAP18-H9654" in reference_genome:
                    predicted_taxon = "Salmonella enterica subsp. enterica serover Enteritidis"
                elif "SAP18-6199" in reference_genome:
                    predicted_taxon = "Salmonella enterica subsp. enterica serover Typhimurium"
                elif "SAP18-8729" in reference_genome:
                    predicted_taxon = "Salmonella enterica subsp. enterica serover Newport"
                elif "LMP18-H2446" in reference_genome:
                    predicted_taxon = "Listeria monocytogenes"
                elif "LMP18-H8393" in reference_genome:
                    predicted_taxon = "Listeria monocytogenes"
                else:
                    raise ValueError("Sample %s not identified as a 2018 PT isolate"%id)

            taxons[id] = predicted_taxon
            reference_genomes[id] = reference_genome
            print("Isolate: %s " %id)
            print("reference genome: %s"%reference_genome)

           # print("Sample: " + id + " FastANI Hit: " + str(predicted_taxon))

        return [taxons, reference_genomes]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="run_fastani.py <input> [options]")
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

    fastani_obj = FastANI(path=path,output_dir=output_dir)
    fastani_obj.fastani()




