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
        db_name = os.path.basename(self.db)

        if not os.path.isdir(fastani_out_dir):
            os.makedirs(fastani_out_dir)
            print("Directory for fastani output made: ", fastani_out_dir)

        for read in self.runfiles.reads:
            # get id
            id = self.runfiles.reads[read].id
            fastani_result = "/fastani_%s.out"%id
            # change self.path to local dir if path is a basemounted dir
            if os.path.isdir(self.path + "/AppResults"):
                self.path = self.output_dir
            assembly = "/spades_output/%s/contigs.fasta" % id

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
                #rd = csv.reader(file, delimiter="\t", quotechar='"')
                for line in file:
                    parts = line.split()
                    predicted_taxon = str(parts[1:2])
                    if "SAP18" in predicted_taxon:
                        predicted_taxon = "Salmonella"
                    elif "LMP18" in predicted_taxon:
                        predicted_taxon = "Listeria"
            print("Sample: " + id + " FastANI Hit: " + str(predicted_taxon))





if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')
    parser = argparse.ArgumentParser(usage="run_fastani.py <input> [options]")
    parser.add_argument("input", type=str, nargs='?', help="path to dir containing read files")
    parser.add_argument("-o", default="", nargs='?', type=str, help="Name of output_dir")
    parser.add_argument("-species", nargs='?', type=str2bool, default=True, help="return dictionary of predicted "
                                                                                 "species (i.e. genus and species of "
                                                                                 "top fastani hits). default: "
                                                                                 "-species=True")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    path = os.path.abspath(args.input)
    output_dir = args.o
    species = args.species

    if not output_dir:
        output_dir = os.getcwd()

    fastani_obj = FastANI(path=path,output_dir=output_dir)
    fastani_obj.fastani()

    # if species:
    #     fastani_obj.fastani_species()


