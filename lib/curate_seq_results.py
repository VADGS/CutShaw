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


class SeqResults:
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

        self.seq_results_out_dir = self.output_dir + "/seq_results_output/"

    def seq_results(self):
        # create output directory
        seq_results_out_dir = self.seq_results_out_dir
        db_name = os.path.basename(self.db)

        if not os.path.isdir(seq_results_out_dir):
            os.makedirs(seq_results_out_dir)
            print("Directory for seq_results output made: ", seq_results_out_dir)
        seq_results_result = "/%s/seq_result.tsv" % (seq_results_out_dir, id)

        seq_results = {}

        for read in self.runfiles.reads:
            # get id
            id = self.runfiles.reads[read].id
            fastani_result = "/fastani_%s.out"%id
            # change self.path to local dir if path is a basemounted dir
            if os.path.isdir(self.path + "/AppResults"):
                self.path = self.output_dir

            fastani_obj = run_fastani.FastANI(path=path, output_dir=output_dir)
            fastani_taxon = fastani_obj.fastani()[0][id]
            fastani_reference = fastani_obj.fastani()[1][id].replace(".fasta", "")


            seq_results[id] = {"SampleID": id, "IsolateID": fastani_reference, "Genus": fastani_taxon,
                               "Sequencer": "Illumina MiSeq", "Machine": "NA", "FlowCell": "NA", "LibKit": "Nextera XT",
                               "Chemistry": "NA", "RunDate": "NA", "SequencedBy": None, "SamplesPerRun": 16,
                               "SeqLength": "NA", "Reads": None, "MeanR1Qual": None, "MeanR2Qual": None,
                               "PercMapped": None, "MeanDepth": None, "CovLT10": None, "SNPs": None, "MeanInsert": None,
                               "NG50": None, "GenomeFraction": None, "Contigs": None, "LengthDelta": None,
                               "UnalignedLength": None, "MostAbundantOrganism": "NA", "Misannotated": None, "Coverage":
                               None}





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