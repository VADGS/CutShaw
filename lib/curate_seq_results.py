#!/usr/bin/env python3

# author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os
import sys
import shutil
import argparse
import re
import csv
import pandas
import datetime
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))
from CutShaw.core import fileparser
from CutShaw.core import calldocker
from CutShaw.lib import run_fastani
from CutShaw.lib import run_cfsansnp
from CutShaw.lib import run_quast


class CutShaw:
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

        self.cutshaw_out_dir = self.output_dir + "/CutShaw_output/"

    def seq_results(self):
        # create output directory
        cutshaw_out_dir = self.cutshaw_out_dir
        db_name = os.path.basename(self.db)

        if not os.path.isdir(cutshaw_out_dir):
            os.makedirs(cutshaw_out_dir)
            print("Directory for CutShaw output made: ", cutshaw_out_dir)

        seq_results = {}

        for read in self.runfiles.reads:
            # get id
            id = self.runfiles.reads[read].id
            fastani_result = "/fastani_%s.out"%id
            # change self.path to local dir if path is a basemounted dir
            if os.path.isdir(self.path + "/AppResults"):
                self.path = self.output_dir

            fastani_obj = run_fastani.FastANI(path=self.path, output_dir=self.output_dir)
            fastani_taxon = fastani_obj.fastani()[0][id]
            fastani_reference = fastani_obj.fastani()[1][id].replace(".fasta", "")

            seq_results[id] = {"SampleID": id, "IsolateID": fastani_reference, "Organism": fastani_taxon, "Genus": fastani_taxon.split()[0],
                               "Sequencer": "Illumina MiSeq", "Machine": None, "FlowCell": None, "LibKit": "Nextera XT",
                               "Chemistry": "NA", "RunDate": "NA", "SequencedBy": "DCLS", "SamplesPerRun": 16,
                               "SeqLength": "NA", "Reads": None, "MeanR1Qual": None, "MeanR2Qual": None,
                               "PercMapped": None, "MeanDepth": None, "CovLT10": 0, "SNPs": None, "MeanInsert": None,
                               "NG50": None, "GenomeFraction": None, "Contigs": None, "LengthDelta": None,
                               "UnalignedLength": None, "MostAbundantOrganism": "NA", "Misannotated": None, "Coverage":
                               None}

            # From CG_pipeline grab: number of reads, Mean R1 and R2. and coverage
            if not os.path.isfile("%s/cg_pipeline_output/%s_readMetrics.tsv" % (output_dir, id)):
                CGPipeline_obj = run_cfsansnp.CGPipeline(path=self.path, output_dir=self.output_dir)
                CGPipeline_obj.read_metrics(from_mash=from_mash)

            with open("%s/cg_pipeline_output/%s_readMetrics.tsv" % (output_dir, id)) as tsv_file:
                tsv_reader = list(csv.DictReader(tsv_file, delimiter="\t"))

                for line in tsv_reader:
                    if any(fwd_format in line["File"] for fwd_format in ["_1.fastq", "_R1.fastq"]):
                        seq_results[id]["MeanR1Qual"] = line["avgQuality"]
                        seq_results[id]["Reads"] = int(line["numReads"])
                        seq_results[id]["Coverage"] = float(line["coverage"])
                    if any(rev_format in line["File"] for rev_format in ["_2.fastq", "_R2.fastq"]):
                        seq_results[id]["MeanR2Qual"] = line["avgQuality"]
                        seq_results[id]["r2_totalBases"] = line["totalBases"]
                        seq_results[id]["Reads"] += int(line["numReads"])
                        seq_results[id]["Coverage"] += float(line["coverage"])

            # From Quast grab: GenomeFraction, Contigs, NG50, Unaligned length
            if not os.path.isfile("%s/quast_output/%s/transposed_report.tsv" %(output_dir, id)):
                quast_obj = run_quast.Quast(path=self.path, output_dir=self.output_dir)
                quast_obj.quast()

            with open("%s/quast_output/%s/transposed_report.tsv" %(output_dir, id)) as tsv_file:
                tsv_reader = list(csv.DictReader(tsv_file, delimiter="\t"))

                seq_results[id]["Contigs"] = tsv_reader[0]["# contigs (>= 1000 bp)"]
                seq_results[id]["GenomeFraction"] = int(round(float((tsv_reader[0]["Genome fraction (%)"]))))
                seq_results[id]["NG50"] = int((tsv_reader[0]["NG50"]))
                seq_results[id]["UnalignedLength"] = (tsv_reader[0]["Unaligned length"])
                seq_results[id]["LengthDelta"] = int(abs(int((tsv_reader[0]["Total length"])) -
                                                  int((tsv_reader[0]["Reference length"]))))

            if seq_results[id]["SampleID"] == seq_results[id]["IsolateID"]:
                seq_results[id]["Misannotated"] = "FALSE"
            else:
                seq_results[id]["Misannotated"] = "TRUE"

            # From CFSAN-SNP grab: MeanDepth, MeanInsert, Percent Mapped,  SNPs, Machine, Flowcell
            if not os.path.isfile("%s/cfsansnp_output/%s/metrics.tsv" %(output_dir, id)):
                cfsansnp_obj = run_cfsansnp.CfsanSnp(path=self.path, output_dir=self.output_dir)
                cfsansnp_obj.cfsansnp()

            with open("%s/cfsansnp_output/%s/metrics.tsv" %(output_dir, id)) as tsv_file:
                tsv_reader = list(csv.DictReader(tsv_file, delimiter="\t"))

                seq_results[id]["MeanDepth"] = tsv_reader[0]["Average_Pileup_Depth"]
                seq_results[id]["MeanInsert"] = int(round(float(tsv_reader[0]["Average_Insert_Size"])))
                seq_results[id]["PercMapped"] = tsv_reader[0]["Percent_of_Reads_Mapped"]
                seq_results[id]["SNPs"] = tsv_reader[0]["Phase2_Preserved_SNPs"]
                seq_results[id]["FlowCell"] = tsv_reader[0]["Flowcell"]
                seq_results[id]["Machine"] = tsv_reader[0]["Machine"]

        seq_results_out = "%s/seq_results.tsv"%cutshaw_out_dir

        # Change data dictionary to dataframe to csv
        df = pandas.DataFrame(seq_results).T[["SampleID", "IsolateID", "Organism", "Genus","Sequencer","Machine", "FlowCell", "LibKit",
                                               "Chemistry", "RunDate", "SequencedBy", "SamplesPerRun", "Reads", "SeqLength",
                                               "MeanR1Qual", "MeanR2Qual", "PercMapped", "MeanDepth", "CovLT10", "SNPs",
                                               "MeanInsert", "NG50", "GenomeFraction", "Contigs", "LengthDelta",
                                               "UnalignedLength", "MostAbundantOrganism", "Misannotated", "Coverage"]]
        df.to_csv(seq_results_out, sep = '\t', index=False)

        print("Seq results curated. Output saved as %s"%seq_results_out)


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

    CutShaw_obj = CutShaw(path=path,output_dir=output_dir)
    CutShaw_obj.seq_results()
