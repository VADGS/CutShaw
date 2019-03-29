#!/usr/bin/env python3

#author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os
import sys
import shutil
import csv
import argparse
import subprocess
from shutil import copyfile
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

from CutShaw.core import fileparser
from CutShaw.lib import curate_seq_results



def main():
    parser = argparse.ArgumentParser(usage="cutshaw.py <input> [options]")
    parser.add_argument("input", type=str, nargs='?', help="path to dir containing read files")
    parser.add_argument("-o", default="", nargs='?', type=str, help="Name of output_dir")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    path=args.input
    output_dir = args.o
    analyst= input("Enter first and last name of WGS analyst: ")
    prefix = analyst.replace(" ", "_")
    prefix = f"{prefix}_CutShaw_report"

    if output_dir.endswith('/'):
        output_dir = output_dir[:-1]

    if '/' in output_dir and "CutShaw_output" not in output_dir:
        project = output_dir.split('/')
        project = project[-1]
    elif "CutShaw_output" in output_dir:
        project = datetime.datetime.today().strftime('%Y-%m-%d')
    else:
        project = output_dir

    create_competency_report = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))[:-4]\
                               + "/report_files/create_competency_report.sh"
    report_template = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))[:-4]\
                      + "/report_files/report_template.Rnw"
    seq_stats = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))[:-4]\
                + "/report_files/seq_stats.tsv"
    seq_results = os.path.abspath(output_dir) + "/CutShaw_output/seq_results.tsv"


    # Curate seq_results file
    print("Initiating quality assessment of read data")
    CutShaw_obj = curate_seq_results.CutShaw(path=path,output_dir=output_dir)
    CutShaw_obj.seq_results()

    # Generate CutShaw report

    command = [create_competency_report, "-p", prefix, "-t", f"\"{analyst}\"", "-T", report_template,
               "-o", os.path.abspath(output_dir)+ "/CutShaw_output", "-s", seq_results, "-S", seq_stats]

    reports_dir = output_dir + "/reports/"

    CutShaw_report = f"{prefix}.pdf"

    if not os.path.isdir(reports_dir):
        os.makedirs(reports_dir)
        print("Directory for WGS reports made:", reports_dir)

    if not os.path.isfile(f"{reports_dir}/{CutShaw_report}"):
        #print(str(command))
        subprocess.call(command)
        copyfile(f"{output_dir}/CutShaw_output/{CutShaw_report}", f"{reports_dir}/{CutShaw_report}")

        print(f"CutShaw is complete! Output saved as {reports_dir}{CutShaw_report}")
    else:
        print(f"CutShaw Report exists in {reports_dir}")

if __name__ == '__main__':
    main()
