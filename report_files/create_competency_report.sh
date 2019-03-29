#!/bin/bash

#----------
# Defaults 
#----------

kit='miseqv2-250'
log_file=''
output_dir='.'
prefix='Competency_report'
quantitative='FALSE'
report_template=''
template='competency_2018_template.Rnw'
title=''

# Input Files
run_stats='run_stats.tsv'
run_results='run_results.tsv'
seq_stats='seq_stats.tsv'
seq_results='seq_results.tsv'

#-------
# Usage 
#-------

display_usage() { 
  echo -e "\nUsage: $0 -o 'output_dir' -p 'prefix' -t 'report_title'\n"
  echo -e "  -k  type of sequencing kit."
  echo -e "      Default: miseqv2-250" 
  echo -e "      Other options: miseqv3-300, miseqv2micro-150, miseqv2nano-250,"
  echo -e "                     miniseqMid-150, miniseqHigh-150" 
  echo -e "  -l  Log file."
  echo -e "      Default: <prefix>.log" 
  echo -e "  -m  message for inclusion in the results section."
  echo -e "      Default: ''" 
  echo -e "  -o  output directory."
  echo -e "      Default: ." 
  echo -e "  -s  Sequencing results file."
  echo -e "      Default: seq_results.tsv" 
  echo -e "  -S  Sequencing stats file."
  echo -e "      Default: seq_stats.tsv" 
  echo -e "  -t  Report title."
  echo -e "      Default: ''" 
  echo -e "  -T  Report template (.Rnw)."
  echo -e "      Default: competency_2018_template.Rnw" 
  echo -e "  -h  Displays this help message\n"
  echo -e "create_competency_report.sh"
  echo -e "This script is part of a competency assessment analysis pipeline and produces "
  echo -e "quantitative and graphical reports.\n"
  exit
} 


#---------
# Options
#---------

while getopts ":k:l:m:o:p:s:S:t:T:h" opt; do
  case $opt in
    k) kit=$OPTARG;;
    l) log_file=$OPTARG;;
    m) message=$OPTARG;;
    o) output_dir=${OPTARG%/};;
    p) prefix=$OPTARG;;
    s) seq_results=$OPTARG;;
    S) seq_stats=$OPTARG;;
    t) title=$OPTARG;;
    T) report_template=$OPTARG;;
    h) display_usage;;  
   \?) #unrecognized option - show help
      echo -e \\n"Option -$OPTARG not allowed."
      display_usage;;
  esac
done

if [ -z $log_file ]; then 
  log_file="${output_dir}/create_competency_report.log"
fi  


#----------------------
# Create report
#----------------------

sweave_report="${output_dir}/${prefix}.Rnw"
tex_report="${output_dir}/${prefix}.tex"
cp $report_template $sweave_report
  
#module load R

echo "Converting Sweave to Tex" | tee $log_file

cd ${output_dir}

R -e "Sweave('$sweave_report')" --args  --kit="$kit" --message="$message" \
  --seq="$seq_results" \
  --seqstats=$"$seq_stats" >>$log_file 2>&1
exit_status=$?    

if [ "$exit_status" -ne 0 ]; then
  echo "Sweave was unable to complete the conversion. Please check the log file." \
    | tee -a $log_file
  exit 1
fi  

#----------------------------------------------------
# Insert title into tex report
#----------------------------------------------------
echo $title
title_clean=$(echo $title | sed 's/"//g')
echo $title_clean
sed -i['.bak'] "s/\"Title goes here\"/ ${title_clean}/" $tex_report || true

echo "Tex conversion successful." | tee -a $log_file

#---------------
# Create report
#---------------

echo "Converting Tex to Pdf" | tee -a $log_file
#module load latex
# Convert it twice to avoid problems
R CMD pdflatex -interaction=nonstopmode $tex_report 2>&1 >/dev/null
R CMD pdflatex -interaction=nonstopmode $tex_report >>$log_file 2>&1
exit_status=$?

if [ "$exit_status" -ne 0 ]; then
  echo "The tex to pdf conversion had a non-zero exit status. Please check the pdf file" \
  " and the log file." | tee -a $log_file
  exit 1
fi  

echo "Report creation successful." | tee -a $log_file
