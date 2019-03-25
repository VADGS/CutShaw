wgs_competency
================

The current scripts provide functionality for creation of quantitative and graphical reports (see example reports) summarizing a bacterial whole genome DNA sequencing competency assessment analysis pipeline (more scripts to come). The included templates may serve as the basis for customization. The files seq_stats.tsv and run_stats.tsv are comprised of summary stats derived from a distribution of metrics for six isolates sequenced during the 2018 PulseNet-GenomeTrakr harmonized Proficiency Testing (PT) exercise as part of a 16-sample 2X250 Illumina MiSeq run (V2 chemistry) using the Nextera XT DNA Library Prep Kit. For valid use in competency assessments, equivalent SOPs should be employed for sequencing and analysis of the isolates listed below. 

## Software Dependencies

#### R:
    3.4 or higher
    
#### R libraries:
    cowplot
    ggplot2
    ggsci
    grid
    gridExtra
    optparse
    xtable

## Usage

In order to create a graphical report - 

    create_competency_report.sh -p <prefix> -t <Report title> -T graphical_template.Rnw

In order to create a quantitative report - 

    create_competency_report.sh -q -p <prefix> -t <Report title> -T quantitative_template.Rnw

## Bacterial Isolates

    SAP18-0432    Salmonella enterica subsp. enterica serovar Enteritidis
    SAP18-H9654   Salmonella enterica subsp. enterica serovar Enteritidis
    SAP18-6199    Salmonella enterica subsp. enterica serovar Typhimurium
    SAP18-8729    Salmonella enterica subsp. enterica serovar Newport
    LMP18-H2446   Listeria monocytogenes
    LMP18-H8393   Listeria monocytogenes

## Analysis Software
  If the provided summary stats are used, it is recommended to perform analysis using equivalent versions of the following software programs. 

    CFSAN SNP Pipeline v1.0
        bowtie2-2.3.4.1-beta
        picard-2.9.2
        samtools-1.3.1
        VarScan v2.3
    FastQC v0.11.5
    SPAdes Genome Assembler v3.11.1
    QUAST v4.5

## Authors

Joseph D. Baugher, Ph.D.

## Maintainer

Joseph D. Baugher, Ph.D., joseph.baugher@fda.hhs.gov
