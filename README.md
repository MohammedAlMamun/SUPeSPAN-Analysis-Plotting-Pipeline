# SUPeSPAN-Analysis-Plotting-Pipeline
Introduction

This repository contains the SUPeSPAN_Analysis R script, a streamlined tool to analyze strand-specific DNA replication data by integrating supernatant sequencing (SUP) with IP-based eSPAN. It runs primary analysis for sequencing data i.e. alignment against reference genome, quality check, binned coverage and peak detection. Furthermore, it does the secondary analysis of averaging the signal around activated replication origin set in a given experiment, calculates watson and crick strand average, lagging and leading compartemnts, nascent and parent fractions around replication origins. The script also simulates similar data for randomly sampled genomic corrdinates and compares it to the replication origin data for noise filtering and allow a comprehensive comparison between different strands and compartments to detect protein binding asymmetries. 

The script was designed to allow lab members to quickly perform analyses without opening or modifying the full workflow. It handles multiple sequencing datasets, combines coverage and peak information, and provides strand-specific output for replication studies.

# Function Overview

SUPeSPAN_Analysis(
  Input_R1, Input_R2,
  BrDU_R1, BrDU_R2,
  ChIP_R1, ChIP_R2,
  eSPAN_R1, eSPAN_R2,
  bSUP_R1, bSUP_R2,
  eSUP_R1, eSUP_R2,
  ExpTitle = "None",
  PeakClass = "brdu",
  Window = "None",
  iterations = 100
)

# Arguments

| Argument      | Description                                                   |
| ------------- | ------------------------------------------------------------- |
| `Input_R1/R2` | Strand-specific input sequencing files (paired)               |
| `BrDU_R1/R2`  | BrdU-labeled newly synthesized DNA (paired)                   |
| `ChIP_R1/R2`  | ChIP IP data (paired)                                         |
| `eSPAN_R1/R2` | Standard eSPAN sequencing (paired)                            |
| `bSUP_R1/R2`  | BrdU supernatant sequencing (paired)                          |
| `eSUP_R1/R2`  | eSPAN supernatant sequencing (paired)                         |
| `ExpTitle`    | Title for the experiment (used in outputs/plots)              |
| `PeakClass`   | Peak type to use: `"brdu"` (default) or `"early"`  or `"late"`|
| `Window`      | Optional genomic window (default: 95% peak-spread)            |
| `iterations`  | Number of iterations for simulation (default: 100)            |

# Quick Start Example

# save the file named SUPeSPAN_Analysis.R in a suitable location in your work station

# open Run_SUPeSPAN.R in Rstudio and edit as needed or open a new R script and then

# Source the analysis script
source("path/to/SUPeSPAN_Analysis.R")

# Run the analysis with example files
          SUPeSPAN_Analysis( Input_R1 = "/path/to/Input_R1_001.fastq.gz",
                             Input_R2 = "/path/to/Input_R2_001.fastq.gz",
                             BrDU_R1 = "/path/to/BrDU_R1_001.fastq.gz",
                             BrDU_R2 = "/path/to/BrDU_R2_001.fastq.gz",
                             ChIP_R1 = "/path/to/ChIP_R1_001.fastq.gz",
                             ChIP_R2 = "/path/to/ChIP_R2_001.fastq.gz",
                             eSPAN_R1 = "/path/to/eSPAN_R1_001.fastq.gz",
                             eSPAN_R2 = "/path/to/eSPAN_R2_001.fastq.gz",
                             bSUP_R1 = "/path/to/bSUP_R1_001.fastq.gz",
                             bSUP_R2 = "/path/to/bSUP_R2_001.fastq.gz",
                             eSUP_R1 = "/path/to/eSUP_R1_001.fastq.gz",
                             eSUP_R2 = "/path/to/eSUP_R2_001.fastq.gz",
                             
                             ExpTitle = "Your-Sample-Title" )

# Note: This script outputs strand-specific coverage, peak analyses, and integrated SUP + eSPAN results and makes downstream plotting and statistical analysis.
