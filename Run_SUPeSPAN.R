
          # Source the analysis script
          source("path/to/SUPeSPAN_Analysis.R")

          # run the analysis
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