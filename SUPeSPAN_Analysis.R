


SUPeSPAN_Analysis <- function( Input_R1, Input_R2, BrDU_R1, BrDU_R2, ChIP_R1, ChIP_R2,
                               eSPAN_R1, eSPAN_R2, bSUP_R1, bSUP_R2, eSUP_R1, eSUP_R2,
                               ExpTitle = "None", PeakClass = "brdu", 
                               Window = "None", iterations = 100 ){
  
  
  ## load packages
  
  packages <- c("basicPlotteR", "plyr", "tidyverse", "dplyr", "plotrix", "rasterpdf", "imager",
                "VennDiagram", "grid", "gridBase", "gridExtra", "ShortRead", "csaw", "shiny",
                "BSgenome.Scerevisiae.UCSC.sacCer3", "Rsubread", "GenomicAlignments", "ORFik",
                "IRanges", "MACSr", "readxl", "data.table")
  
  suppressWarnings(suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE)))
  
  ## load basic files 
  
  All_Ori <- read.table("/Applications/ngsAnalyser.app/Contents/Resources/app/OriginList_Full.bed", header = TRUE, quote = "/t")
  All_Ori_Link <- "/Applications/ngsAnalyser.app/Contents/Resources/app/OriginList_Full.bed"
  
  E_Ori <- read.table("/Applications/ngsAnalyser.app/Contents/Resources/app/E_Rep.bed", header = F, quote = "/t")[ ,1:4]
  L_Ori <- read.table("/Applications/ngsAnalyser.app/Contents/Resources/app/L_Rep.bed", header = F, quote = "/t")[ ,1:4]
  
  image <- load.image("/Applications/ngsAnalyser.app/Contents/Resources/app/ForkPic_1.jpeg")
  
  colnames(E_Ori) <- c("chrom", "chromStart", "chromEnd", "name"); E_Ori$mid <- round((E_Ori$chromStart + E_Ori$chromEnd)/2)
  colnames(L_Ori) <- c("chrom", "chromStart", "chromEnd", "name"); L_Ori$mid <- round((L_Ori$chromStart + L_Ori$chromEnd)/2)
  colnames(All_Ori) <- c("chrom", "chromStart", "chromEnd", "name"); All_Ori$mid <- round((All_Ori$chromStart + All_Ori$chromEnd)/2)
  
  # Sequencing Alignment & Binned Coverage Calculation
  
  Alignment_Mapping_CoverageRatio <- function(Input_R1, Input_R2, BrDU_R1, BrDU_R2, ChIP_R1, ChIP_R2,
                                              eSPAN_R1, eSPAN_R2, bSUP_R1, bSUP_R2, eSUP_R1, eSUP_R2,
                                              slidingWindow = "YES") {
    
    
    
    #
    useDef <- function(a,d) ifelse(isTruthy(a), a,d)
    
    ExpTitle = useDef(ExpTitle, "None")
    
    
    if(ExpTitle == "None"){
      Pro_1 <- unlist(strsplit(basename(Input_R1), split='_', fixed=TRUE))[[1]]
    } else {
      Pro_1 <- ExpTitle
    }
    
    #
    message(paste0("Experiment: ", Pro_1))
    #
    
    suppressWarnings(dir.create(paste0("~/Desktop/", Pro_1)))   #create directory named with the protein in the Desktop
    
    ## Quality check of fastqs'
    #
    message("Running QC ...")
    #
    if(!dir.exists(paste0("~/Desktop/", Pro_1, "/", Pro_1, "_", "QR", ".html"))){
      
      fls = c(Input_R1, Input_R2, ChIP_R1, ChIP_R2, BrDU_R1, BrDU_R2, 
              bSUP_R1, bSUP_R2, eSPAN_R1, eSPAN_R2, eSUP_R1, eSUP_R2)
      
      names(fls) = sub(".fastq", "", basename(fls))
      
      qas = lapply(seq_along(fls),
                   function(i, fls) qa(readFastq(fls[i]), names(fls)[i]),
                   fls)
      qa = do.call(rbind, qas)
      rpt = report(qa, dest = paste0("~/Desktop/", Pro_1, "/", Pro_1, "_", "QR", ".html"))
      
    }
    
    ## Run alignment
    #
    message("Running alignments ...")
    message("Reference yeast genome : S288C")
    #
    if(!dir.exists(paste0("~/Desktop/", Pro_1, "/", "Bam"))){
      
      RunAlign <- function(File_R1, File_R2, SampName){
        
        tempdir(check = TRUE)
        
        Sam <- tempfile(fileext = ".sam")
        Bam <- tempfile(fileext = ".bam")
        nmCollate <- tempfile(fileext = ".bam")
        fixMat <- tempfile(fileext = ".bam")
        SrtBam <- tempfile(fileext = ".bam")
        
        
        Pro_1 <- Pro_1
        Pro_2 <- SampName
        
        suppressWarnings(dir.create(paste0("~/Desktop/", Pro_1, "/", "Bam")))
        
        AlnLog <- paste0("~/Desktop/", Pro_1, "/", "Bam", "/", Pro_1, "_", Pro_2, ".log")
        SFBam <- paste0("~/Desktop/", Pro_1, "/", "Bam", "/", Pro_1, "_", Pro_2, ".bam")
        
        #read the indexed reference genome for the alignment of sequenced data
        ref_index <- "/Applications/ngsAnalyser.app/Contents/Resources/app/bowtie2-2.4.4-macos-x86_64/indexes/S288C_Ref"
        
        #following commands will run the alignemnt, check quality, sort, filter and index the resultant bam file 
        
        system(sprintf("(/Applications/ngsAnalyser.app/Contents/Resources/app/bowtie2-2.4.4-macos-x86_64/bowtie2 -p 8  --no-discordant --fr -x %s -1 %s -2 %s -S %s) 2> %s", 
                       ref_index, File_R1, File_R2, Sam, AlnLog))
        
        system(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools view -bS -@ 15 -q 30 -f 2 %s > %s", Sam, Bam))
        
        system(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools collate -@ 15 -o %s %s", nmCollate, Bam))
        
        system(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools fixmate -@ 15 -m %s %s", nmCollate, fixMat))
        
        system(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools sort -l 9 -@ 15 -m 1024M  -O bam -o %s %s", SrtBam, fixMat))
        
        system(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools markdup -@ 15 %s %s", SrtBam, SFBam))
        
        system(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools index -@ 15 %s", SFBam))
        
        unlink(c(Sam, Bam, nmCollate, fixMat, SrtBam), recursive = T, force = T)
        
      }
      
      RunAlign(Input_R1, Input_R2, "Input")
      RunAlign(ChIP_R1, ChIP_R2, "ChIP")
      RunAlign(BrDU_R1, BrDU_R2, "BrDU")
      RunAlign(bSUP_R1, bSUP_R2, "bSUP")
      RunAlign(eSPAN_R1, eSPAN_R2, "eSPAN")
      RunAlign(eSUP_R1, eSUP_R2, "eSUP")
      
    }
    #
    ## Calculate genome-wide binned coverage
    #
    message("Calculating read coverage ...")
    #
    bamFiles <- Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Bam", "/", "*", ".bam"))
    
     if(!dir.exists(paste0("~/Desktop/", Pro_1, "/", "Coverage"))){
      
      #coverage
      BamCoverage <- function(bamFile, binSize = 300, stepSize = 10, slidingWindow = "YES", byReads_5p = T){
        
        Pro_1 <- unlist(strsplit(tools::file_path_sans_ext(basename(bamFile)), split='_', fixed=TRUE))[[1]] #extract protein name
        Pro_2 <- unlist(strsplit(tools::file_path_sans_ext(basename(bamFile)), split='_', fixed=TRUE))[[2]] #extract sample name
        
        suppressWarnings(dir.create(paste0("~/Desktop/", Pro_1, "/", "Coverage"))) 
        
        
        tempdir(check = TRUE)
        
        GenomFile <- tempfile(fileext = ".txt")
        binFile <- tempfile(fileext = ".bed")
        
        command_1 <- "/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools idxstats %s | awk 'BEGIN {OFS=\"\\t\"} {if ($2>0) print ($1,$2)}' >  %s"
        system(sprintf(command_1, bamFile, GenomFile))
        
        if(slidingWindow=="YES"){
          command_2 <- "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools makewindows -g %s -w %s -s %s | sort -k1,1V -k2,2n > %s"
          system(sprintf(command_2, GenomFile, binSize, stepSize, binFile))
        } else {
          command_2 <- "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools makewindows -g %s -w %s | sort -k1,1V -k2,2n > %s"
          system(sprintf(command_2, GenomFile, binSize, binFile))
        }
        
        pncFiles_watson <- tempfile(fileext = ".bed")
        pncFiles_crick <- tempfile(fileext = ".bed")
        
        if(byReads_5p == TRUE){

          message(paste0("Calculating coverage with 5' ends of first-mate reads for", " ", 
                         tools::file_path_sans_ext(basename(bamFile)) ))
          
          # calculate coverage at watson strand by 5' end of the first mate reads
          command_3 <- paste0(
            "/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools view -h -@ 8 -q 30 -F 3840 -f 64 -L %s %s |",
            "grep -v XS:i: |",
            "/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools view -@ 8 -b - |",
            "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools genomecov -ibam stdin -strand + -d -5 |",
            "awk 'BEGIN {OFS=\"\\t\"} {if ($3>0) print $1,$2,$2,\"%s\",$3}' | sort -k1,1V -k2,2n > %s"
          )

          # calculate coverage at crick strand by 5' end of the first mate reads
          command_4 <- paste0(
            "/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools view -h -@ 8 -q 30 -F 3840 -f 64 -L %s %s |",
            "grep -v XS:i: |",
            "/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools view -@ 8 -b - |",
            "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools genomecov -ibam stdin -strand - -d -5 |",
            "awk 'BEGIN {OFS=\"\\t\"} {if ($3>0) print $1,$2,$2,\"%s\",$3}' | sort -k1,1V -k2,2n > %s"
          )
        }
        
        #
        # Watson (+)
        system(sprintf(command_3, binFile, bamFile,
                       paste0(tools::file_path_sans_ext(basename(bamFile)), "_watson"),
                       pncFiles_watson))
        
        # Crick (-)
        system(sprintf(command_4, binFile, bamFile,
                       paste0(tools::file_path_sans_ext(basename(bamFile)), "_crick"),
                       pncFiles_crick))
        
        #
        # sum the counts per bin for watson and crick separately and store in finFiles
        
        finFiles_watson <- paste0("~/Desktop/", Pro_1, "/", "Coverage", "/", Pro_1, "_", Pro_2, "_", "watson.bed")
        finFiles_crick <- paste0("~/Desktop/", Pro_1, "/", "Coverage", "/", Pro_1, "_", Pro_2, "_", "crick.bed")
        
        command_5 <- paste(
          "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools map",
          "-a %s -b %s -null 0 -o sum |",
          "awk 'BEGIN {OFS=\"\\t\"} {if ($4>=0) print $1,$2,$3,\"%s\",$4}' > %s"
        )
        
        # Watson (+)
        system(sprintf(command_5, binFile, pncFiles_watson,
                       paste0(tools::file_path_sans_ext(basename(bamFile)), "_watson"),
                       finFiles_watson))
        # Crick (-)
        system(sprintf(command_5, binFile, pncFiles_crick,
                       paste0(tools::file_path_sans_ext(basename(bamFile)), "_crick"),
                       finFiles_crick))
        
        #
        
        unlink(c(GenomFile, binFile, pncFiles_watson, pncFiles_crick), recursive = T, force = T)
        
      }
      
      for(i in 1:length(bamFiles)){
        BamCoverage(bamFile = bamFiles[i])
      }
      
     }
    
    ## Calculate Ratio
    #
    message("Calculating enrichment ratios!")
    #
    if(!dir.exists(paste0("~/Desktop/", Pro_1, "/", "Ratios"))){
      
      CoverageFiles <- Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Coverage", "/", "*", ".bed"))
      
      CalculateRatio <- function(IP_coverage, Input_coverage){
        
        IP.df = read.table(IP_coverage, header = F)
        In.df = read.table(Input_coverage, header = F)
        
        suppressWarnings(dir.create(paste0("~/Desktop/", Pro_1, "/", "Ratios"))) 
        
        
        IP_Sum <- sum(as.numeric(IP.df[,5]))
        In_Sum <- sum(as.numeric(In.df[,5]))
        corrFactor <- IP_Sum/In_Sum
        Ratio <- round(IP.df[,5]/In.df[,5]/corrFactor, 4)
        In.score.norm <- round(In.df[,5]*corrFactor)
        Ratio[!is.finite(Ratio)] <- 0
        
        strand <- cbind.data.frame(IP.df[,1], IP.df[,2], IP.df[,3], IP.df[,4], IP.df[,5], In.score.norm, Ratio)
        
        chroms <- unique(strand[,1])
        s_strand <- NULL
        for(i in 1:length(chroms)){
          chr <- strand[strand[,1]==chroms[i], ]
          x <- chr[,2]
          y <- chr$Ratio
          splineObject <- smooth.spline(x, y)
          chr$splineSmooth <- round(as.numeric(splineObject$y), 3)
          s_strand <- rbind.data.frame(s_strand, chr)
        }
        
        s_strand$ppois <- ppois(  q=s_strand[,5] - 1, 
                                  lambda=s_strand$`In.score.norm`, 
                                  lower.tail=FALSE, log=FALSE      )
        
        
        colnames(s_strand) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'ip.score', 'in.score', 'ratio', 'smooth', 'pvalue')
        
        ###
        
        write.table(s_strand, paste0("~/Desktop/", Pro_1, "/", "Ratios", "/", s_strand$name[1], ".bed"),
                    quote=FALSE, row.names=FALSE, sep="\t")
        
      }
      
      CalculateRatio(CoverageFiles[6], CoverageFiles[6])
      CalculateRatio(CoverageFiles[4], CoverageFiles[6])
      CalculateRatio(CoverageFiles[2], CoverageFiles[6])
      CalculateRatio(CoverageFiles[8], CoverageFiles[6])
      CalculateRatio(CoverageFiles[10], CoverageFiles[2])
      CalculateRatio(CoverageFiles[12], CoverageFiles[8])
      
      CalculateRatio(CoverageFiles[5], CoverageFiles[5])
      CalculateRatio(CoverageFiles[3], CoverageFiles[5])
      CalculateRatio(CoverageFiles[1], CoverageFiles[5])
      CalculateRatio(CoverageFiles[7], CoverageFiles[5])
      CalculateRatio(CoverageFiles[9], CoverageFiles[1])
      CalculateRatio(CoverageFiles[11], CoverageFiles[7])
      
    }
    
    ## Define and process peaks
    #
    message("Processing peaks ...")
    #
    if(!dir.exists(paste0("~/Desktop/", Pro_1, "/", "Peaks"))){
      
      bamFiles <- Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Bam", "/", "*", ".bam"))
      
      ProcessPeaks <- function(bamFiles){
        
        PeakFinder <- function(IPBam, InBam){
          
          Input_combined <- read.table(pipe(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools bamtobed -i %s", InBam)) )
          Input_Plus <- read.table(pipe(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools bamtobed -i %s | awk '{OFS=\"\\t\"} {if($6==\"+\") print $0}'", InBam)) )
          Input_Minus <- read.table(pipe(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools bamtobed -i %s | awk '{OFS=\"\\t\"} {if($6==\"-\") print $0}'", InBam)) )
          
          IP_combined <- read.table(pipe(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools bamtobed -i %s", IPBam)) )
          IP_Plus <- read.table(pipe(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools bamtobed -i %s | awk '{OFS=\"\\t\"} {if($6==\"+\") print $0}'", IPBam)) )
          IP_Minus <- read.table(pipe(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools bamtobed -i %s | awk '{OFS=\"\\t\"} {if($6==\"-\") print $0}'", IPBam)) )
          
          OriginPeaks <- function(IP_DF, Input_DF){
            
            tempdir(check = TRUE)
            IP <- tempfile(fileext=".bed")
            Input <- tempfile(fileext=".bed")
            outDir <- tempdir()
            peakFile <- tempfile(fileext = ".bed")
            ColHeads <- "\"chrom\\tpeakStart\\tpeakEnd\\tpeakLength\\tpeakSummit\\toriName\\toriStart\\toriEnd\""
            
            write.table(IP_DF, file = IP, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
            write.table(Input_DF, file = Input, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
            
            system(sprintf("macs2 callpeak -t %s -c %s -f BED -g 12157105 -p 10e-6 --nomodel -n %s --outdir %s 2> /dev/null", IP, Input, "Peak", outDir))
            allPeaks <- read.delim2(paste0(outDir, "/Peak_peaks.xls"), comment.char="#")
            
            #allPeaks <- suppressMessages(read.delim2(callpeak(tfile = IP, cfile = Input, gsize = 12157105, format = "BED", pvalue = 10e-6, nomodel = F, 
            #                                                  outdir = outDir, name = Pro_1)$outputs[3], comment.char="#"))
            
            write.table(allPeaks, file = peakFile, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
            
            Peaks_at_Origins <- read.table(pipe(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools intersect -wa -wb -a %s -b %s | awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$4,$5,$14,$12,$13}'", 
                                                        peakFile, All_Ori_Link, ColHeads)), header = TRUE ) 
            Peaks_at_Origins <- Peaks_at_Origins[!duplicated(Peaks_at_Origins$oriName), ]
            
            Peaks_at_Origins$oriCenter <- round((Peaks_at_Origins$oriStart + Peaks_at_Origins$oriEnd)/2)
            
            return(Peaks_at_Origins)
            
            unlink(c(IP, Input, peakFile))
            unlink(outDir, recursive = TRUE)
            
          }
          
          Watson <- OriginPeaks(IP_DF=IP_Plus, Input_DF=Input_Plus)
          Crick <- OriginPeaks(IP_DF=IP_Minus, Input_DF=Input_Minus)
          Combo <- OriginPeaks(IP_DF=IP_combined, Input_DF=Input_combined)
          
          ResPeaks <- list(Watson, Crick, Combo)
          names(ResPeaks) <- c("watsonPeaks", "crickPeaks", "comboPeaks")
          
          return(ResPeaks)
          
        }
        
        BrDU_Input <- PeakFinder(IPBam=bamFiles[1], InBam=bamFiles[3])
        
        #Define overlapping peaks
        
        tempdir(check = TRUE)
        
        pF_1 <- tempfile(fileext = ".bed")
        pF_2 <- tempfile(fileext = ".bed")
        
        BrDU_ColHeads <- "\"chrom\\tBWpStart\\tBWpEnd\\tBWpSummit\\tBCpStart\\tBCpEnd\\tBCpSummit\\toriName\\toriCenter\""
        
        #Overlaps between both watson and crick strands
        
        #BrDU
        write.table(BrDU_Input$watsonPeaks, file = pF_1, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
        write.table(BrDU_Input$crickPeaks, file = pF_2, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
        Overlapping_BrDU_Peaks <- read.table(pipe(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools intersect -wa -wb -a %s -b %s | awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$5,$11,$12,$14,$15,$18}'", 
                                                          pF_1, pF_2, BrDU_ColHeads)), header = TRUE ) 
        Overlapping_BrDU_Peaks <- Overlapping_BrDU_Peaks[!duplicated(Overlapping_BrDU_Peaks$oriName), ]
        
        
        unlink(c(pF_1, pF_2))
        
        #add name column to peakLists
        if(dim(Overlapping_BrDU_Peaks)[1]>0){
          Overlapping_BrDU_Peaks$name <- "BrDU"
        }
        
        
        ###Calculate Fork Locations from BrDU data
        SignalRanges <- function(PeakList){
          
          
          Chroms <- paste0("chr", as.roman(1:16))
          
          ForkPos <- NULL
          
          
          for(i in 1:16){
            
            PeakPos <- NULL
            
            i <- i
            
            Chr_Peaks <- PeakList[PeakList$chrom == Chroms[i], ]
            
            
            if(length(Chr_Peaks$chrom)==0) next 
            
            
            for(y in 1:length(Chr_Peaks$chrom)){
              
              y <- y
              
              if(Chr_Peaks$name[1]=="BrDU"){
                WatsonLeft <- round((Chr_Peaks$BWpSummit[y] - Chr_Peaks$BWpStart[y]))
                WatsonRight <- round((Chr_Peaks$BWpEnd[y] - Chr_Peaks$BWpSummit[y]))
                
                CrickLeft <- round((Chr_Peaks$BCpSummit[y] - Chr_Peaks$BCpStart[y]))
                CrickRight <- round((Chr_Peaks$BCpEnd[y] - Chr_Peaks$BCpSummit[y]))
                
                Peaks <- cbind.data.frame(WatsonLeft=WatsonLeft, WatsonRight=WatsonRight, 
                                          CrickLeft=CrickLeft, CrickRight=CrickRight, oriName=Chr_Peaks$oriName[y])
              }
              
              if(Chr_Peaks$name[1]=="BrDU_ChIP"){
                WatsonLeft <- round((Chr_Peaks$BWpSummit[y] - Chr_Peaks$CWpStart[y]))
                WatsonRight <- round((Chr_Peaks$CWpEnd[y] - Chr_Peaks$BWpSummit[y]))
                
                CrickLeft <- round((Chr_Peaks$BCpSummit[y] - Chr_Peaks$CCpStart[y]))
                CrickRight <- round((Chr_Peaks$CCpEnd[y] - Chr_Peaks$BCpSummit[y]))
                
                Peaks <- cbind.data.frame(WatsonLeft=WatsonLeft, WatsonRight=WatsonRight, 
                                          CrickLeft=CrickLeft, CrickRight=CrickRight, oriName=Chr_Peaks$oriName[y])
              }
              
              if(Chr_Peaks$name[1]=="BrDU_ChIP_eSPAN"){
                WatsonLeft <- round((Chr_Peaks$BWpSummit[y] - Chr_Peaks$EWpStart[y]))
                WatsonRight <- round((Chr_Peaks$EWpEnd[y] - Chr_Peaks$BWpSummit[y]))
                
                CrickLeft <- round((Chr_Peaks$BCpSummit[y] - Chr_Peaks$ECpStart[y]))
                CrickRight <- round((Chr_Peaks$ECpEnd[y] - Chr_Peaks$BCpSummit[y]))
                
                Peaks <- cbind.data.frame(WatsonLeft=WatsonLeft, WatsonRight=WatsonRight, 
                                          CrickLeft=CrickLeft, CrickRight=CrickRight, oriName=Chr_Peaks$oriName[y])
              }
              
              PeakPos <- rbind.data.frame(PeakPos, Peaks)
              PeakPos[PeakPos < 0] <- 0
              
            }
            
            ForkPos <- rbind.data.frame(ForkPos, PeakPos)
            
          }
          return(ForkPos)
        }
        
        #BrDU
        if(dim(Overlapping_BrDU_Peaks)[1]>0){
          ##Lagging and Leading strand lengths synthesised 
          ForkPos <- SignalRanges(Overlapping_BrDU_Peaks)
          
          LaggLeadSynthesis <- cbind.data.frame(lagging=ForkPos$CrickLeft+ForkPos$WatsonRight, 
                                                leading=ForkPos$WatsonLeft+ForkPos$CrickRight, 
                                                oriName=ForkPos$oriName)
          
          
          LeadingAverage <- round(mean(LaggLeadSynthesis$leading))
          LaggingAverage <- round(mean(LaggLeadSynthesis$lagging))
          LeadingSd <- round(sd(LaggLeadSynthesis$leading))
          LaggingSd <- round(sd(LaggLeadSynthesis$lagging))
          
          ###Estimate Averaging Window
          
          #Remove extreme outliers from the data
          Left <- c(ForkPos$WatsonLeft, ForkPos$CrickLeft)
          Right <- c(ForkPos$WatsonRight, ForkPos$CrickRight)
          
          Q <- quantile(Left, probs=0.99, na.rm = T)
          I <- IQR(Left)
          up  <-  Q + 1.5*I # Upper Range
          NewLeft <- Left[which(Left < up)]
          
          Q <- quantile(Right, probs=0.99, na.rm = T)
          I <- IQR(Right)
          up  <-  Q + 1.5*I # Upper Range
          NewRight <- Right[which(Right < up)]
          
          #Range for Averaging Window
          LeftLim <- max(NewLeft) #Left pan
          RightLim <- max(NewRight) #Right pan
          
          #Round the window
          LeftSide <- round(LeftLim/1000+0.5)*1000
          RightSide <- round(RightLim/1000+0.5)*1000
          
          if(LeftSide==RightSide){
            AveragingWindow <- LeftSide
          } else {
            AveragingWindow <- max(c(LeftSide, RightSide))
          } 
        }
        
        #
        binSize <- paste0( "300_500" )
        stepSize <- paste0( "10_50" )
        
        StatDat <- rbind.data.frame(LeftLim, RightLim, paste0(LeadingAverage, "±", LeadingSd), paste0(LaggingAverage, "±", LaggingSd), 
                                    AveragingWindow, binSize, stepSize)
        
        rownames(StatDat) <- c('Left', 'Right', "LeadSynthesis", "LaggSynthesis", "AveragingWindow", 'bin', 'slide')
        colnames(StatDat) <- " "
        
        ##
        dir.create(paste0("~/Desktop/", Pro_1, "/", "Peaks"), showWarnings = FALSE)
        
        write.table(BrDU_Input$comboPeaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "Primary_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
        write.table(Overlapping_BrDU_Peaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "BrDU_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
        write.table(StatDat, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "Synthesis.bed"), quote=FALSE, sep="\t")
        
        E_Ori <- read.table("~/Desktop/ngsAnalyser1.1.4/E_Rep.bed", header = F, quote = "/t")[ ,1:4]
        L_Ori <- read.table("~/Desktop/ngsAnalyser1.1.4/L_Rep.bed", header = F, quote = "/t")[ ,1:4]
        
        colnames(E_Ori) <- c("chrom", "chromStart", "chromEnd", "oriName"); E_Ori$oriCenter <- round((E_Ori$chromStart + E_Ori$chromEnd)/2); E_Ori$name <- "BrDU"
        colnames(L_Ori) <- c("chrom", "chromStart", "chromEnd", "oriName"); L_Ori$oriCenter <- round((L_Ori$chromStart + L_Ori$chromEnd)/2); L_Ori$name <- "BrDU"
        
        write.table(E_Ori, paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", "E_Ori", ".bed"), quote=FALSE, row.names=FALSE, sep="\t")
        write.table(L_Ori, paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", "L_Ori", ".bed"), quote=FALSE, row.names=FALSE, sep="\t")
        
      }
      
      ProcessPeaks(bamFiles)
      
      
    }
    
    ###
    
    rm(list=ls())
    gc()
    
    #
    message("Alignment & Primary Analysis complete!")
    #
  }
  
  Simulation_Calculation_Plotting <- function( NumOfSamples = 1, 
                                               Sample_1 ){
    
    # Read strand-wise ratio files from the previously generated data
    read_sample_files <- function(NumOfSamples,
                                  coverage_folder = "Coverage",
                                  peak_folder = "Peaks") {
      
      file_types <- c("Input", "BrDU", "ChIP", "eSPAN", "bSUP", "eSUP")
      strands <- c("watson", "crick")
      
      # Total file count (no versions anymore)
      total <- NumOfSamples * length(file_types) * length(strands)
      total <- total + NumOfSamples * 3  # three peak/info files per sample
      counter <- 0
      
      # === Loop over samples ===
      for (i in 1:NumOfSamples) {
        sample_path <- get(paste0("Sample_", i))
        sample_base <- basename(sample_path)
        
        # --- Load strand-wise coverage files ---
        for (ft in file_types) {
          for (strand in strands) {
            counter <- counter + 1
            
            file_name <- paste0(sample_base, "_", ft, "_", strand, ".bed")
            full_path <- file.path(sample_path, coverage_folder, file_name)
            var_name <- paste0("S", i, "_", ft, "_", strand)   # version removed
            
            if (file.exists(full_path)) {
              assign(var_name,
                     read.table(full_path, header = FALSE,
                                col.names = c("chrom", "chromStart", "chromEnd", "name", "ip.score")),
                     envir = .GlobalEnv)
            } else {
              warning(paste("Missing file:", full_path))
              assign(var_name, NULL, envir = .GlobalEnv)
            }
            
            if (counter %% 10 == 0 || counter == total)
              message(sprintf("Processed %d / %d files", counter, total))
          }
        }
        
        # --- Load peak files ---
        peak_dir <- file.path(sample_path, peak_folder)
        peak_base <- paste0(sample_base, "_")
        
        peak_files <- list(
          Peaks   = list(file = paste0(peak_base, "BrDU_Peaks.bed"),      header = TRUE),
          brPeaks = list(file = paste0(peak_base, "Primary_Peaks.bed"), header = TRUE),
          Info    = list(file = paste0(peak_base, "Synthesis.bed"),      header = FALSE)
        )
        
        for (pf in names(peak_files)) {
          counter <- counter + 1
          
          full_path <- file.path(peak_dir, peak_files[[pf]]$file)
          var_name <- paste0("S", i, "_", pf)
          
          if (file.exists(full_path)) {
            assign(var_name,
                   read.table(full_path, header = peak_files[[pf]]$header),
                   envir = .GlobalEnv)
          } else {
            warning(paste("Missing peak file:", full_path))
            assign(var_name, NULL, envir = .GlobalEnv)
          }
          
          if (counter %% 20 == 0 || counter == total)
            message(sprintf("Processed %d / %d files", counter, total))
        }
      }
      
      message("✅ All Coverage and Peak files loaded successfully.")
    }
    
    
    # Read coverage + peak data for all samples
    read_sample_files(NumOfSamples, coverage_folder = "Coverage", peak_folder = "Peaks")
    
    
    # Calculations
    
    Calc_data_simulation <- function(Sample){
      
      combine_sample_data <- function(sample_id) {
        file_types <- c("Input", "BrDU", "ChIP", "eSPAN", "bSUP", "eSUP")
        strands <- c("watson", "crick")
        
        sample_data <- list()
        
        # --- Combine coverage data (no versions anymore) ---
        for (ft in file_types) {
          for (strand in strands) {
            
            var_name <- paste0(sample_id, "_", ft, "_", strand)   # version removed
            
            if (exists(var_name, envir = .GlobalEnv)) {
              sample_data[[paste0(ft, "_", strand)]] <- get(var_name, envir = .GlobalEnv)
            } else {
              warning(paste("Missing:", var_name))
              sample_data[[paste0(ft, "_", strand)]] <- NULL
            }
          }
        }
        
        # --- Include peak and info files ---
        peak_names <- c("Peaks", "brPeaks", "Info")
        
        for (pf in peak_names) {
          var_name <- paste0(sample_id, "_", pf)
          
          if (exists(var_name, envir = .GlobalEnv)) {
            sample_data[[pf]] <- get(var_name, envir = .GlobalEnv)
          } else {
            warning(paste("Missing:", var_name))
            sample_data[[pf]] <- NULL
          }
        }
        
        message(paste("✅ Combined all data for", sample_id))
        return(sample_data)
      }
      
      assign(Sample, combine_sample_data(Sample))
      
      #fix averaging window
      useDef <- function(a,d) ifelse(isTruthy(a), a,d)
      
      Window = useDef(Window, "None")
      
      if(Window == "None"){
        Window = as.numeric(get(Sample)$Info[5, 2])
      } else {
        Window <- as.numeric(Window)
      }
      #
      if(PeakClass == "early"){ 
        PeakList = E_Ori
        names(PeakList)[names(PeakList) == "mid"] <- "oriCenter"
        names(PeakList)[names(PeakList) == "name"] <- "oriName"
        PeakList$name <- "Early-Origins"}
      if(PeakClass == "late"){ 
        PeakList = L_Ori
        names(PeakList)[names(PeakList) == "mid"] <- "oriCenter"
        names(PeakList)[names(PeakList) == "name"] <- "oriName"
        PeakList$name <- "Late-Origins" }
      if(PeakClass == "brdu"){ 
        PeakList = get(Sample)$Peaks 
        brPeaks <- get(Sample)$brPeaks
        ###
        RemList <- anti_join(brPeaks, PeakList, by = c("chrom", "oriName", "oriCenter"))
        NamList <- anti_join(brPeaks, RemList, by = c("chrom", "oriName", "oriCenter"))
        NamList <- NamList[!duplicated(NamList$oriCenter), ]
        PeakList$BrDUSummit <- NamList$peakSummit
        # #change oricenter to peaksummit for brdu peaks
        PeakList$oriCenter <- PeakList$BrDUSummit
      }
      #
      #
      message(paste0("calculating data and simulation results around fired replication origins"))
      #
      
      Data_and_RanSamp <- function(expSamp){
        
        V <- "ip.score"
        
        rowStat <- function(DF){
          
          Quantiles <- apply(DF, 1, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
          Ses <- apply(DF, 1, std.error)
          Means <- apply(DF, 1, mean)
          Sds <- apply(DF, 1, sd)
          
          return(list(q.25 = Quantiles["25%",],
                      Median = Quantiles["50%",],
                      q.75 = Quantiles["75%",],
                      Std.err = Ses,
                      Mean = Means,
                      Sd = Sds))
        }
        
        #
        ###
        PeakList$AvBstart <- PeakList$oriCenter - Window
        PeakList$AvBend <- PeakList$oriCenter + Window
        ####
        
        
        ###data
        
        data <- {
          
          
          chrS <- paste0("chr", as.roman(1:16))
          ##
          IP_T <- NULL
          IP_Tw <- NULL
          IP_Tc <- NULL
          IP_LagLead <- NULL
          
          for(i in 1:length(chrS)){
            
            i <- i
            
            #origins at given chromosomes
            OriList_ROs <- PeakList[PeakList$chrom == chrS[i], ]
            
            #300 bp bin with 10 bp step
            Crick_ROs <- get(Sample)[[paste0(expSamp, "_crick")]]
            Crick_ROs <- Crick_ROs[Crick_ROs$chrom == chrS[i], ]
            
            Watson_ROs <- get(Sample)[[paste0(expSamp, "_watson")]]
            Watson_ROs <- Watson_ROs[Watson_ROs$chrom == chrS[i], ]
            
            if(length(OriList_ROs$chrom)==0) next 
            
            IP_C <- NULL
            IP_Cwat <- NULL
            IP_Ccri <- NULL
            IP_Claglead <- NULL
            
            for(y in 1:length(OriList_ROs$chrom)){
              
              y <- y
              
              #ratio, wat_val, cri_val
              
              Crick_S <- Crick_ROs[Crick_ROs$chromStart>=OriList_ROs$AvBstart[y] & Crick_ROs$chromStart<=OriList_ROs$AvBend[y], ]
              Watson_S <- Watson_ROs[Watson_ROs$chromStart>=OriList_ROs$AvBstart[y] & Watson_ROs$chromStart<=OriList_ROs$AvBend[y], ]
              
              #check and define stepsize
              step_c <- Crick_S$chromStart[2]-Crick_S$chromStart[1]
              step_w <- Watson_S$chromStart[2]-Watson_S$chromStart[1]
              
              if(step_c == step_w){ 
                step <- step_c } else {
                  warning(print("stepSize not equal"))
                }
              
              IP_Z <- log2(Watson_S[,V] / Crick_S[,V]); IP_Z[!is.finite(IP_Z)] <- 0
              
              IP_wat <- Watson_S[,V]; IP_wat[!is.finite(IP_wat)] <- 0
              IP_cri <- Crick_S[,V]; IP_cri[!is.finite(IP_cri)] <- 0
              
              if(length(IP_Z)==round(2*Window/step)){
                Ratios <- IP_Z
                Wat_Val <- IP_wat
                Cri_Val <- IP_cri
              } 
              
              if(length(IP_Z) < round(2*Window/step)){
                if(length(1:(length(IP_Z)/2)) < round(2*Window/step)/2){
                  Ratios <- c(rep(0, (round(2*Window/step))-length(Crick_S$chromStart)), IP_Z )
                  Wat_Val <- c(rep(0, (round(2*Window/step))-length(Watson_S$chromStart)), IP_wat )
                  Cri_Val <- c(rep(0, (round(2*Window/step))-length(Crick_S$chromStart)), IP_cri )
                }
                if(length((length(IP_Z)/2+1):length(IP_Z)) < round(2*Window/step)/2){
                  Ratios <- c(IP_Z, rep(0, (round(2*Window/step))-length(Crick_S$chromStart)) )
                  Wat_Val <- c(IP_wat, rep(0, (round(2*Window/step))-length(Watson_S$chromStart)) )
                  Cri_Val <- c(IP_cri, rep(0, (round(2*Window/step))-length(Crick_S$chromStart)) )
                }
              } 
              
              if(length(IP_Z) > round(2*Window/step)){
                Ratios <- IP_Z[-c((round(2*Window/step)+1):length(IP_Z))]
                Wat_Val <- IP_wat[-c((round(2*Window/step)+1):length(IP_wat))]
                Cri_Val <- IP_cri[-c((round(2*Window/step)+1):length(IP_cri))]
              }
              
              #left-right/lagg-lead - 
              
              IP_Z_wat_left <- Watson_ROs[Watson_ROs$chromStart>=OriList_ROs$AvBstart[y] & Watson_ROs$chromStart<=OriList_ROs$oriCenter[y], ]
              IP_Z_wat_right <- Watson_ROs[Watson_ROs$chromStart>=OriList_ROs$oriCenter[y] & Watson_ROs$chromStart<=OriList_ROs$AvBend[y], ]
              
              IP_Z_cri_left <- Crick_ROs[Crick_ROs$chromStart>=OriList_ROs$AvBstart[y] & Crick_ROs$chromStart<=OriList_ROs$oriCenter[y], ]
              IP_Z_cri_right <- Crick_ROs[Crick_ROs$chromStart>=OriList_ROs$oriCenter[y] & Crick_ROs$chromStart<=OriList_ROs$AvBend[y], ]
              
              if(expSamp == "bSUP" || expSamp == "eSUP"){
                IP_Z_CL <- sum(IP_Z_cri_left[,V], na.rm = TRUE)
                IP_Z_CR <- sum(IP_Z_cri_right[,V], na.rm = TRUE)
                IP_Z_WL <- sum(IP_Z_wat_left[,V], na.rm = TRUE)
                IP_Z_WR <- sum(IP_Z_wat_right[,V], na.rm = TRUE)
              } else {
                IP_Z_CL <- sum(IP_Z_wat_left[,V], na.rm = TRUE)
                IP_Z_CR <- sum(IP_Z_wat_right[,V], na.rm = TRUE)
                IP_Z_WL <- sum(IP_Z_cri_left[,V], na.rm = TRUE)
                IP_Z_WR <- sum(IP_Z_cri_right[,V], na.rm = TRUE)
              }
              
              IP_Z_CL_WR <- IP_Z_CL + IP_Z_WR; IP_Z_CL_WR[!is.finite(IP_Z_CL_WR)] <- 0
              IP_Z_WL_CR <- IP_Z_WL + IP_Z_CR; IP_Z_WL_CR[!is.finite(IP_Z_WL_CR)] <- 0
              
              IP_Z_Pvalue <- binom.test(round(c(IP_Z_CL_WR, IP_Z_WL_CR)))$p.value
              
              IP_Z_LaLe <- IP_Z_CL_WR/IP_Z_WL_CR
              IP_Z_LaLe <- log2((abs(IP_Z_LaLe))^(sign(IP_Z_LaLe)))
              IP_Z_LaLe[!is.finite(IP_Z_LaLe)] <- 0
              
              #
              Rat_r <- matrix(Ratios, ncol = 1)
              Rat_w <- matrix(Wat_Val, ncol = 1)
              Rat_c <- matrix(Cri_Val, ncol = 1)
              LagLead <- cbind.data.frame(round(IP_Z_CL_WR), round(IP_Z_WL_CR), IP_Z_Pvalue, IP_Z_LaLe)
              #
              IP_C <- cbind(IP_C, Rat_r)
              IP_Cwat <- cbind(IP_Cwat, Rat_w)
              IP_Ccri <- cbind(IP_Ccri, Rat_c)
              IP_Claglead <- rbind.data.frame(IP_Claglead, LagLead)
              
            }
            
            IP_Claglead <- cbind.data.frame(OriList_ROs$chrom, OriList_ROs$oriName, IP_Claglead)
            colnames(IP_Claglead) <- c('chrom', 'name', "Lagg.sum", "Lead.sum", "p_value", "Bias")
            
            IP_T <- cbind(IP_T, IP_C); IP_T <- as.data.frame(IP_T)
            IP_Tw <- cbind(IP_Tw, IP_Cwat); IP_Tw <- as.data.frame(IP_Tw)
            IP_Tc <- cbind(IP_Tc, IP_Ccri); IP_Tc <- as.data.frame(IP_Tc)
            IP_LagLead <- rbind(IP_LagLead, IP_Claglead); IP_LagLead <- as.data.frame(IP_LagLead)
            
          }
          
          colnames(IP_T) <- c(1:length(PeakList$chrom))
          colnames(IP_Tw) <- c(1:length(PeakList$chrom))
          colnames(IP_Tc) <- c(1:length(PeakList$chrom))
          
          
          W2C_Median_data <- cbind.data.frame(median = rowStat(IP_T)$Median)
          #rownames(W2C_Median_data) <- paste0("bin", 1:length(W2C_Median_data[,1]))
          
          Wat_Median_data <- cbind.data.frame(median = rowStat(IP_Tw)$Median)
          #rownames(Wat_Median_data) <- paste0("bin", 1:length(Wat_Median_data[,1]))
          
          Cri_Median_data <- cbind.data.frame(median = rowStat(IP_Tc)$Median)
          #rownames(Cri_Median_data) <- paste0("bin", 1:length(Cri_Median_data[,1]))
          
          Lagg_Lead_Data <- IP_LagLead
          WatCri_Data <- cbind.data.frame(watson = Wat_Median_data$median, crick = Cri_Median_data$median, wat2cri = W2C_Median_data$median)
          rownames(WatCri_Data) <- paste0("bin", 1:length(WatCri_Data[,1]))
          
        }
        
        #Random Sampling
        
        simulation <- {
          
          
          Total_NumberOfelements <- length(PeakList$chrom)
          
          GetDistances <- function(Positions){
            return(Positions[2: length(Positions)] - Positions[1:(length(Positions)-1)])
          }
          
          sampleByMinDistance <- function(Vect, Size, MinDistance){
            
            Reps <- 0
            
            Samp <- TRUE
            while(Samp){
              
              Reps <- Reps + 1
              # if(Reps %% 1000 == 1)
              #   print(paste("Working on", Reps))
              
              x1 <- sort(round(runif(Size, min=0+MinDistance, max=Vect-MinDistance)))
              xdist <- GetDistances(x1)
              
              while(min(xdist)<MinDistance){
                # print('working' )
                x2 <- which.min(xdist) + sample(x = c(0,1), size = 1)
                x1 <- x1[-x2]
                x1 <- sort(c(x1, round(runif(1, min=0+MinDistance, max=Vect-MinDistance))))
                xdist <- GetDistances(x1)
              }
              
              #x1 <- c(0, x1, Vect)
              
              return(x1)
            }
          }
          
          ChromLengths <- seqlengths(Scerevisiae)[-17]
          
          ChromNames <- seqnames(Scerevisiae)[-17]
          
          Chrom_NumberOfelements <- c()
          for (i in 1:length(ChromNames)){
            Chrom_NumberOfelements <- c(Chrom_NumberOfelements, length(PeakList[PeakList$chrom==ChromNames[i], ]$chrom))
          }
          
          MinDists <- c()
          for (i in 1:length(ChromNames)){
            if(Chrom_NumberOfelements[i] <= 1){
              Dist <- 0
            } else {
              Dist <- min(GetDistances(PeakList[PeakList$chrom==ChromNames[i], ]$oriCenter))
            }
            MinDists <- c(MinDists, Dist)
          }
          
          
          NumberOfRepetitions <- iterations
          
          Averages_r <- as.data.frame(matrix(0, ncol = NumberOfRepetitions, nrow = (Window*2)/10))
          Averages_w <- as.data.frame(matrix(0, ncol = NumberOfRepetitions, nrow = (Window*2)/10))
          Averages_c <- as.data.frame(matrix(0, ncol = NumberOfRepetitions, nrow = (Window*2)/10))
          Averages_lagg <- as.data.frame(matrix(0, ncol = NumberOfRepetitions, nrow = Total_NumberOfelements))
          Averages_lead <- as.data.frame(matrix(0, ncol = NumberOfRepetitions, nrow = Total_NumberOfelements))
          Averages_bias <- as.data.frame(matrix(0, ncol = NumberOfRepetitions, nrow = Total_NumberOfelements))
          
          
          for(a in 1: NumberOfRepetitions){
            
            Sampled.Elements <- list()
            
            for(j in 1:length(ChromLengths)){
              
              print(paste("SAMPLE", "-", expSamp, ",", "Element", "-", unique(PeakList$name), ",", "Running simulation", a, "for chromosome", j))
              
              if(Chrom_NumberOfelements[j] == 0){
                SampElements <- NA
              }
              if(Chrom_NumberOfelements[j] == 1){
                SampElements <- sample(ChromLengths[j], 1)
              }
              if(Chrom_NumberOfelements[j] > 1){
                
                SampElements <- sampleByMinDistance(ChromLengths[j], Chrom_NumberOfelements[j], MinDists[j])
                UniVec <- length(unique(SampElements))
                
                if(length(SampElements)==UniVec){
                  Rep <- FALSE
                } else {
                  print("Repeating the sampling")
                }
                
              }
              
              Sampled.Elements[[j]] <- SampElements
            }
            names(Sampled.Elements) <- seqnames(Scerevisiae)[-17]
            
            ElementsList <- NULL
            for(k in 1:length(names(Sampled.Elements))){
              List <- data.frame(chrom=names(Sampled.Elements)[k], elements = Sampled.Elements[[k]])
              ElementsList <- rbind.data.frame(ElementsList, List)
            }
            
            ElementsList <- ElementsList[complete.cases(ElementsList), ]
            
            ElementsList$AvBstart <- ElementsList$elements - Window
            ElementsList$AvBend <- ElementsList$elements + Window
            
            #
            
            IP_T <- NULL
            IP_Tw <- NULL
            IP_Tc <- NULL
            IP_LagLead <- NULL
            
            for(i in 1:length(names(Sampled.Elements))){
              
              i <- i
              
              Elements <- ElementsList[ElementsList$chrom == names(Sampled.Elements)[i], ]
              
              #300 bp bin with 10 bp step
              Crick_ROs <- get(Sample)[[paste0(expSamp, "_crick")]]
              Crick_ROs <- Crick_ROs[Crick_ROs$chrom == chrS[i], ]
              
              Watson_ROs <- get(Sample)[[paste0(expSamp, "_watson")]]
              Watson_ROs <- Watson_ROs[Watson_ROs$chrom == chrS[i], ]
              
              
              if(length(Elements$chrom)==0) next
              
              IP_C <- NULL
              IP_Cwat <- NULL
              IP_Ccri <- NULL
              IP_Claglead <- NULL
              
              for(y in 1:length(Elements$chrom)){
                
                y <- y
                
                #ratio, wat_val, cri_val
                
                Crick_S <- Crick_ROs[Crick_ROs$chromStart>=Elements$AvBstart[y] & Crick_ROs$chromStart<=Elements$AvBend[y], ]
                Watson_S <- Watson_ROs[Watson_ROs$chromStart>=Elements$AvBstart[y] & Watson_ROs$chromStart<=Elements$AvBend[y], ]
                
                #check and define stepsize
                step_c <- Crick_S$chromStart[2]-Crick_S$chromStart[1]
                step_w <- Watson_S$chromStart[2]-Watson_S$chromStart[1]
                
                if(step_c == step_w){ 
                  step <- step_c } else {
                    warning(print("stepSize not equal"))
                  }
                
                
                IP_Z <- log2(Watson_S[,V] / Crick_S[,V]); IP_Z[!is.finite(IP_Z)] <- 0
                
                IP_wat <- Watson_S[,V]; IP_wat[!is.finite(IP_wat)] <- 0
                IP_cri <- Crick_S[,V]; IP_cri[!is.finite(IP_cri)] <- 0
                
                if(length(IP_Z)==round(2*Window/step)){
                  Ratios <- IP_Z
                  Wat_Val <- IP_wat
                  Cri_Val <- IP_cri
                } 
                
                if(length(IP_Z) < round(2*Window/step)){
                  if(length(1:(length(IP_Z)/2)) < round(2*Window/step)/2){
                    Ratios <- c(rep(0, (round(2*Window/step))-length(Crick_S$chromStart)), IP_Z )
                    Wat_Val <- c(rep(0, (round(2*Window/step))-length(Watson_S$chromStart)), IP_wat )
                    Cri_Val <- c(rep(0, (round(2*Window/step))-length(Crick_S$chromStart)), IP_cri )
                  }
                  if(length((length(IP_Z)/2+1):length(IP_Z)) < round(2*Window/step)/2){
                    Ratios <- c(IP_Z, rep(0, (round(2*Window/step))-length(Crick_S$chromStart)) )
                    Wat_Val <- c(IP_wat, rep(0, (round(2*Window/step))-length(Watson_S$chromStart)) )
                    Cri_Val <- c(IP_cri, rep(0, (round(2*Window/step))-length(Crick_S$chromStart)) )
                  }
                } 
                
                if(length(IP_Z) > round(2*Window/step)){
                  Ratios <- IP_Z[-c((round(2*Window/step)+1):length(IP_Z))]
                  Wat_Val <- IP_wat[-c((round(2*Window/step)+1):length(IP_wat))]
                  Cri_Val <- IP_cri[-c((round(2*Window/step)+1):length(IP_cri))]
                }
                
                #left-right/lagg-lead 
                
                IP_Z_wat_left <- Watson_ROs[Watson_ROs$chromStart>=Elements$AvBstart[y] & Watson_ROs$chromStart<=Elements$elements[y], ]
                IP_Z_wat_right <- Watson_ROs[Watson_ROs$chromStart>=Elements$elements[y] & Watson_ROs$chromStart<=Elements$AvBend[y], ]
                
                IP_Z_cri_left <- Crick_ROs[Crick_ROs$chromStart>=Elements$AvBstart[y] & Crick_ROs$chromStart<=Elements$elements[y], ]
                IP_Z_cri_right <- Crick_ROs[Crick_ROs$chromStart>=Elements$elements[y] & Crick_ROs$chromStart<=Elements$AvBend[y], ]
                
                if(expSamp == "bSUP" || expSamp == "eSUP"){
                  IP_Z_CL <- sum(IP_Z_cri_left[,V], na.rm = TRUE)
                  IP_Z_CR <- sum(IP_Z_cri_right[,V], na.rm = TRUE)
                  IP_Z_WL <- sum(IP_Z_wat_left[,V], na.rm = TRUE)
                  IP_Z_WR <- sum(IP_Z_wat_right[,V], na.rm = TRUE)
                } else {
                  IP_Z_CL <- sum(IP_Z_wat_left[,V], na.rm = TRUE)
                  IP_Z_CR <- sum(IP_Z_wat_right[,V], na.rm = TRUE)
                  IP_Z_WL <- sum(IP_Z_cri_left[,V], na.rm = TRUE)
                  IP_Z_WR <- sum(IP_Z_cri_right[,V], na.rm = TRUE)
                }
                
                IP_Z_CL_WR <- IP_Z_CL + IP_Z_WR; IP_Z_CL_WR[!is.finite(IP_Z_CL_WR)] <- 0
                IP_Z_WL_CR <- IP_Z_WL + IP_Z_CR; IP_Z_WL_CR[!is.finite(IP_Z_WL_CR)] <- 0
                
                IP_Z_Pvalue <- binom.test(round(c(IP_Z_CL_WR, IP_Z_WL_CR)))$p.value
                
                
                IP_Z_LaLe <- IP_Z_CL_WR/IP_Z_WL_CR
                IP_Z_LaLe <- log2((abs(IP_Z_LaLe))^(sign(IP_Z_LaLe)))
                IP_Z_LaLe[!is.finite(IP_Z_LaLe)] <- 0
                
                #
                
                Rat_r <- matrix(Ratios, ncol = 1)
                Rat_w <- matrix(Wat_Val, ncol = 1)
                Rat_c <- matrix(Cri_Val, ncol = 1)
                LagLead <- cbind.data.frame(round(IP_Z_CL_WR), round(IP_Z_WL_CR), IP_Z_Pvalue, IP_Z_LaLe)
                #
                
                IP_C <- cbind(IP_C, Rat_r)
                IP_Cwat <- cbind(IP_Cwat, Rat_w)
                IP_Ccri <- cbind(IP_Ccri, Rat_c)
                IP_Claglead <- rbind.data.frame(IP_Claglead, LagLead)
                
              }
              
              IP_Claglead <- cbind.data.frame(Elements$chrom, Elements$elements, IP_Claglead)
              colnames(IP_Claglead) <- c('chrom', 'center', "Lagg.sum", "Lead.sum", "p_value", "Bias")
              
              
              IP_T <- cbind(IP_T, IP_C); IP_T <- as.data.frame(IP_T)
              IP_Tw <- cbind(IP_Tw, IP_Cwat); IP_Tw <- as.data.frame(IP_Tw)
              IP_Tc <- cbind(IP_Tc, IP_Ccri); IP_Tc <- as.data.frame(IP_Tc)
              IP_LagLead <- rbind(IP_LagLead, IP_Claglead); IP_LagLead <- as.data.frame(IP_LagLead)
              
            }
            
            colnames(IP_T) <- c(1:length(PeakList$chrom))
            colnames(IP_Tw) <- c(1:length(PeakList$chrom))
            colnames(IP_Tc) <- c(1:length(PeakList$chrom))
            
            #w2c calculated ratio
            W2C_Median <- cbind.data.frame(median = rowStat(IP_T)$Median)
            Averages_r <- cbind.data.frame(Averages_r, W2C_Median)
            Averages_r <- Averages_r[,(dim(Averages_r)[2]-NumberOfRepetitions+1):dim(Averages_r)[2]]
            
            #watson enrichment
            Wat_Median <- cbind.data.frame(median = rowStat(IP_Tw)$Median)
            Averages_w <- cbind.data.frame(Averages_w, Wat_Median)
            Averages_w <- Averages_w[,(dim(Averages_w)[2]-NumberOfRepetitions+1):dim(Averages_w)[2]]
            
            #crick enrichment
            Cri_Median <- cbind.data.frame(median = rowStat(IP_Tc)$Median)
            Averages_c <- cbind.data.frame(Averages_c, Cri_Median)
            Averages_c <- Averages_c[,(dim(Averages_c)[2]-NumberOfRepetitions+1):dim(Averages_c)[2]]
            
            #lagging sum for elements sampled
            Lagg_sums <- cbind.data.frame(lagging = IP_LagLead$Lagg.sum)
            Averages_lagg <- cbind.data.frame(Averages_lagg, Lagg_sums)
            Averages_lagg <- Averages_lagg[,(dim(Averages_lagg)[2]-NumberOfRepetitions+1):dim(Averages_lagg)[2]]
            
            #leading sum for elements sampled
            Lead_sums <- cbind.data.frame(lagging = IP_LagLead$Lead.sum)
            Averages_lead <- cbind.data.frame(Averages_lead, Lead_sums)
            Averages_lead <- Averages_lead[,(dim(Averages_lead)[2]-NumberOfRepetitions+1):dim(Averages_lead)[2]]
            
            #lagg2lead calculated bias for elements sampled
            La2Le_bias <- cbind.data.frame(la2le_bias = IP_LagLead$Bias)
            Averages_bias <- cbind.data.frame(Averages_bias, La2Le_bias)
            Averages_bias <- Averages_bias[,(dim(Averages_bias)[2]-NumberOfRepetitions+1):dim(Averages_bias)[2]]
            
          }
          
          #w2c
          W2C_RanSamp <- rowMeans(Averages_r)
          W2C_RanSamp <- cbind.data.frame(w2c = W2C_RanSamp)
          rownames(W2C_RanSamp) <- paste0("bin", 1:length(W2C_RanSamp[,1]))
          
          #watson
          Wat_RanSamp <- rowMeans(Averages_w)
          Wat_RanSamp <- cbind.data.frame(watson = Wat_RanSamp)
          rownames(Wat_RanSamp) <- paste0("bin", 1:length(Wat_RanSamp[,1]))
          
          #crick
          Cri_RanSamp <- rowMeans(Averages_c)
          Cri_RanSamp <- cbind.data.frame(crick = Cri_RanSamp)
          rownames(Cri_RanSamp) <- paste0("bin", 1:length(Cri_RanSamp[,1]))
          
          #lagging
          Lagg_RanSamp <- rowMeans(Averages_lagg)
          Lagg_RanSamp <- cbind.data.frame(lagging = Lagg_RanSamp)
          #rownames(Lagg_RanSamp) <- paste0("element", 1:length(Lagg_RanSamp[,1]))
          
          #leading
          Lead_RanSamp <- rowMeans(Averages_lead)
          Lead_RanSamp <- cbind.data.frame(leading = Lead_RanSamp)
          #rownames(Lead_RanSamp) <- paste0("element", 1:length(Lead_RanSamp[,1]))
          
          #bias
          La2Le_RanSamp <- rowMeans(Averages_bias)
          La2Le_RanSamp <- cbind.data.frame(bias = La2Le_RanSamp)
          #rownames(La2Le_RanSamp) <- paste0("element", 1:length(La2Le_RanSamp[,1]))
          
          Lagg_Lead_RanSamp <- cbind.data.frame(chrom = IP_LagLead$chrom, lagging = Lagg_RanSamp$lagging, leading = Lead_RanSamp$leading, bias = La2Le_RanSamp$bias)
          WatCri_RanSamp <- cbind.data.frame(watson = Wat_RanSamp$watson, crick = Cri_RanSamp$crick, wat2cri = W2C_RanSamp$w2c)
          rownames(WatCri_RanSamp) <- paste0("bin", 1:length(WatCri_RanSamp[,1]))
          
          
        }
        
        ###
        
        Results <- list(LaLe_Data = Lagg_Lead_Data, WatCri_Data = WatCri_Data,
                        LaLe_Simu = Lagg_Lead_RanSamp, WatCri_Simu = WatCri_RanSamp)
        
        return(Results)
        
      }
      
      Input_results <- Data_and_RanSamp( expSamp = "Input")
      ChIP_results <- Data_and_RanSamp( expSamp = "ChIP")
      BrDU_results <- Data_and_RanSamp( expSamp = "BrDU")
      eSPAN_results <- Data_and_RanSamp( expSamp = "eSPAN")
      bSUP_results <- Data_and_RanSamp( expSamp = "bSUP")
      eSUP_results <- Data_and_RanSamp( expSamp = "eSUP")
      
      ####
      
      Res <- list(Input_results, ChIP_results, BrDU_results, eSPAN_results, bSUP_results, eSUP_results)
      names(Res) <- list('Input_results', 'ChIP_results', 'BrDU_results', 'eSPAN_results', 'bSUP_results', 'eSUP_results')
      
      return(Res)
      
    }
    
    for(i in 1:NumOfSamples){
      assign(paste0("S", i, "R"), Calc_data_simulation(Sample = paste0("S", i)))
    }
    
    
    # Plotting
    
    plot_results <- function(Sample){
      
      # read peak files, fix window and decide peaklist for plotting 
      combine_sample_peaks <- function(sample_id) {
        peak_names <- c("Peaks", "brPeaks", "Info")
        sample_data <- list()
        
        for (pf in peak_names) {
          var_name <- paste0(sample_id, "_", pf)
          
          if (exists(var_name, envir = .GlobalEnv)) {
            sample_data[[pf]] <- get(var_name, envir = .GlobalEnv)
          } else {
            warning(paste("⚠️ Missing:", var_name))
            sample_data[[pf]] <- NULL
          }
        }
        
        message(paste("✅ Combined peak files for", sample_id))
        return(sample_data)
      }
      
      assign(Sample, combine_sample_peaks(Sample))
      
      #fix averaging window
      useDef <- function(a,d) ifelse(isTruthy(a), a,d)
      
      Window = useDef(Window, "None")
      
      if(Window == "None"){
        Window = as.numeric(get(Sample)$Info[5, 2])
      } else {
        Window <- as.numeric(Window)
      }
      #
      #
      if(PeakClass == "early"){ 
        PeakList = E_Ori
        names(PeakList)[names(PeakList) == "mid"] <- "oriCenter"
        names(PeakList)[names(PeakList) == "name"] <- "oriName"
        PeakList$name <- "Early-Origins"}
      if(PeakClass == "late"){ 
        PeakList = L_Ori
        names(PeakList)[names(PeakList) == "mid"] <- "oriCenter"
        names(PeakList)[names(PeakList) == "name"] <- "oriName"
        PeakList$name <- "Late-Origins" }
      if(PeakClass == "brdu"){ 
        PeakList = get(Sample)$Peaks 
        brPeaks <- get(Sample)$brPeaks
        ###
        RemList <- anti_join(brPeaks, PeakList, by = c("chrom", "oriName", "oriCenter"))
        NamList <- anti_join(brPeaks, RemList, by = c("chrom", "oriName", "oriCenter"))
        NamList <- NamList[!duplicated(NamList$oriCenter), ]
        PeakList$BrDUSummit <- NamList$peakSummit
        # #change oricenter to peaksummit for brdu peaks
        PeakList$oriCenter <- PeakList$BrDUSummit
      }
      #
      
      # now define plotting device --- name the pdf output
      
      pdf(paste0(Sample_1, "/", basename(Sample_1), "_", "SUPeSPAN.pdf"), width = 12, height = 10)
      
      ## first page
      
      page_1 <- {
        
        #first page
        par(mar = c(5, 4, 4, 2) + 0.1)
        PlotMat <- {matrix(c(  
          0,0,0,0,0,0,0,0,0,
          0,1,1,1,1,1,1,1,0,
          0,1,1,1,1,1,1,1,0,
          0,2,2,2,2,3,3,3,0,
          0,2,2,2,2,3,3,3,0,
          0,2,2,2,2,3,3,3,0,
          0,4,4,4,4,4,4,4,0), 
          
          7,9,byrow=TRUE)}
        layout(PlotMat, c(1,1), c(1,1), TRUE)
        
        plot(NULL, xlim=c(0,2), ylim=c(0,2), ylab=" ", xlab=" ", yaxt="n", xaxt="n", bty = "n")
        txt <- paste0("Experiment: ", ExpTitle, "\n", date())
        text(1, 1, labels = txt, cex = 2.5, font = 2)
        plot(image, axes = FALSE)
        plot_venn <- {
          
          x <- list(
            BrDU = PeakList$oriName
          )
          
          futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
          VennObject <- venn.diagram(x, filename = NULL, main = "Origin Firing Events",
                                     main.fontface = "bold", main.cex = 1.5, 
                                     main.pos = c(0.5, 1.25),
                                     category.names = c("(BrDU_Peaks)"),
                                     # Circles
                                     lwd = 2,
                                     lty = 'blank',
                                     fill = adjustcolor("cornflowerblue", alpha.f = 0.45),
                                     # Numbers
                                     cex = 1.6,
                                     fontface = "italic",
                                     # Set names
                                     cat.cex = 1.25,
                                     cat.fontface = "bold",
                                     cat.default.pos = "outer",
                                     cat.dist = c(0.1))
          
          # Grid regions of current base plot (ie from frame)
          plot.new()  
          vps <- baseViewports()
          pushViewport(vps$inner, vps$figure, vps$plot)
          grid.draw(VennObject)
          popViewport(3)
          
          
        }
        plot(NULL, xlim=c(0,2), ylim=c(0,2), ylab=" ", xlab=" ", yaxt="n", xaxt="n", bty = "n")
        txt <- paste0("Aligner: Bowtie2", "; ", "Mode: first-mate read 5' end", "; ", "Centering: Peak-Summit")
        text(1, 1, labels = txt, cex = 1.25, font = 4)
        
        mtext("Page 1", side = 1, line = - 2, outer = TRUE, font = 3, cex = 1)
        
        
      }
      
      ##plot stranded enrichments
      
      page_2_and_3 <- {
        
        PlotEnrichments <- function(DataFile, PlotHeader, stepSize){
          
          Watson <- smooth.spline(1:length(DataFile$watson), DataFile$watson, spar = 0.5)$y
          Crick <- smooth.spline(1:length(DataFile$crick), DataFile$crick, spar = 0.5)$y
          
          Y <- max(round(abs(range(c(Watson, Crick*(-1))))+0.5))
          
          plot(Watson,
               ylim = c(-Y, Y),
               main = PlotHeader,
               ylab = "Average Enrichment", cex.main=1, xlab = "Distance from Origin Center (Kbp)", 
               xaxt = "n", col = 'brown3', type = 'l', lwd = 2, bty = 'n', las = 2, xaxs = 'i', yaxs = 'i')
          
          # lines(Wat25, lwd=0.5, col = adjustcolor("red", alpha.f = 0.2))
          # lines(Wat75, lwd=0.5, col = adjustcolor("red", alpha.f = 0.2))
          # 
          # polygon(x = c(1:length(Watson), rev(1:length(Watson))), 
          #         y = c(Wat25, rev(Wat75)), 
          #         col = adjustcolor("red", alpha.f = 0.2), border = NA)
          
          lines((Crick)*(-1), lwd=2, col = 'cornflowerblue', type = 'l')
          # lines((Cri25)*(-1), lwd=0.5, col = adjustcolor("cornflowerblue", alpha.f = 0.2))
          # lines((Cri75)*(-1), lwd=0.5, col = adjustcolor("cornflowerblue", alpha.f = 0.2))
          
          # polygon(x = c(1:length(Crick*(-1)), rev(1:length(Crick*(-1)))), 
          #         y = c(Cri25*(-1), rev(Cri75*(-1))), 
          #         col = adjustcolor("cornflowerblue", alpha.f = 0.2), border = NA)
          
          abline(h=0,lwd=0.4); abline(v=(length(Watson))/2,lwd=0.4)
          axisLabels <- seq(-Window,
                            +Window,
                            length.out = 9)
          axisLabels[c(2,4,6,8)] <- NA
          At <- (Window/stepSize)*seq(0,2,0.25); At[1] <- 1
          axis(1, at=At, labels = signif(axisLabels/1000, 2))
          
        }
        PlotEnrichments_norm <- function(DataFile, SimuFile, PlotHeader, stepSize){
          
          Wat_data <- DataFile$watson; Wat_simu <- SimuFile$watson
          Wat_simnorm <- Wat_data/Wat_simu; Wat_simnorm[!is.finite(Wat_simnorm)] <- 0
          
          Cri_data <- DataFile$crick; Cri_simu <- SimuFile$crick
          Cri_simnorm <- Cri_data/Cri_simu; Cri_simnorm[!is.finite(Cri_simnorm)] <- 0
          
          Watson <- smooth.spline(1:length(Wat_simnorm), Wat_simnorm, spar = 0.5)$y
          Crick <- smooth.spline(1:length(Cri_simnorm), Cri_simnorm, spar = 0.5)$y
          
          Y <- max(round(abs(range(c(Watson, Crick*(-1))))+0.5))
          
          plot(Watson,
               ylim = c(-Y, Y),
               main = PlotHeader,
               ylab = "Average Enrichment", cex.main=1, xlab = "Distance from Origin Center (Kbp)", 
               xaxt = "n", col = 'brown3', type = 'l', lwd = 2, bty = 'n', las = 2, xaxs = 'i', yaxs = 'i')
          
          # lines(Wat25, lwd=0.5, col = adjustcolor("red", alpha.f = 0.2))
          # lines(Wat75, lwd=0.5, col = adjustcolor("red", alpha.f = 0.2))
          # 
          # polygon(x = c(1:length(Watson), rev(1:length(Watson))), 
          #         y = c(Wat25, rev(Wat75)), 
          #         col = adjustcolor("red", alpha.f = 0.2), border = NA)
          
          lines((Crick)*(-1), lwd=2, col = 'cornflowerblue', type = 'l')
          # lines((Cri25)*(-1), lwd=0.5, col = adjustcolor("cornflowerblue", alpha.f = 0.2))
          # lines((Cri75)*(-1), lwd=0.5, col = adjustcolor("cornflowerblue", alpha.f = 0.2))
          
          # polygon(x = c(1:length(Crick*(-1)), rev(1:length(Crick*(-1)))), 
          #         y = c(Cri25*(-1), rev(Cri75*(-1))), 
          #         col = adjustcolor("cornflowerblue", alpha.f = 0.2), border = NA)
          
          abline(h=0,lwd=0.4); abline(v=(length(Watson))/2,lwd=0.4)
          axisLabels <- seq(-Window,
                            +Window,
                            length.out = 9)
          axisLabels[c(2,4,6,8)] <- NA
          At <- (Window/stepSize)*seq(0,2,0.25); At[1] <- 1
          axis(1, at=At, labels = signif(axisLabels/1000, 2))
          
        }
        PlotEnrichments_dnorm <- function(DataFile_ip, SimuFile_ip, DataFile_in, SimuFile_in, PlotHeader, stepSize){
          
          #
          Wat_data_ip <- DataFile_ip$watson; Wat_simu_ip <- SimuFile_ip$watson
          Wat_simnorm_ip <- Wat_data_ip/Wat_simu_ip; Wat_simnorm_ip[!is.finite(Wat_simnorm_ip)] <- 0
          Cri_data_ip <- DataFile_ip$crick; Cri_simu_ip <- SimuFile_ip$crick
          Cri_simnorm_ip <- Cri_data_ip/Cri_simu_ip; Cri_simnorm_ip[!is.finite(Cri_simnorm_ip)] <- 0
          #
          Wat_data_in <- DataFile_in$watson; Wat_simu_in <- SimuFile_in$watson
          Wat_simnorm_in <- Wat_data_in/Wat_simu_in; Wat_simnorm_in[!is.finite(Wat_simnorm_in)] <- 0
          Cri_data_in <- DataFile_in$crick; Cri_simu_in <- SimuFile_in$crick
          Cri_simnorm_in <- Cri_data_in/Cri_simu_in; Cri_simnorm_in[!is.finite(Cri_simnorm_in)] <- 0
          #
          Wat_ip_in <- Wat_simnorm_ip/Wat_simnorm_in; Wat_ip_in[!is.finite(Wat_ip_in)] <- 0
          Cri_ip_in <- Cri_simnorm_ip/Cri_simnorm_in; Cri_ip_in[!is.finite(Cri_ip_in)] <- 0
          
          ##
          Watson <- smooth.spline(1:length(Wat_ip_in), Wat_ip_in, spar = 0.5)$y
          Crick <- smooth.spline(1:length(Cri_ip_in), Cri_ip_in, spar = 0.5)$y
          
          Y <- max(round(abs(range(c(Watson, Crick*(-1))))+0.5))
          
          plot(Watson,
               ylim = c(-Y, Y),
               main = PlotHeader,
               ylab = "Average Enrichment", cex.main=1, xlab = "Distance from Origin Center (Kbp)", 
               xaxt = "n", col = 'brown3', type = 'l', lwd = 2, bty = 'n', las = 2, xaxs = 'i', yaxs = 'i')
          
          # lines(Wat25, lwd=0.5, col = adjustcolor("red", alpha.f = 0.2))
          # lines(Wat75, lwd=0.5, col = adjustcolor("red", alpha.f = 0.2))
          # 
          # polygon(x = c(1:length(Watson), rev(1:length(Watson))), 
          #         y = c(Wat25, rev(Wat75)), 
          #         col = adjustcolor("red", alpha.f = 0.2), border = NA)
          
          lines((Crick)*(-1), lwd=2, col = 'cornflowerblue', type = 'l')
          # lines((Cri25)*(-1), lwd=0.5, col = adjustcolor("cornflowerblue", alpha.f = 0.2))
          # lines((Cri75)*(-1), lwd=0.5, col = adjustcolor("cornflowerblue", alpha.f = 0.2))
          
          # polygon(x = c(1:length(Crick*(-1)), rev(1:length(Crick*(-1)))), 
          #         y = c(Cri25*(-1), rev(Cri75*(-1))), 
          #         col = adjustcolor("cornflowerblue", alpha.f = 0.2), border = NA)
          
          abline(h=0,lwd=0.4); abline(v=(length(Watson))/2,lwd=0.4)
          axisLabels <- seq(-Window,
                            +Window,
                            length.out = 9)
          axisLabels[c(2,4,6,8)] <- NA
          At <- (Window/stepSize)*seq(0,2,0.25); At[1] <- 1
          axis(1, at=At, labels = signif(axisLabels/1000, 2))
          
        }
        #second page
        
        PlotMat <- {matrix(c(     
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # 
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16), 
          
          20,20,byrow=TRUE)}
        layout(PlotMat, c(1,1), c(1,1), TRUE)
        
        #input
        PlotEnrichments(DataFile = S1R$Input_results$WatCri_Data, PlotHeader = "Input_data", stepSize = 10)
        PlotEnrichments(DataFile = S1R$Input_results$WatCri_Simu, PlotHeader = "Input_simu", stepSize = 10)
        PlotEnrichments_norm(DataFile = S1R$Input_results$WatCri_Data, SimuFile = S1R$Input_results$WatCri_Simu, PlotHeader = "Input_data_simu", stepSize = 10)
        PlotEnrichments_dnorm(DataFile_ip = S1R$Input_results$WatCri_Data, SimuFile_ip = S1R$Input_results$WatCri_Simu, 
                              DataFile_in = S1R$Input_results$WatCri_Data, SimuFile_in = S1R$Input_results$WatCri_Simu,
                              PlotHeader = "Input_Input", stepSize = 10)
        #chip
        PlotEnrichments(DataFile = S1R$ChIP_results$WatCri_Data, PlotHeader = "ChIP_data", stepSize = 10)
        PlotEnrichments(DataFile = S1R$ChIP_results$WatCri_Simu, PlotHeader = "ChIP_simu", stepSize = 10)
        PlotEnrichments_norm(DataFile = S1R$ChIP_results$WatCri_Data, SimuFile = S1R$ChIP_results$WatCri_Simu, PlotHeader = "ChIP_data_simu", stepSize = 10)
        PlotEnrichments_dnorm(DataFile_ip = S1R$ChIP_results$WatCri_Data, SimuFile_ip = S1R$ChIP_results$WatCri_Simu, 
                              DataFile_in = S1R$Input_results$WatCri_Data, SimuFile_in = S1R$Input_results$WatCri_Simu,
                              PlotHeader = "ChIP_Input", stepSize = 10)
        #BrDU
        PlotEnrichments(DataFile = S1R$BrDU_results$WatCri_Data, PlotHeader = "BrDU_data", stepSize = 10)
        PlotEnrichments(DataFile = S1R$BrDU_results$WatCri_Simu, PlotHeader = "BrDU_simu", stepSize = 10)
        PlotEnrichments_norm(DataFile = S1R$BrDU_results$WatCri_Data, SimuFile = S1R$BrDU_results$WatCri_Simu, PlotHeader = "BrDU_data_simu", stepSize = 10)
        PlotEnrichments_dnorm(DataFile_ip = S1R$BrDU_results$WatCri_Data, SimuFile_ip = S1R$BrDU_results$WatCri_Simu, 
                              DataFile_in = S1R$Input_results$WatCri_Data, SimuFile_in = S1R$Input_results$WatCri_Simu,
                              PlotHeader = "BrDU_Input", stepSize = 10)
        
        
        mtext("Read enrichments around fired origins by strands", side = 3, line = - 4, outer = TRUE, font = 2, cex = 2)
        
        mtext("Page 2", side = 1, line = - 2, outer = TRUE, font = 3, cex = 1)
        
        
        #third page
        
        #espan
        PlotEnrichments(DataFile = S1R$eSPAN_results$WatCri_Data, PlotHeader = "eSPAN_data", stepSize = 10)
        PlotEnrichments(DataFile = S1R$eSPAN_results$WatCri_Simu, PlotHeader = "eSPAN_simu", stepSize = 10)
        PlotEnrichments_norm(DataFile = S1R$eSPAN_results$WatCri_Data, SimuFile = S1R$eSPAN_results$WatCri_Simu, PlotHeader = "eSPAN_data_simu", stepSize = 10)
        PlotEnrichments_dnorm(DataFile_ip = S1R$eSPAN_results$WatCri_Data, SimuFile_ip = S1R$eSPAN_results$WatCri_Simu, 
                              DataFile_in = S1R$BrDU_results$WatCri_Data, SimuFile_in = S1R$BrDU_results$WatCri_Simu,
                              PlotHeader = "eSPAN_BrDU", stepSize = 10)
        
        #bsup
        PlotEnrichments(DataFile = S1R$bSUP_results$WatCri_Data, PlotHeader = "bSUP_data", stepSize = 10)
        PlotEnrichments(DataFile = S1R$bSUP_results$WatCri_Simu, PlotHeader = "bSUP_simu", stepSize = 10)
        PlotEnrichments_norm(DataFile = S1R$bSUP_results$WatCri_Data, SimuFile = S1R$bSUP_results$WatCri_Simu, PlotHeader = "bSUP_data_simu", stepSize = 10)
        PlotEnrichments_dnorm(DataFile_ip = S1R$bSUP_results$WatCri_Data, SimuFile_ip = S1R$bSUP_results$WatCri_Simu, 
                              DataFile_in = S1R$Input_results$WatCri_Data, SimuFile_in = S1R$Input_results$WatCri_Simu,
                              PlotHeader = "bSUP_Input", stepSize = 10)
        
        #esup
        PlotEnrichments(DataFile = S1R$eSUP_results$WatCri_Data, PlotHeader = "eSUP_data", stepSize = 10)
        PlotEnrichments(DataFile = S1R$eSUP_results$WatCri_Simu, PlotHeader = "eSUP_simu", stepSize = 10)
        PlotEnrichments_norm(DataFile = S1R$eSUP_results$WatCri_Data, SimuFile = S1R$eSUP_results$WatCri_Simu, PlotHeader = "eSUP_data_simu", stepSize = 10)
        PlotEnrichments_dnorm(DataFile_ip = S1R$eSUP_results$WatCri_Data, SimuFile_ip = S1R$eSUP_results$WatCri_Simu, 
                              DataFile_in = S1R$bSUP_results$WatCri_Data, SimuFile_in = S1R$bSUP_results$WatCri_Simu,
                              PlotHeader = "eSUP_bSUP", stepSize = 10)
        
        mtext("Read enrichments around fired origins by strands", side = 3, line = - 4, outer = TRUE, font = 2, cex = 2)
        
        mtext("Page 3", side = 1, line = - 2, outer = TRUE, font = 3, cex = 1)
        
        
      }
      
      ##plot strand biases
      
      page_4_and_5 <- {
        
        PlotAverages <- function(DataFile, PlotHeader, stepSize){
          
          w2c <- DataFile$wat2cri
          
          w2c <- smooth.spline(1:length(w2c), w2c, spar = 0.5)$y
          
          Y <- max(round(abs(range(w2c))+0.5))
          
          #Y <- 0.5
          
          plot(w2c,
               ylim = c(-Y, +Y),
               main = PlotHeader,
               ylab = "log2 watson/crick", cex.main=1, xlab = "Distance from Origin Center (Kbp)", 
               xaxt = "n", col = 'red', type = 'l', lwd = 2, bty = 'n', las = 2, xaxs='i', yaxs='i')
          
          # polygon(x = c(1:length(Med), rev(1:length(Med))), 
          #         y = c(q25, rev(q75)), 
          #         col = adjustcolor("blue", alpha.f = 0.2), border = NA)
          
          abline(h=0, lwd=0.4); abline(v=(length(w2c))/2,lwd=0.4)
          axisLabels <- seq(-Window,
                            +Window,
                            length.out = 9)
          axisLabels[c(2,4,6,8)] <- NA
          At <- (Window/stepSize)*seq(0,2,0.25); At[1] <- 1
          axis(1, at=At, labels = signif(axisLabels/1000, 2))
          
        }
        PlotAverages_inNorm <- function(DataFile_ip, DataFile_in, PlotHeader, stepSize){
          
          eps <- 1  # pseudocount to avoid division by zero
          
          log_ip_w <- log2((DataFile_ip$watson + eps) / (DataFile_in$watson + eps))
          log_ip_c <- log2((DataFile_ip$crick  + eps) / (DataFile_in$crick  + eps))
          
          w2c <- log_ip_w - log_ip_c
          
          w2c <- smooth.spline(1:length(w2c), w2c, spar = 0.6)$y
          
          
          Y <- max(round(abs(range(w2c))+0.5))
          
          
          #Y <- 0.5
          
          
          plot(w2c,
               ylim = c(-Y, +Y),
               main = PlotHeader,
               ylab = "log2 watson/crick", cex.main=1, xlab = "Distance from Origin Center (Kbp)", 
               xaxt = "n", col = 'red', type = 'l', lwd = 2, bty = 'n', las = 2, xaxs='i', yaxs='i')
          
          # polygon(x = c(1:length(Med), rev(1:length(Med))), 
          #         y = c(q25, rev(q75)), 
          #         col = adjustcolor("blue", alpha.f = 0.2), border = NA)
          
          abline(h=0, lwd=0.4); abline(v=(length(w2c))/2,lwd=0.4)
          axisLabels <- seq(-Window,
                            +Window,
                            length.out = 9)
          axisLabels[c(2,4,6,8)] <- NA
          At <- (Window/stepSize)*seq(0,2,0.25); At[1] <- 1
          axis(1, at=At, labels = signif(axisLabels/1000, 2))
          
        }
        PlotAverages_simNorm <- function(DataFile, SimuFile, PlotHeader, stepSize){
          
          eps <- 1  # small constant to avoid divide-by-zero
          
          log_data_w <- log2(DataFile$watson + eps)
          log_data_c <- log2(DataFile$crick  + eps)
          log_simu_w <- log2(SimuFile$watson + eps)
          log_simu_c <- log2(SimuFile$crick  + eps)
          
          # Watson–Crick difference, experiment vs simulation
          w2c <- (log_data_w - log_data_c) - (log_simu_w - log_simu_c)
          
          # Smooth for plotting
          w2c <- round(smooth.spline(seq_along(w2c), w2c, spar = 0.6)$y, 2)
          
          Y <- max(round(abs(range(w2c))+0.5))
          
          #Y <- 0.5
          
          plot(w2c,
               ylim = c(-Y, +Y),
               main = PlotHeader,
               ylab = "log2 watson/crick", cex.main=1, xlab = "Distance from Origin Center (Kbp)", 
               xaxt = "n", col = 'red', type = 'l', lwd = 2, bty = 'n', las = 2, xaxs='i', yaxs='i')
          
          # polygon(x = c(1:length(Med), rev(1:length(Med))), 
          #         y = c(q25, rev(q75)), 
          #         col = adjustcolor("blue", alpha.f = 0.2), border = NA)
          
          abline(h=0, lwd=0.4); abline(v=(length(w2c))/2,lwd=0.4)
          axisLabels <- seq(-Window,
                            +Window,
                            length.out = 9)
          axisLabels[c(2,4,6,8)] <- NA
          At <- (Window/stepSize)*seq(0,2,0.25); At[1] <- 1
          axis(1, at=At, labels = signif(axisLabels/1000, 2))
          
        }
        PlotAverages_dNorm <- function(DataFile_ip, SimuFile_ip, DataFile_in, SimuFile_in, PlotHeader, stepSize){
          
          eps <- 1  # small pseudocount to stabilize divisions
          
          # Compute log2 values for all tracks
          log_ip_w  <- log2(DataFile_ip$watson + eps)
          log_ip_c  <- log2(DataFile_ip$crick  + eps)
          log_in_w  <- log2(DataFile_in$watson + eps)
          log_in_c  <- log2(DataFile_in$crick  + eps)
          
          log_sim_ip_w <- log2(SimuFile_ip$watson + eps)
          log_sim_ip_c <- log2(SimuFile_ip$crick  + eps)
          log_sim_in_w <- log2(SimuFile_in$watson + eps)
          log_sim_in_c <- log2(SimuFile_in$crick  + eps)
          
          # Compute Watson–Crick asymmetry (W - C) for each dataset
          log_w2c_ip  <- (log_ip_w  - log_ip_c)  - (log_sim_ip_w  - log_sim_ip_c)
          log_w2c_in  <- (log_in_w  - log_in_c)  - (log_sim_in_w  - log_sim_in_c)
          
          # IP vs Input difference (equivalent to ratio of ratios)
          w2c <- log_w2c_ip - log_w2c_in
          
          # Smooth
          w2c <- round(smooth.spline(seq_along(w2c), w2c, spar = 0.6)$y, 2)
          
          # Clean up any infinities (shouldn't occur with pseudocount)
          w2c[!is.finite(w2c)] <- 0
          
          Y <- max(round(abs(range(w2c))+0.5))
          
          #Y <- 0.5
          
          plot(w2c,
               ylim = c(-Y, +Y),
               main = PlotHeader,
               ylab = "log2 watson/crick", cex.main=1, xlab = "Distance from Origin Center (Kbp)", 
               xaxt = "n", col = 'red', type = 'l', lwd = 2, bty = 'n', las = 2, xaxs='i', yaxs='i')
          
          # polygon(x = c(1:length(Med), rev(1:length(Med))), 
          #         y = c(q25, rev(q75)), 
          #         col = adjustcolor("blue", alpha.f = 0.2), border = NA)
          
          abline(h=0, lwd=0.4); abline(v=(length(w2c))/2,lwd=0.4)
          axisLabels <- seq(-Window,
                            +Window,
                            length.out = 9)
          axisLabels[c(2,4,6,8)] <- NA
          At <- (Window/stepSize)*seq(0,2,0.25); At[1] <- 1
          axis(1, at=At, labels = signif(axisLabels/1000, 2))
          
        }
        
        #fourth page
        
        PlotMat <- {matrix(c(     
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # 
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16), 
          
          20,20,byrow=TRUE)}
        layout(PlotMat, c(1,1), c(1,1), TRUE)
        
        #input
        PlotAverages(DataFile = S1R$Input_results$WatCri_Data, PlotHeader = "Input_data", stepSize = 10)
        PlotAverages_inNorm(DataFile_ip = S1R$Input_results$WatCri_Data, DataFile_in = S1R$Input_results$WatCri_Data, PlotHeader = "Input_Input_data", stepSize = 10)
        PlotAverages_simNorm(DataFile = S1R$Input_results$WatCri_Data, SimuFile = S1R$Input_results$WatCri_Simu, PlotHeader = "Input_data_simu", stepSize = 10)
        PlotAverages_dNorm(DataFile_ip = S1R$Input_results$WatCri_Data, SimuFile_ip = S1R$Input_results$WatCri_Simu, 
                           DataFile_in = S1R$Input_results$WatCri_Data, SimuFile_in = S1R$Input_results$WatCri_Simu,
                           PlotHeader = "Input_Input_data_simu", stepSize = 10)
        #chip
        PlotAverages(DataFile = S1R$ChIP_results$WatCri_Data, PlotHeader = "ChIP_data", stepSize = 10)
        PlotAverages_inNorm(DataFile_ip = S1R$ChIP_results$WatCri_Data, DataFile_in = S1R$Input_results$WatCri_Data, PlotHeader = "ChIP_Input_data", stepSize = 10)
        PlotAverages_simNorm(DataFile = S1R$ChIP_results$WatCri_Data, SimuFile = S1R$ChIP_results$WatCri_Simu, PlotHeader = "ChIP_data_simu", stepSize = 10)
        PlotAverages_dNorm(DataFile_ip = S1R$ChIP_results$WatCri_Data, SimuFile_ip = S1R$ChIP_results$WatCri_Simu, 
                           DataFile_in = S1R$Input_results$WatCri_Data, SimuFile_in = S1R$Input_results$WatCri_Simu,
                           PlotHeader = "ChIP_Input_data_simu", stepSize = 10)
        #BrDU
        PlotAverages(DataFile = S1R$BrDU_results$WatCri_Data, PlotHeader = "BrDU_data", stepSize = 10)
        PlotAverages_inNorm(DataFile_ip = S1R$BrDU_results$WatCri_Data, DataFile_in = S1R$Input_results$WatCri_Data, PlotHeader = "BrDU_Input_data", stepSize = 10)
        PlotAverages_simNorm(DataFile = S1R$BrDU_results$WatCri_Data, SimuFile = S1R$BrDU_results$WatCri_Simu, PlotHeader = "BrDU_data_simu", stepSize = 10)
        PlotAverages_dNorm(DataFile_ip = S1R$BrDU_results$WatCri_Data, SimuFile_ip = S1R$BrDU_results$WatCri_Simu, 
                           DataFile_in = S1R$Input_results$WatCri_Data, SimuFile_in = S1R$Input_results$WatCri_Simu,
                           PlotHeader = "BrDU_Input_data_simu", stepSize = 10)
        
        
        mtext("Watson over Crick biases around fired origins", side = 3, line = - 4, outer = TRUE, font = 2, cex = 2)
        mtext("Page 4", side = 1, line = - 2, outer = TRUE, font = 3, cex = 1)
        
        
        #fifth page
        
        #espan
        PlotAverages(DataFile = S1R$eSPAN_results$WatCri_Data, PlotHeader = "eSPAN_data", stepSize = 10)
        PlotAverages_inNorm(DataFile_ip = S1R$eSPAN_results$WatCri_Data, DataFile_in = S1R$BrDU_results$WatCri_Data, PlotHeader = "eSPAN_BrDU_data", stepSize = 10)
        PlotAverages_simNorm(DataFile = S1R$eSPAN_results$WatCri_Data, SimuFile = S1R$eSPAN_results$WatCri_Simu, PlotHeader = "eSPAN_data_simu", stepSize = 10)
        PlotAverages_dNorm(DataFile_ip = S1R$eSPAN_results$WatCri_Data, SimuFile_ip = S1R$eSPAN_results$WatCri_Simu, 
                           DataFile_in = S1R$BrDU_results$WatCri_Data, SimuFile_in = S1R$BrDU_results$WatCri_Simu,
                           PlotHeader = "eSPAN_BrDU_data_simu", stepSize = 10)
        
        #bsup
        PlotAverages(DataFile = S1R$bSUP_results$WatCri_Data, PlotHeader = "bSUP_data", stepSize = 10)
        PlotAverages_inNorm(DataFile_ip = S1R$bSUP_results$WatCri_Data, DataFile_in = S1R$Input_results$WatCri_Data, PlotHeader = "bSUP_Input_data", stepSize = 10)
        PlotAverages_simNorm(DataFile = S1R$bSUP_results$WatCri_Data, SimuFile = S1R$bSUP_results$WatCri_Simu, PlotHeader = "bSUP_data_simu", stepSize = 10)
        PlotAverages_dNorm(DataFile_ip = S1R$bSUP_results$WatCri_Data, SimuFile_ip = S1R$bSUP_results$WatCri_Simu, 
                           DataFile_in = S1R$Input_results$WatCri_Data, SimuFile_in = S1R$Input_results$WatCri_Simu,
                           PlotHeader = "bSUP_Input_data_simu", stepSize = 10)
        
        #esup
        PlotAverages(DataFile = S1R$eSUP_results$WatCri_Data, PlotHeader = "eSUP_data", stepSize = 10)
        PlotAverages_inNorm(DataFile_ip = S1R$eSUP_results$WatCri_Data, DataFile_in = S1R$bSUP_results$WatCri_Data, PlotHeader = "eSUP_bSUP_data", stepSize = 10)
        PlotAverages_simNorm(DataFile = S1R$eSUP_results$WatCri_Data, SimuFile = S1R$eSUP_results$WatCri_Simu, PlotHeader = "eSUP_data_simu", stepSize = 10)
        PlotAverages_dNorm(DataFile_ip = S1R$eSUP_results$WatCri_Data, SimuFile_ip = S1R$eSUP_results$WatCri_Simu, 
                           DataFile_in = S1R$bSUP_results$WatCri_Data, SimuFile_in = S1R$bSUP_results$WatCri_Simu,
                           PlotHeader = "eSUP_bSUP_data_simu", stepSize = 10)
        
        mtext("Watson over Crick biases around fired origins", side = 3, line = - 4, outer = TRUE, font = 2, cex = 2)
        mtext("Page 5", side = 1, line = - 2, outer = TRUE, font = 3, cex = 1)
        
        
      }
      
      ##dotplots
      
      page_6_and_7 <- {
        
        PlotDots <- function(DataFile, PlotHeader){
          
          dotV <- DataFile$Bias; dotV[!is.finite(dotV)] <- 0
          dotP <- DataFile$p_value
          
          Y <- max(round(abs(range(dotV))+0.5))
          
          boxplot(dotV, ylim = c(-Y,Y), las = 2, ylab = 'log2 lagging/leading', 
                  cex.main=1, xlab = " ", bty = 'n', main=PlotHeader, col=adjustcolor("grey", alpha.f = 0.25), lwd=0.5)
          
          if(length(dotV[which(dotP > 10e-6)])>0){
            spreadPoints(values=dotV[which(dotP > 10e-6)], position=1.0, pointCex=0.65, col="blue", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          if(length(dotV[which(dotV >=0 & dotP <= 10e-6)])>0){
            spreadPoints(values=dotV[which(dotV >=0 & dotP <= 10e-6)], position=1.0, pointCex=0.65, col="orange", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          if(length(dotV[which(dotV < 0 & dotP <= 10e-6)])>0){
            spreadPoints(values=dotV[which(dotV < 0 & dotP <= 10e-6)], position=1.0, pointCex=0.65, col="purple", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          
          
          Lagg <- length(dotV[which(dotV >=0 & dotP <= 10e-6)])
          Lead <- length(dotV[which(dotV < 0 & dotP <= 10e-6)])
          Indt <- length(dotV[which(dotP > 10e-6)])
          
          legend("bottomright", legend = c(paste0("lagg", " (", Lagg, ")"), 
                                           paste0("lead", " (", Lead, ")"), 
                                           paste0("inde", " (", Indt, ")")), 
                 col = c("orange", "purple", "blue"), pch = 19, pt.cex=0.65, bty = "n", cex = 0.65)
          
          
          InName <- PlotHeader
          
          #calculate p value
          #Decision Tree
          biasP <- binom.test(round(c(Lagg+Lead, Indt)))$p.value
          
          if(biasP <= 10e-3 & (Lagg+Lead) > Indt){
            
            biasQ <- binom.test(round(c(Lagg, Lead)))$p.value
            
            if(biasQ <= 10e-6){
              if(Lagg > Lead){
                conc <- paste0("Strong bias for lagging synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
              if(Lead > Lagg){
                conc <- paste0("Strong bias for leading synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-6 & biasQ <= 10e-4){
              if(Lagg > Lead){
                conc <- paste0("Weak bias for lagging synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
              if(Lead > Lagg){
                conc <- paste0("Weak bias for leading synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-4  & biasQ <= 10e-2){
              if(Lagg > Lead){
                conc <- paste0("Very weak bias for lagging synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
              if(Lead > Lagg){
                conc <- paste0("Very weak bias for leading synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-2){
              conc <- paste0("No significant strandedness in DNA synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
            }
          } else {
            conc <- paste0("No decisive bias pattern detected", "\n","too many origins show indeterminable bias")
          }
          
          mtext(conc, side = 1, line = 2, cex = 0.65)
        }
        PlotDots_inNorm <- function(DataFile_ip, DataFile_in, PlotHeader){
          
          lagg_ip <- DataFile_ip$Lagg.sum; lagg_in <- DataFile_in$Lagg.sum
          lagg_ip_in <- lagg_ip/lagg_in; lagg_ip_in[!is.finite(lagg_ip_in)] <- 0
          
          lead_ip <- DataFile_ip$Lead.sum; lead_in <- DataFile_in$Lead.sum
          lead_ip_in <- lead_ip/lead_in; lead_ip_in[!is.finite(lead_ip_in)] <- 0
          
          dotZ <- (lagg_ip_in - lead_ip_in) / sqrt( (lagg_ip_in / lagg_in) + (lead_ip_in / lead_in) )
          dotP <- 2 * pnorm(-abs(dotZ))
          
          dotV <- lagg_ip_in/lead_ip_in; dotV <- log2((abs(dotV))^(sign(dotV)))
          
          Y <- max(round(abs(range(dotV))+0.5))
          
          boxplot(dotV, ylim = c(-Y,Y), las = 2, ylab = 'log2 lagging/leading', 
                  cex.main=1, xlab = " ", bty = 'n', main=PlotHeader, col=adjustcolor("grey", alpha.f = 0.25), lwd=0.5)
          
          if(length(dotV[which(dotP > 10e-6)])>0){
            spreadPoints(values=dotV[which(dotP > 10e-6)], position=1.0, pointCex=0.65, col="blue", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          if(length(dotV[which(dotV >=0 & dotP <= 10e-6)])>0){
            spreadPoints(values=dotV[which(dotV >=0 & dotP <= 10e-6)], position=1.0, pointCex=0.65, col="orange", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          if(length(dotV[which(dotV < 0 & dotP <= 10e-6)])>0){
            spreadPoints(values=dotV[which(dotV < 0 & dotP <= 10e-6)], position=1.0, pointCex=0.65, col="purple", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          
          
          Lagg <- length(dotV[which(dotV >=0 & dotP <= 10e-6)])
          Lead <- length(dotV[which(dotV < 0 & dotP <= 10e-6)])
          Indt <- length(dotV[which(dotP > 10e-6)])
          
          legend("bottomright", legend = c(paste0("lagg", " (", Lagg, ")"), 
                                           paste0("lead", " (", Lead, ")"), 
                                           paste0("inde", " (", Indt, ")")), 
                 col = c("orange", "purple", "blue"), pch = 19, pt.cex=0.65, bty = "n", cex = 0.65)
          
          
          InName <- PlotHeader
          
          #calculate p value
          #Decision Tree
          biasP <- binom.test(round(c(Lagg+Lead, Indt)))$p.value
          
          if(biasP <= 10e-3 & (Lagg+Lead) > Indt){
            
            biasQ <- binom.test(round(c(Lagg, Lead)))$p.value
            
            if(biasQ <= 10e-6){
              if(Lagg > Lead){
                conc <- paste0("Strong bias for lagging synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
              if(Lead > Lagg){
                conc <- paste0("Strong bias for leading synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-6 & biasQ <= 10e-4){
              if(Lagg > Lead){
                conc <- paste0("Weak bias for lagging synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
              if(Lead > Lagg){
                conc <- paste0("Weak bias for leading synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-4  & biasQ <= 10e-2){
              if(Lagg > Lead){
                conc <- paste0("Very weak bias for lagging synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
              if(Lead > Lagg){
                conc <- paste0("Very weak bias for leading synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-2){
              conc <- paste0("No significant strandedness in DNA synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
            }
          } else {
            conc <- paste0("No decisive bias pattern detected", "\n","too many origins show indeterminable bias")
          }
          
          mtext(conc, side = 1, line = 2, cex = 0.65)
        }
        PlotDots_simNorm <- function(DataFile, SimuFile, PlotHeader){
          
          lagg_data <- DataFile$Lagg.sum; lagg_simu <- SimuFile$lagging
          lagg_data_simu <- lagg_data/lagg_simu; lagg_data_simu[!is.finite(lagg_data_simu)] <- 0
          
          lead_data <- DataFile$Lead.sum; lead_simu <- SimuFile$leading
          lead_data_simu <- lead_data/lead_simu; lead_data_simu[!is.finite(lead_data_simu)] <- 0
          
          dotZ <- (lagg_data_simu - lead_data_simu) / sqrt( (lagg_data_simu / lagg_simu) + (lead_data_simu / lead_simu) )
          dotP <- 2 * pnorm(-abs(dotZ))
          
          dotV <- lagg_data_simu/lead_data_simu; dotV <- log2((abs(dotV))^(sign(dotV)))
          
          Y <- max(round(abs(range(dotV))+0.5))
          
          boxplot(dotV, ylim = c(-Y,Y), las = 2, ylab = 'log2 lagging/leading', 
                  cex.main=1, xlab = " ", bty = 'n', main=PlotHeader, col=adjustcolor("grey", alpha.f = 0.25), lwd=0.5)
          
          if(length(dotV[which(dotP > 10e-6)])>0){
            spreadPoints(values=dotV[which(dotP > 10e-6)], position=1.0, pointCex=0.65, col="blue", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          if(length(dotV[which(dotV >=0 & dotP <= 10e-6)])>0){
            spreadPoints(values=dotV[which(dotV >=0 & dotP <= 10e-6)], position=1.0, pointCex=0.65, col="orange", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          if(length(dotV[which(dotV < 0 & dotP <= 10e-6)])>0){
            spreadPoints(values=dotV[which(dotV < 0 & dotP <= 10e-6)], position=1.0, pointCex=0.65, col="purple", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          
          
          Lagg <- length(dotV[which(dotV >=0 & dotP <= 10e-6)])
          Lead <- length(dotV[which(dotV < 0 & dotP <= 10e-6)])
          Indt <- length(dotV[which(dotP > 10e-6)])
          
          legend("bottomright", legend = c(paste0("lagg", " (", Lagg, ")"), 
                                           paste0("lead", " (", Lead, ")"), 
                                           paste0("inde", " (", Indt, ")")), 
                 col = c("orange", "purple", "blue"), pch = 19, pt.cex=0.65, bty = "n", cex = 0.65)
          
          
          InName <- PlotHeader
          
          #calculate p value
          #Decision Tree
          biasP <- binom.test(round(c(Lagg+Lead, Indt)))$p.value
          
          if(biasP <= 10e-3 & (Lagg+Lead) > Indt){
            
            biasQ <- binom.test(round(c(Lagg, Lead)))$p.value
            
            if(biasQ <= 10e-6){
              if(Lagg > Lead){
                conc <- paste0("Strong bias for lagging synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
              if(Lead > Lagg){
                conc <- paste0("Strong bias for leading synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-6 & biasQ <= 10e-4){
              if(Lagg > Lead){
                conc <- paste0("Weak bias for lagging synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
              if(Lead > Lagg){
                conc <- paste0("Weak bias for leading synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-4  & biasQ <= 10e-2){
              if(Lagg > Lead){
                conc <- paste0("Very weak bias for lagging synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
              if(Lead > Lagg){
                conc <- paste0("Very weak bias for leading synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-2){
              conc <- paste0("No significant strandedness in DNA synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
            }
          } else {
            conc <- paste0("No decisive bias pattern detected", "\n","too many origins show indeterminable bias")
          }
          
          mtext(conc, side = 1, line = 2, cex = 0.65)
        }
        PlotDots_dNorm <- function(DataFile_ip, SimuFile_ip, DataFile_in, SimuFile_in, PlotHeader){
          
          #ip
          lagg_data_ip <- DataFile_ip$Lagg.sum; lagg_simu_ip <- SimuFile_ip$lagging
          lagg_data_simu_ip <- lagg_data_ip/lagg_simu_ip; lagg_data_simu_ip[!is.finite(lagg_data_simu_ip)] <- 0
          
          lead_data_ip <- DataFile_ip$Lead.sum; lead_simu_ip <- SimuFile_ip$leading
          lead_data_simu_ip <- lead_data_ip/lead_simu_ip; lead_data_simu_ip[!is.finite(lead_data_simu_ip)] <- 0
          #in
          lagg_data_in <- DataFile_in$Lagg.sum; lagg_simu_in <- SimuFile_in$lagging
          lagg_data_simu_in <- lagg_data_in/lagg_simu_in; lagg_data_simu_in[!is.finite(lagg_data_simu_in)] <- 0
          
          lead_data_in <- DataFile_in$Lead.sum; lead_simu_in <- SimuFile_in$leading
          lead_data_simu_in <- lead_data_in/lead_simu_in; lead_data_simu_in[!is.finite(lead_data_simu_in)] <- 0
          #
          lagg_ip_in <- lagg_data_simu_ip/lagg_data_simu_in; lagg_ip_in[!is.finite(lagg_ip_in)] <- 0
          lead_ip_in <- lead_data_simu_ip/lead_data_simu_in; lead_ip_in[!is.finite(lead_ip_in)] <- 0
          #
          dotZ <- (lagg_ip_in - lead_ip_in) / sqrt( (lagg_ip_in / lagg_simu_in) + (lead_ip_in / lead_simu_in) )
          dotP <- 2 * pnorm(-abs(dotZ))
          
          dotV <- lagg_ip_in/lead_ip_in; dotV <- log2((abs(dotV))^(sign(dotV)))
          
          Y <- max(round(abs(range(dotV))+0.5))
          
          boxplot(dotV, ylim = c(-Y,Y), las = 2, ylab = 'log2 lagging/leading', 
                  cex.main=1, xlab = " ", bty = 'n', main=PlotHeader, col=adjustcolor("grey", alpha.f = 0.25), lwd=0.5)
          
          if(length(dotV[which(dotP > 10e-6)])>0){
            spreadPoints(values=dotV[which(dotP > 10e-6)], position=1.0, pointCex=0.65, col="blue", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          if(length(dotV[which(dotV >=0 & dotP <= 10e-6)])>0){
            spreadPoints(values=dotV[which(dotV >=0 & dotP <= 10e-6)], position=1.0, pointCex=0.65, col="orange", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          if(length(dotV[which(dotV < 0 & dotP <= 10e-6)])>0){
            spreadPoints(values=dotV[which(dotV < 0 & dotP <= 10e-6)], position=1.0, pointCex=0.65, col="purple", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          
          
          Lagg <- length(dotV[which(dotV >=0 & dotP <= 10e-6)])
          Lead <- length(dotV[which(dotV < 0 & dotP <= 10e-6)])
          Indt <- length(dotV[which(dotP > 10e-6)])
          
          legend("bottomright", legend = c(paste0("lagg", " (", Lagg, ")"), 
                                           paste0("lead", " (", Lead, ")"), 
                                           paste0("inde", " (", Indt, ")")), 
                 col = c("orange", "purple", "blue"), pch = 19, pt.cex=0.65, bty = "n", cex = 0.65)
          
          
          InName <- PlotHeader
          
          #calculate p value
          #Decision Tree
          biasP <- binom.test(round(c(Lagg+Lead, Indt)))$p.value
          
          if(biasP <= 10e-3 & (Lagg+Lead) > Indt){
            
            biasQ <- binom.test(round(c(Lagg, Lead)))$p.value
            
            if(biasQ <= 10e-6){
              if(Lagg > Lead){
                conc <- paste0("Strong bias for lagging synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
              if(Lead > Lagg){
                conc <- paste0("Strong bias for leading synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-6 & biasQ <= 10e-4){
              if(Lagg > Lead){
                conc <- paste0("Weak bias for lagging synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
              if(Lead > Lagg){
                conc <- paste0("Weak bias for leading synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-4  & biasQ <= 10e-2){
              if(Lagg > Lead){
                conc <- paste0("Very weak bias for lagging synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
              if(Lead > Lagg){
                conc <- paste0("Very weak bias for leading synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-2){
              conc <- paste0("No significant strandedness in DNA synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
            }
          } else {
            conc <- paste0("No decisive bias pattern detected", "\n","too many origins show indeterminable bias")
          }
          
          mtext(conc, side = 1, line = 2, cex = 0.65)
        }
        
        
        #sixth page
        
        PlotMat <- {matrix(c(     
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # 
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16), 
          
          20,20,byrow=TRUE)}
        layout(PlotMat, c(1,1), c(1,1), TRUE)
        
        #input
        PlotDots(DataFile = S1R$Input_results$LaLe_Data, PlotHeader = "Input_data")
        PlotDots_inNorm(DataFile_ip = S1R$Input_results$LaLe_Data, DataFile_in = S1R$Input_results$LaLe_Data, PlotHeader = "Input_Input_data")
        PlotDots_simNorm(DataFile = S1R$Input_results$LaLe_Data, SimuFile = S1R$Input_results$LaLe_Simu, PlotHeader = "Input_data_simu")
        PlotDots_dNorm(DataFile_ip = S1R$Input_results$LaLe_Data, SimuFile_ip = S1R$Input_results$LaLe_Simu, 
                       DataFile_in = S1R$Input_results$LaLe_Data, SimuFile_in = S1R$Input_results$LaLe_Simu,
                       PlotHeader = "Input_Input_data_simu")
        #chip
        PlotDots(DataFile = S1R$ChIP_results$LaLe_Data, PlotHeader = "ChIP_data")
        PlotDots_inNorm(DataFile_ip = S1R$ChIP_results$LaLe_Data, DataFile_in = S1R$Input_results$LaLe_Data, PlotHeader = "ChIP_Input_data")
        PlotDots_simNorm(DataFile = S1R$ChIP_results$LaLe_Data, SimuFile = S1R$ChIP_results$LaLe_Simu, PlotHeader = "ChIP_data_simu")
        PlotDots_dNorm(DataFile_ip = S1R$ChIP_results$LaLe_Data, SimuFile_ip = S1R$ChIP_results$LaLe_Simu, 
                       DataFile_in = S1R$Input_results$LaLe_Data, SimuFile_in = S1R$Input_results$LaLe_Simu,
                       PlotHeader = "ChIP_Input_data_simu")
        #BrDU
        PlotDots(DataFile = S1R$BrDU_results$LaLe_Data, PlotHeader = "BrDU_data")
        PlotDots_inNorm(DataFile_ip = S1R$BrDU_results$LaLe_Data, DataFile_in = S1R$Input_results$LaLe_Data, PlotHeader = "BrDU_Input_data")
        PlotDots_simNorm(DataFile = S1R$BrDU_results$LaLe_Data, SimuFile = S1R$BrDU_results$LaLe_Simu, PlotHeader = "BrDU_data_simu")
        PlotDots_dNorm(DataFile_ip = S1R$BrDU_results$LaLe_Data, SimuFile_ip = S1R$BrDU_results$LaLe_Simu, 
                       DataFile_in = S1R$Input_results$LaLe_Data, SimuFile_in = S1R$Input_results$LaLe_Simu,
                       PlotHeader = "BrDU_Input_data_simu")
        
        
        mtext("Lagging & Leading bias pattern at individual origins", side = 3, line = - 4, outer = TRUE, font = 2, cex = 2)
        mtext("Page 6", side = 1, line = - 2, outer = TRUE, font = 3, cex = 1)
        
        
        #seventh page
        
        #espan
        PlotDots(DataFile = S1R$eSPAN_results$LaLe_Data, PlotHeader = "eSPAN_data")
        PlotDots_inNorm(DataFile_ip = S1R$eSPAN_results$LaLe_Data, DataFile_in = S1R$BrDU_results$LaLe_Data, PlotHeader = "eSPAN_BrDU_data")
        PlotDots_simNorm(DataFile = S1R$eSPAN_results$LaLe_Data, SimuFile = S1R$eSPAN_results$LaLe_Simu, PlotHeader = "eSPAN_data_simu")
        PlotDots_dNorm(DataFile_ip = S1R$eSPAN_results$LaLe_Data, SimuFile_ip = S1R$eSPAN_results$LaLe_Simu, 
                       DataFile_in = S1R$BrDU_results$LaLe_Data, SimuFile_in = S1R$BrDU_results$LaLe_Simu,
                       PlotHeader = "eSPAN_BrDU_data_simu")
        
        #bsup
        PlotDots(DataFile = S1R$bSUP_results$LaLe_Data, PlotHeader = "bSUP_data")
        PlotDots_inNorm(DataFile_ip = S1R$bSUP_results$LaLe_Data, DataFile_in = S1R$Input_results$LaLe_Data, PlotHeader = "bSUP_Input_data")
        PlotDots_simNorm(DataFile = S1R$bSUP_results$LaLe_Data, SimuFile = S1R$bSUP_results$LaLe_Simu, PlotHeader = "bSUP_data_simu")
        PlotDots_dNorm(DataFile_ip = S1R$bSUP_results$LaLe_Data, SimuFile_ip = S1R$bSUP_results$LaLe_Simu, 
                       DataFile_in = S1R$Input_results$LaLe_Data, SimuFile_in = S1R$Input_results$LaLe_Simu,
                       PlotHeader = "bSUP_Input_data_simu")
        
        #esup
        PlotDots(DataFile = S1R$eSUP_results$LaLe_Data, PlotHeader = "eSUP_data")
        PlotDots_inNorm(DataFile_ip = S1R$eSUP_results$LaLe_Data, DataFile_in = S1R$bSUP_results$LaLe_Data, PlotHeader = "eSUP_bSUP_data")
        PlotDots_simNorm(DataFile = S1R$eSUP_results$LaLe_Data, SimuFile = S1R$eSUP_results$LaLe_Simu, PlotHeader = "eSUP_data_simu")
        PlotDots_dNorm(DataFile_ip = S1R$eSUP_results$LaLe_Data, SimuFile_ip = S1R$eSUP_results$LaLe_Simu, 
                       DataFile_in = S1R$bSUP_results$LaLe_Data, SimuFile_in = S1R$bSUP_results$LaLe_Simu,
                       PlotHeader = "eSUP_bSUP_data_simu")
        
        mtext("Lagging & Leading bias pattern at individual origins", side = 3, line = - 4, outer = TRUE, font = 2, cex = 2)
        mtext("Page 7", side = 1, line = - 2, outer = TRUE, font = 3, cex = 1)
        
        
      }
      
      ##boxplots
      
      page_8_and_9 <- {
        
        PlotBoxes <- function(DataFile, PlotHeader){
          
          bx1 <- DataFile$Lagg.sum; bx2 <- DataFile$Lead.sum
          bx1 <- log2((abs(bx1))^(sign(bx1))); bx2 <- log2((abs(bx2))^(sign(bx2)))
          
          bx1[!is.finite(bx1)] <- 0; bx2[!is.finite(bx2)] <- 0
          
          Y1 <- max(round(abs(range(c(bx1, bx2))) + 0.5))
          Y2 <- min(round(abs(range(c(bx1, bx2))) - 0.5))
          
          boxplot(bx1, bx2, outline = TRUE,
                  main = PlotHeader,
                  ylim = c(Y2, Y1),
                  names = c("lagging", "leading"), las = 1, 
                  col = c(adjustcolor("orange", alpha.f = 0.8, red.f = 1), adjustcolor("purple", alpha.f = 0.8, green.f = 1)), 
                  ylab = "log2 of read_erichment")
          abline(h = 0, lwd = 0.5, lty = 2, col = adjustcolor('black', alpha.f = 0.75))
          
          if(sum(bx1) == 0 || sum(bx2) == 0){
            p_value <- "NA"
          } else {
            p_value <- signif(wilcox.test(bx1, bx2, paired = TRUE)$p.value, 2)
          }
          mtext(side = 1, line = 1, cex = 0.6, paste0(" ", "\n", p_value))
          
        }
        PlotBoxes_inNorm <- function(DataFile_ip, DataFile_in, PlotHeader){
          
          lagg_ip <- DataFile_ip$Lagg.sum; lead_ip <- DataFile_ip$Lead.sum
          lagg_in <- DataFile_in$Lagg.sum; lead_in <- DataFile_in$Lead.sum
          
          bx1 <- lagg_ip/lagg_in; bx2 <- lead_ip/lead_in
          
          bx1 <- log2((abs(bx1))^(sign(bx1))); bx2 <- log2((abs(bx2))^(sign(bx2)))
          
          bx1[!is.finite(bx1)] <- 0; bx2[!is.finite(bx2)] <- 0
          #bx1[bx1 == 0] <- 1; bx2[bx2 == 0] <- 1
          Y <- max(round(abs(range(c(bx1, bx2))) + 0.5))
          
          boxplot(bx1, bx2, outline = TRUE,
                  main = PlotHeader,
                  ylim = c(-Y, Y),
                  names = c("lagging", "leading"), las = 1, 
                  col = c(adjustcolor("orange", alpha.f = 0.8, red.f = 1), adjustcolor("purple", alpha.f = 0.8, green.f = 1)), 
                  ylab = "log2 of read_erichment")
          
          abline(h = 0, lwd = 0.75, lty = 3, col = adjustcolor('black', alpha.f = 0.75))
          
          if(sum(bx1) == 0 || sum(bx2) == 0){
            p_value <- "NA"
          } else {
            p_value <- signif(wilcox.test(bx1, bx2, paired = TRUE)$p.value, 2)
          }
          mtext(side = 1, line = 1, cex = 0.6, paste0(" ", "\n", p_value))
          
        }
        PlotBoxes_simNorm <- function(DataFile, SimuFile, PlotHeader){
          
          lagg_data <- DataFile$Lagg.sum; lagg_simu <- SimuFile$lagging
          lagg_data_simu <- lagg_data/lagg_simu; lagg_data_simu[!is.finite(lagg_data_simu)] <- 0
          
          lead_data <- DataFile$Lead.sum; lead_simu <- SimuFile$leading
          lead_data_simu <- lead_data/lead_simu; lead_data_simu[!is.finite(lead_data_simu)] <- 0
          
          bx1 <- lagg_data_simu; bx2 <- lead_data_simu
          
          bx1 <- log2((abs(bx1))^(sign(bx1))); bx2 <- log2((abs(bx2))^(sign(bx2)))
          bx1[!is.finite(bx1)] <- 0; bx2[!is.finite(bx2)] <- 0
          #bx1[bx1 == 0] <- 1; bx2[bx2 == 0] <- 1
          
          Y <- max(round(abs(range(c(bx1, bx2))) + 0.5))
          
          boxplot(bx1, bx2, outline = TRUE,
                  main = PlotHeader,
                  ylim = c(-Y, Y),
                  names = c("lagging", "leading"), las = 1, 
                  col = c(adjustcolor("orange", alpha.f = 0.8, red.f = 1), adjustcolor("purple", alpha.f = 0.8, green.f = 1)), 
                  ylab = "log2 of read_erichment")
          
          
          abline(h = 0, lwd = 0.5, lty = 2, col = adjustcolor('black', alpha.f = 0.75))
          
          
          if(sum(bx1) == 0 || sum(bx2) == 0){
            p_value <- "NA"
          } else {
            p_value <- signif(wilcox.test(bx1, bx2, paired = TRUE)$p.value, 2)
          }
          mtext(side = 1, line = 1, cex = 0.6, paste0(" ", "\n", p_value))
        }
        PlotBoxes_dNorm <- function(DataFile_ip, SimuFile_ip, DataFile_in, SimuFile_in, PlotHeader){
          
          #ip
          lagg_data_ip <- DataFile_ip$Lagg.sum; lagg_simu_ip <- SimuFile_ip$lagging
          lagg_data_simu_ip <- lagg_data_ip/lagg_simu_ip; lagg_data_simu_ip[!is.finite(lagg_data_simu_ip)] <- 0
          
          lead_data_ip <- DataFile_ip$Lead.sum; lead_simu_ip <- SimuFile_ip$leading
          lead_data_simu_ip <- lead_data_ip/lead_simu_ip; lead_data_simu_ip[!is.finite(lead_data_simu_ip)] <- 0
          #in
          lagg_data_in <- DataFile_in$Lagg.sum; lagg_simu_in <- SimuFile_in$lagging
          lagg_data_simu_in <- lagg_data_in/lagg_simu_in; lagg_data_simu_in[!is.finite(lagg_data_simu_in)] <- 0
          
          lead_data_in <- DataFile_in$Lead.sum; lead_simu_in <- SimuFile_in$leading
          lead_data_simu_in <- lead_data_in/lead_simu_in; lead_data_simu_in[!is.finite(lead_data_simu_in)] <- 0
          #
          lagg_ip_in <- lagg_data_simu_ip/lagg_data_simu_in; lagg_ip_in[!is.finite(lagg_ip_in)] <- 0
          lead_ip_in <- lead_data_simu_ip/lead_data_simu_in; lead_ip_in[!is.finite(lead_ip_in)] <- 0
          
          #
          bx1 <- lagg_ip_in; bx2 <- lead_ip_in
          
          bx1 <- log2((abs(bx1))^(sign(bx1))); bx2 <- log2((abs(bx2))^(sign(bx2)))
          bx1[!is.finite(bx1)] <- 0; bx2[!is.finite(bx2)] <- 0
          #bx1[bx1 == 0] <- 1; bx2[bx2 == 0] <- 1
          
          
          Y <- max(round(abs(range(c(bx1, bx2))) + 0.5))
          
          boxplot(bx1, bx2, outline = TRUE,
                  main = PlotHeader,
                  ylim = c(-Y, Y),
                  names = c("lagging", "leading"), las = 1, 
                  col = c(adjustcolor("orange", alpha.f = 0.8, red.f = 1), adjustcolor("purple", alpha.f = 0.8, green.f = 1)), 
                  ylab = "log2 of read_erichment")
          
          
          abline(h = 0, lwd = 0.5, lty = 2, col = adjustcolor('black', alpha.f = 0.75))
          
          
          if(sum(bx1) == 0 || sum(bx2) == 0){
            p_value <- "NA"
          } else {
            p_value <- signif(wilcox.test(bx1, bx2, paired = TRUE)$p.value, 2)
          }
          mtext(side = 1, line = 1, cex = 0.6, paste0(" ", "\n", p_value))
        }
        
        #eigth page
        
        PlotMat <- {matrix(c(     
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # 
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16), 
          
          20,20,byrow=TRUE)}
        layout(PlotMat, c(1,1), c(1,1), TRUE)
        
        #input
        PlotBoxes(DataFile = S1R$Input_results$LaLe_Data, PlotHeader = "Input_data")
        PlotBoxes_inNorm(DataFile_ip = S1R$Input_results$LaLe_Data, DataFile_in = S1R$Input_results$LaLe_Data, PlotHeader = "Input_Input_data")
        PlotBoxes_simNorm(DataFile = S1R$Input_results$LaLe_Data, SimuFile = S1R$Input_results$LaLe_Simu, PlotHeader = "Input_data_simu")
        PlotBoxes_dNorm(DataFile_ip = S1R$Input_results$LaLe_Data, SimuFile_ip = S1R$Input_results$LaLe_Simu, 
                        DataFile_in = S1R$Input_results$LaLe_Data, SimuFile_in = S1R$Input_results$LaLe_Simu,
                        PlotHeader = "Input_Input_data_simu")
        #chip
        PlotBoxes(DataFile = S1R$ChIP_results$LaLe_Data, PlotHeader = "ChIP_data")
        PlotBoxes_inNorm(DataFile_ip = S1R$ChIP_results$LaLe_Data, DataFile_in = S1R$Input_results$LaLe_Data, PlotHeader = "ChIP_Input_data")
        PlotBoxes_simNorm(DataFile = S1R$ChIP_results$LaLe_Data, SimuFile = S1R$ChIP_results$LaLe_Simu, PlotHeader = "ChIP_data_simu")
        PlotBoxes_dNorm(DataFile_ip = S1R$ChIP_results$LaLe_Data, SimuFile_ip = S1R$ChIP_results$LaLe_Simu, 
                        DataFile_in = S1R$Input_results$LaLe_Data, SimuFile_in = S1R$Input_results$LaLe_Simu,
                        PlotHeader = "ChIP_Input_data_simu")
        #BrDU
        PlotBoxes(DataFile = S1R$BrDU_results$LaLe_Data, PlotHeader = "BrDU_data")
        PlotBoxes_inNorm(DataFile_ip = S1R$BrDU_results$LaLe_Data, DataFile_in = S1R$Input_results$LaLe_Data, PlotHeader = "BrDU_Input_data")
        PlotBoxes_simNorm(DataFile = S1R$BrDU_results$LaLe_Data, SimuFile = S1R$BrDU_results$LaLe_Simu, PlotHeader = "BrDU_data_simu")
        PlotBoxes_dNorm(DataFile_ip = S1R$BrDU_results$LaLe_Data, SimuFile_ip = S1R$BrDU_results$LaLe_Simu, 
                        DataFile_in = S1R$Input_results$LaLe_Data, SimuFile_in = S1R$Input_results$LaLe_Simu,
                        PlotHeader = "BrDU_Input_data_simu")
        
        
        mtext("Lagging & Leading read enrichment distributions", side = 3, line = - 4, outer = TRUE, font = 2, cex = 2)
        mtext("Page 8", side = 1, line = - 2, outer = TRUE, font = 3, cex = 1)
        
        
        #ninth page
        
        #espan
        PlotBoxes(DataFile = S1R$eSPAN_results$LaLe_Data, PlotHeader = "eSPAN_data")
        PlotBoxes_inNorm(DataFile_ip = S1R$eSPAN_results$LaLe_Data, DataFile_in = S1R$BrDU_results$LaLe_Data, PlotHeader = "eSPAN_BrDU_data")
        PlotBoxes_simNorm(DataFile = S1R$eSPAN_results$LaLe_Data, SimuFile = S1R$eSPAN_results$LaLe_Simu, PlotHeader = "eSPAN_data_simu")
        PlotBoxes_dNorm(DataFile_ip = S1R$eSPAN_results$LaLe_Data, SimuFile_ip = S1R$eSPAN_results$LaLe_Simu, 
                        DataFile_in = S1R$BrDU_results$LaLe_Data, SimuFile_in = S1R$BrDU_results$LaLe_Simu,
                        PlotHeader = "eSPAN_BrDU_data_simu")
        
        #bsup
        PlotBoxes(DataFile = S1R$bSUP_results$LaLe_Data, PlotHeader = "bSUP_data")
        PlotBoxes_inNorm(DataFile_ip = S1R$bSUP_results$LaLe_Data, DataFile_in = S1R$Input_results$LaLe_Data, PlotHeader = "bSUP_Input_data")
        PlotBoxes_simNorm(DataFile = S1R$bSUP_results$LaLe_Data, SimuFile = S1R$bSUP_results$LaLe_Simu, PlotHeader = "bSUP_data_simu")
        PlotBoxes_dNorm(DataFile_ip = S1R$bSUP_results$LaLe_Data, SimuFile_ip = S1R$bSUP_results$LaLe_Simu, 
                        DataFile_in = S1R$Input_results$LaLe_Data, SimuFile_in = S1R$Input_results$LaLe_Simu,
                        PlotHeader = "bSUP_Input_data_simu")
        
        #esup
        PlotBoxes(DataFile = S1R$eSUP_results$LaLe_Data, PlotHeader = "eSUP_data")
        PlotBoxes_inNorm(DataFile_ip = S1R$eSUP_results$LaLe_Data, DataFile_in = S1R$bSUP_results$LaLe_Data, PlotHeader = "eSUP_bSUP_data")
        PlotBoxes_simNorm(DataFile = S1R$eSUP_results$LaLe_Data, SimuFile = S1R$eSUP_results$LaLe_Simu, PlotHeader = "eSUP_data_simu")
        PlotBoxes_dNorm(DataFile_ip = S1R$eSUP_results$LaLe_Data, SimuFile_ip = S1R$eSUP_results$LaLe_Simu, 
                        DataFile_in = S1R$bSUP_results$LaLe_Data, SimuFile_in = S1R$bSUP_results$LaLe_Simu,
                        PlotHeader = "eSUP_bSUP_data_simu")
        
        mtext("Lagging & Leading read enrichment distributions", side = 3, line = - 4, outer = TRUE, font = 2, cex = 2)
        mtext("Page 9", side = 1, line = - 2, outer = TRUE, font = 3, cex = 1)
        
        
      }
      
      ##dotplots - parental and nascents
      
      page_10 <- {
        
        #DNA synthesis at origin
        Lagging_Synthesis_NP_Data <- data.frame(chrom = S1R$BrDU_results$LaLe_Data$chrom, name = S1R$BrDU_results$LaLe_Data$name,
                                                nascent = S1R$BrDU_results$LaLe_Data$Lagg.sum, parent = S1R$bSUP_results$LaLe_Data$Lagg.sum)
        Leading_Synthesis_NP_Data <- data.frame(chrom = S1R$BrDU_results$LaLe_Data$chrom, name = S1R$BrDU_results$LaLe_Data$name,
                                                nascent = S1R$BrDU_results$LaLe_Data$Lead.sum, parent = S1R$bSUP_results$LaLe_Data$Lead.sum)
        Lagging_Synthesis_NP_Simu <- data.frame(chrom = S1R$BrDU_results$LaLe_Simu$chrom,
                                                nascent = S1R$BrDU_results$LaLe_Simu$lagging, parent = S1R$bSUP_results$LaLe_Simu$lagging)
        Leading_Synthesis_NP_Simu <- data.frame(chrom = S1R$BrDU_results$LaLe_Simu$chrom, 
                                                nascent = S1R$BrDU_results$LaLe_Simu$leading, parent = S1R$bSUP_results$LaLe_Simu$leading)
        
        #Protein binding signal at origins between nascents and parents
        Lagging_Signal_NP_Data <- data.frame(chrom = S1R$eSPAN_results$LaLe_Data$chrom, name = S1R$eSPAN_results$LaLe_Data$name,
                                             nascent = S1R$eSPAN_results$LaLe_Data$Lagg.sum, parent = S1R$eSUP_results$LaLe_Data$Lagg.sum)
        Leading_Signal_NP_Data <- data.frame(chrom = S1R$eSPAN_results$LaLe_Data$chrom, name = S1R$eSPAN_results$LaLe_Data$name,
                                             nascent = S1R$eSPAN_results$LaLe_Data$Lead.sum, parent = S1R$eSUP_results$LaLe_Data$Lead.sum)
        Lagging_Signal_NP_Simu <- data.frame(chrom = S1R$eSPAN_results$LaLe_Simu$chrom,
                                             nascent = S1R$eSPAN_results$LaLe_Simu$lagging, parent = S1R$eSUP_results$LaLe_Simu$lagging)
        Leading_Signal_NP_Simu <- data.frame(chrom = S1R$eSPAN_results$LaLe_Simu$chrom, 
                                             nascent = S1R$eSPAN_results$LaLe_Simu$leading, parent = S1R$eSUP_results$LaLe_Simu$leading)
        
      
        #tenth page
        
        PlotMat <- {matrix(c(     
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,
          
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,
          
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          9,9,9,9,9,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,
          
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # 
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
          # 13,13,13,13,13,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16), 
          
          20,20,byrow=TRUE)}
        layout(PlotMat, c(1,1), c(1,1), TRUE)
        
        PlotDots_NPS <- function(DataFile, PlotHeader){
          
          nascent <- DataFile$nascent
          parent <- DataFile$parent
          nascent_over_parent <- nascent/parent
          nascent_over_parent[!is.finite(nascent_over_parent)] <- 0
          
          dotZ <- scale(nascent - parent) 
          dotP <- 2 * pnorm(-abs(dotZ))
          
          dotV <- nascent_over_parent
          dotV <- log2((abs(dotV))^(sign(dotV)))
          
          # global_z <- (mean(nascent) - mean(parent)) / sqrt(var(nascent)/length(nascent) + var(parent)/length(parent))
          # global_p <- 2 * pnorm(-abs(global_z))
          
          Y <- max(round(abs(range(dotV))+0.5))
          
          boxplot(dotV, ylim = c(-Y,Y), las = 2, ylab = 'log2 nascent/parent', 
                  cex.main=1, xlab = " ", bty = 'n', main=PlotHeader, col=adjustcolor("grey", alpha.f = 0.25), lwd=0.5)
          
          if(length(dotV[which(dotP > 10e-6)])>0){
            spreadPoints(values=dotV[which(dotP > 10e-6)], position=1.0, pointCex=0.65, col="blue", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          if(length(dotV[which(dotV >=0 & dotP <= 10e-6)])>0){
            spreadPoints(values=dotV[which(dotV >=0 & dotP <= 10e-6)], position=1.0, pointCex=0.65, col="#FF7F0E", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          if(length(dotV[which(dotV < 0 & dotP <= 10e-6)])>0){
            spreadPoints(values=dotV[which(dotV < 0 & dotP <= 10e-6)], position=1.0, pointCex=0.65, col="#6A0DAD", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          
          
          Nascent <- length(dotV[which(dotV >=0 & dotP <= 10e-6)])
          Parent <- length(dotV[which(dotV < 0 & dotP <= 10e-6)])
          Indeter <- length(dotV[which(dotP > 10e-6)])
          
          legend("bottomright", legend = c(paste0("nascent", " (", Nascent, ")"), 
                                           paste0("parent", " (", Parent, ")"), 
                                           paste0("indeter", " (", Indeter, ")")), 
                 col = c("#FF7F0E", "#6A0DAD", "blue"), pch = 19, pt.cex=0.65, bty = "n", cex = 0.65)
          
          
          InName <- PlotHeader
          
          #calculate p value
          #Decision Tree
          biasP <- binom.test(round(c(Nascent+Parent, Indeter)))$p.value
          
          if(biasP <= 10e-3 & (Nascent+Parent) > Indeter){
            
            biasQ <- binom.test(round(c(Nascent, Parent)))$p.value
            
            if(biasQ <= 10e-6){
              if(Nascent > Parent){
                conc <- paste0("Strong bias for Nascent strands", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
              if(Parent > Nascent){
                conc <- paste0("Strong bias for Parent templates", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-6 & biasQ <= 10e-4){
              if(Nascent > Parent){
                conc <- paste0("Weak bias for Nascent strands", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
              if(Parent > Nascent){
                conc <- paste0("Weak bias for Parent templates", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-4  & biasQ <= 10e-2){
              if(Nascent > Parent){
                conc <- paste0("Very weak bias for Nascent strands", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
              if(Parent > Nascent){
                conc <- paste0("Very weak bias for Parent templates", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-2){
              conc <- paste0("No significant preference for protein binding", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
            }
          } else {
            conc <- paste0("No decisive bias pattern detected", "\n","too many origins show indeterminable bias")
          }
          
          mtext(conc, side = 1, line = 2, cex = 0.65)
        }
        PlotDots_inNorm_NPS <- function(DataFile_ip, DataFile_in, PlotHeader){
          
          nascent_ip <- DataFile_ip$nascent; nascent_in <- DataFile_in$nascent
          nascent_ip_in <- nascent_ip/nascent_in; nascent_ip_in[!is.finite(nascent_ip_in)] <- 0
          
          parent_ip <- DataFile_ip$parent; parent_in <- DataFile_in$parent
          parent_ip_in <- parent_ip/parent_in; parent_ip_in[!is.finite(parent_ip_in)] <- 0
          
          dotZ <- (nascent_ip_in - parent_ip_in) / sqrt( (nascent_ip_in / nascent_in) + (parent_ip_in / parent_in) )
          dotP <- 2 * pnorm(-abs(dotZ))
          
          dotV <- nascent_ip_in/parent_ip_in; dotV <- log2((abs(dotV))^(sign(dotV)))
          
          Y <- max(round(abs(range(dotV))+0.5))
          
          boxplot(dotV, ylim = c(-Y,Y), las = 2, ylab = 'log2 nascent/parent', 
                  cex.main=1, xlab = " ", bty = 'n', main=PlotHeader, col=adjustcolor("grey", alpha.f = 0.25), lwd=0.5)
          
          if(length(dotV[which(dotP > 10e-6)])>0){
            spreadPoints(values=dotV[which(dotP > 10e-6)], position=1.0, pointCex=0.65, col="blue", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          if(length(dotV[which(dotV >=0 & dotP <= 10e-6)])>0){
            spreadPoints(values=dotV[which(dotV >=0 & dotP <= 10e-6)], position=1.0, pointCex=0.65, col="#FF7F0E", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          if(length(dotV[which(dotV < 0 & dotP <= 10e-6)])>0){
            spreadPoints(values=dotV[which(dotV < 0 & dotP <= 10e-6)], position=1.0, pointCex=0.65, col="#6A0DAD", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          
          
          Nascent <- length(dotV[which(dotV >=0 & dotP <= 10e-6)])
          Parent <- length(dotV[which(dotV < 0 & dotP <= 10e-6)])
          Indeter <- length(dotV[which(dotP > 10e-6)])
          
          legend("bottomright", legend = c(paste0("nascent", " (", Nascent, ")"), 
                                           paste0("parent", " (", Parent, ")"), 
                                           paste0("indeter", " (", Indeter, ")")), 
                 col = c("#FF7F0E", "#6A0DAD", "blue"), pch = 19, pt.cex=0.65, bty = "n", cex = 0.65)
          
          
          InName <- PlotHeader
          
          #calculate p value
          #Decision Tree
          biasP <- binom.test(round(c(Nascent+Parent, Indeter)))$p.value
          
          if(biasP <= 10e-3 & (Nascent+Parent) > Indeter){
            
            biasQ <- binom.test(round(c(Nascent, Parent)))$p.value
            
            if(biasQ <= 10e-6){
              if(Nascent > Parent){
                conc <- paste0("Strong bias for Nascent strands", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
              if(Parent > Nascent){
                conc <- paste0("Strong bias for Parent templates", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-6 & biasQ <= 10e-4){
              if(Nascent > Parent){
                conc <- paste0("Weak bias for Nascent strands", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
              if(Parent > Nascent){
                conc <- paste0("Weak bias for Parent templates", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-4  & biasQ <= 10e-2){
              if(Nascent > Parent){
                conc <- paste0("Very weak bias for Nascent strands", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
              if(Parent > Nascent){
                conc <- paste0("Very weak bias for Parent templates", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-2){
              conc <- paste0("No significant preference for protein binding", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
            }
          } else {
            conc <- paste0("No decisive bias pattern detected", "\n","too many origins show indeterminable bias")
          }
          
          mtext(conc, side = 1, line = 2, cex = 0.65)
        }
        PlotDots_simNorm_NPS <- function(DataFile, SimuFile, PlotHeader){
          
          nascent_data <- DataFile$nascent; nascent_simu <- SimuFile$nascent
          nascent_data_simu <- nascent_data/nascent_simu; nascent_data_simu[!is.finite(nascent_data_simu)] <- 0
          
          parent_data <- DataFile$parent; parent_simu <- SimuFile$parent
          parent_data_simu <- parent_data/parent_simu; parent_data_simu[!is.finite(parent_data_simu)] <- 0
          
          dotZ <- (nascent_data_simu - parent_data_simu) / sqrt( (nascent_data_simu / nascent_simu) + (parent_data_simu / parent_simu) )
          dotP <- 2 * pnorm(-abs(dotZ))
          
          dotV <- nascent_data_simu/parent_data_simu; dotV <- log2((abs(dotV))^(sign(dotV)))
          
          Y <- max(round(abs(range(dotV))+0.5))
          
          boxplot(dotV, ylim = c(-Y,Y), las = 2, ylab = 'log2 nascent/parent', 
                  cex.main=1, xlab = " ", bty = 'n', main=PlotHeader, col=adjustcolor("grey", alpha.f = 0.25), lwd=0.5)
          
          if(length(dotV[which(dotP > 10e-6)])>0){
            spreadPoints(values=dotV[which(dotP > 10e-6)], position=1.0, pointCex=0.65, col="blue", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          if(length(dotV[which(dotV >=0 & dotP <= 10e-6)])>0){
            spreadPoints(values=dotV[which(dotV >=0 & dotP <= 10e-6)], position=1.0, pointCex=0.65, col="#FF7F0E", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          if(length(dotV[which(dotV < 0 & dotP <= 10e-6)])>0){
            spreadPoints(values=dotV[which(dotV < 0 & dotP <= 10e-6)], position=1.0, pointCex=0.65, col="#6A0DAD", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          
          
          Nascent <- length(dotV[which(dotV >=0 & dotP <= 10e-6)])
          Parent <- length(dotV[which(dotV < 0 & dotP <= 10e-6)])
          Indeter <- length(dotV[which(dotP > 10e-6)])
          
          legend("bottomright", legend = c(paste0("nascent", " (", Nascent, ")"), 
                                           paste0("parent", " (", Parent, ")"), 
                                           paste0("indeter", " (", Indeter, ")")), 
                 col = c("#FF7F0E", "#6A0DAD", "blue"), pch = 19, pt.cex=0.65, bty = "n", cex = 0.65)
          
          
          InName <- PlotHeader
          
          #calculate p value
          #Decision Tree
          biasP <- binom.test(round(c(Nascent+Parent, Indeter)))$p.value
          
          if(biasP <= 10e-3 & (Nascent+Parent) > Indeter){
            
            biasQ <- binom.test(round(c(Nascent, Parent)))$p.value
            
            if(biasQ <= 10e-6){
              if(Nascent > Parent){
                conc <- paste0("Strong bias for Nascent strands", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
              if(Parent > Nascent){
                conc <- paste0("Strong bias for Parent templates", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-6 & biasQ <= 10e-4){
              if(Nascent > Parent){
                conc <- paste0("Weak bias for Nascent strands", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
              if(Parent > Nascent){
                conc <- paste0("Weak bias for Parent templates", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-4  & biasQ <= 10e-2){
              if(Nascent > Parent){
                conc <- paste0("Very weak bias for Nascent strands", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
              if(Parent > Nascent){
                conc <- paste0("Very weak bias for Parent templates", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-2){
              conc <- paste0("No significant preference for protein binding", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
            }
          } else {
            conc <- paste0("No decisive bias pattern detected", "\n","too many origins show indeterminable bias")
          }
          
          mtext(conc, side = 1, line = 2, cex = 0.65)
        }
        PlotDots_dNorm_NPS <- function(DataFile_ip, SimuFile_ip, DataFile_in, SimuFile_in, PlotHeader){
          
          #ip
          nascent_data_ip <- DataFile_ip$nascent; nascent_simu_ip <- SimuFile_ip$nascent
          nascent_data_simu_ip <- nascent_data_ip/nascent_simu_ip; nascent_data_simu_ip[!is.finite(nascent_data_simu_ip)] <- 0
          
          parent_data_ip <- DataFile_ip$parent; parent_simu_ip <- SimuFile_ip$parent
          parent_data_simu_ip <- parent_data_ip/parent_simu_ip; parent_data_simu_ip[!is.finite(parent_data_simu_ip)] <- 0
          #in
          nascent_data_in <- DataFile_in$nascent; nascent_simu_in <- SimuFile_in$nascent
          nascent_data_simu_in <- nascent_data_in/nascent_simu_in; nascent_data_simu_in[!is.finite(nascent_data_simu_in)] <- 0
          
          parent_data_in <- DataFile_in$parent; parent_simu_in <- SimuFile_in$parent
          parent_data_simu_in <- parent_data_in/parent_simu_in; parent_data_simu_in[!is.finite(parent_data_simu_in)] <- 0
          #
          nascent_ip_in <- nascent_data_simu_ip/nascent_data_simu_in; nascent_ip_in[!is.finite(nascent_ip_in)] <- 0
          parent_ip_in <- parent_data_simu_ip/parent_data_simu_in; parent_ip_in[!is.finite(parent_ip_in)] <- 0
          #
          dotZ <- (nascent_ip_in - parent_ip_in) / sqrt( (nascent_ip_in / nascent_simu_in) + (parent_ip_in / parent_simu_in) )
          dotP <- 2 * pnorm(-abs(dotZ))
          
          dotV <- nascent_ip_in/parent_ip_in; dotV <- log2((abs(dotV))^(sign(dotV)))
          
          Y <- max(round(abs(range(dotV))+0.5))
          
          boxplot(dotV, ylim = c(-Y,Y), las = 2, ylab = 'log2 nascent/parent', 
                  cex.main=1, xlab = " ", bty = 'n', main=PlotHeader, col=adjustcolor("grey", alpha.f = 0.25), lwd=0.5)
          
          if(length(dotV[which(dotP > 10e-6)])>0){
            spreadPoints(values=dotV[which(dotP > 10e-6)], position=1.0, pointCex=0.65, col="blue", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          if(length(dotV[which(dotV >=0 & dotP <= 10e-6)])>0){
            spreadPoints(values=dotV[which(dotV >=0 & dotP <= 10e-6)], position=1.0, pointCex=0.65, col="#FF7F0E", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          if(length(dotV[which(dotV < 0 & dotP <= 10e-6)])>0){
            spreadPoints(values=dotV[which(dotV < 0 & dotP <= 10e-6)], position=1.0, pointCex=0.65, col="#6A0DAD", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
          }
          
          
          Nascent <- length(dotV[which(dotV >=0 & dotP <= 10e-6)])
          Parent <- length(dotV[which(dotV < 0 & dotP <= 10e-6)])
          Indeter <- length(dotV[which(dotP > 10e-6)])
          
          legend("bottomright", legend = c(paste0("nascent", " (", Nascent, ")"), 
                                           paste0("parent", " (", Parent, ")"), 
                                           paste0("indeter", " (", Indeter, ")")), 
                 col = c("#FF7F0E", "#6A0DAD", "blue"), pch = 19, pt.cex=0.65, bty = "n", cex = 0.65)
          
          
          InName <- PlotHeader
          
          #calculate p value
          #Decision Tree
          biasP <- binom.test(round(c(Nascent+Parent, Indeter)))$p.value
          
          if(biasP <= 10e-3 & (Nascent+Parent) > Indeter){
            
            biasQ <- binom.test(round(c(Nascent, Parent)))$p.value
            
            if(biasQ <= 10e-6){
              if(Nascent > Parent){
                conc <- paste0("Strong bias for Nascent strands", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
              if(Parent > Nascent){
                conc <- paste0("Strong bias for Parent templates", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-6 & biasQ <= 10e-4){
              if(Nascent > Parent){
                conc <- paste0("Weak bias for Nascent strands", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
              if(Parent > Nascent){
                conc <- paste0("Weak bias for Parent templates", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-4  & biasQ <= 10e-2){
              if(Nascent > Parent){
                conc <- paste0("Very weak bias for Nascent strands", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
              if(Parent > Nascent){
                conc <- paste0("Very weak bias for Parent templates", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
              }
            }
            
            if(biasQ > 10e-2){
              conc <- paste0("No significant preference for protein binding", "\n","p-value (Nascent~Parent): ", signif(biasQ, 3))
            }
          } else {
            conc <- paste0("No decisive bias pattern detected", "\n","too many origins show indeterminable bias")
          }
          
          mtext(conc, side = 1, line = 2, cex = 0.65)
        }
        
        #Lagging_signal
        PlotDots_NPS(DataFile = Lagging_Signal_NP_Data, PlotHeader = "Lagg_NP_Signal_data")
        PlotDots_inNorm_NPS(DataFile_ip = Lagging_Signal_NP_Data, DataFile_in = Lagging_Synthesis_NP_Data, PlotHeader = "Lagg_NP_Signal_Synthesis_data")
        PlotDots_simNorm_NPS(DataFile = Lagging_Signal_NP_Data, SimuFile = Lagging_Signal_NP_Simu, PlotHeader = "Lagg_NP_Signal_data_simu")
        PlotDots_dNorm_NPS(DataFile_ip = Lagging_Signal_NP_Data, SimuFile_ip = Lagging_Signal_NP_Simu, 
                           DataFile_in = Lagging_Synthesis_NP_Data, SimuFile_in = Lagging_Synthesis_NP_Simu,
                           PlotHeader = "Lagg_NP_Signal_Synthesis_data_simu")
        
        #Leading_signal
        PlotDots_NPS(DataFile = Leading_Signal_NP_Data, PlotHeader = "Lead_NP_Signal_data")
        PlotDots_inNorm_NPS(DataFile_ip = Leading_Signal_NP_Data, DataFile_in = Leading_Synthesis_NP_Data, PlotHeader = "Lead_NP_Signal_Synthesis_data")
        PlotDots_simNorm_NPS(DataFile = Leading_Signal_NP_Data, SimuFile = Leading_Signal_NP_Simu, PlotHeader = "Lead_NP_Signal_data_simu")
        PlotDots_dNorm_NPS(DataFile_ip = Leading_Signal_NP_Data, SimuFile_ip = Leading_Signal_NP_Simu, 
                           DataFile_in = Leading_Synthesis_NP_Data, SimuFile_in = Leading_Synthesis_NP_Simu,
                           PlotHeader = "Lead_NP_Signal_Synthesis_data_simu")
        
        
        #####boxplots
    
        PlotBoxes_inNorm_NPS <- function(DataFile_ip, DataFile_in, PlotHeader){
          
          nascent_ip <- DataFile_ip$nascent; nascent_in <- DataFile_in$nascent
          nascent_ip_in <- nascent_ip/nascent_in; nascent_ip_in[!is.finite(nascent_ip_in)] <- 0
          
          parent_ip <- DataFile_ip$parent; parent_in <- DataFile_in$parent
          parent_ip_in <- parent_ip/parent_in; parent_ip_in[!is.finite(parent_ip_in)] <- 0
          
          bx1 <- nascent_ip_in; bx2 <- parent_ip_in
          
          bx1 <- log2((abs(bx1))^(sign(bx1))); bx2 <- log2((abs(bx2))^(sign(bx2)))
          
          bx1[!is.finite(bx1)] <- 0; bx2[!is.finite(bx2)] <- 0
          #bx1[bx1 == 0] <- 1; bx2[bx2 == 0] <- 1
          Y <- max(round(abs(range(c(bx1, bx2))) + 0.5))
          
          boxplot(bx1, bx2, outline = TRUE,
                  main = PlotHeader,
                  ylim = c(-Y, Y),
                  names = c("nascent", "parent"), las = 1, 
                  col = c(adjustcolor("#FF7F0E", alpha.f = 0.8, red.f = 1), adjustcolor("#6A0DAD", alpha.f = 0.8, green.f = 1)), 
                  ylab = "log2 of read_erichment")
          
          abline(h = 0, lwd = 0.75, lty = 3, col = adjustcolor('black', alpha.f = 0.75))
          
          if(sum(bx1) == 0 || sum(bx2) == 0){
            p_value <- "NA"
          } else {
            p_value <- signif(wilcox.test(bx1, bx2, paired = TRUE)$p.value, 2)
          }
          mtext(side = 1, line = 1, cex = 0.6, paste0(" ", "\n", p_value))
          
        }
        PlotBoxes_dNorm_NPS <- function(DataFile_ip, SimuFile_ip, DataFile_in, SimuFile_in, PlotHeader){
          
          #ip
          nascent_data_ip <- DataFile_ip$nascent; nascent_simu_ip <- SimuFile_ip$nascent
          nascent_data_simu_ip <- nascent_data_ip/nascent_simu_ip; nascent_data_simu_ip[!is.finite(nascent_data_simu_ip)] <- 0
          
          parent_data_ip <- DataFile_ip$parent; parent_simu_ip <- SimuFile_ip$parent
          parent_data_simu_ip <- parent_data_ip/parent_simu_ip; parent_data_simu_ip[!is.finite(parent_data_simu_ip)] <- 0
          #in
          nascent_data_in <- DataFile_in$nascent; nascent_simu_in <- SimuFile_in$nascent
          nascent_data_simu_in <- nascent_data_in/nascent_simu_in; nascent_data_simu_in[!is.finite(nascent_data_simu_in)] <- 0
          
          parent_data_in <- DataFile_in$parent; parent_simu_in <- SimuFile_in$parent
          parent_data_simu_in <- parent_data_in/parent_simu_in; parent_data_simu_in[!is.finite(parent_data_simu_in)] <- 0
          #
          nascent_ip_in <- nascent_data_simu_ip/nascent_data_simu_in; nascent_ip_in[!is.finite(nascent_ip_in)] <- 0
          parent_ip_in <- parent_data_simu_ip/parent_data_simu_in; parent_ip_in[!is.finite(parent_ip_in)] <- 0
          
          #
          bx1 <- nascent_ip_in; bx2 <- parent_ip_in
          
          bx1 <- log2((abs(bx1))^(sign(bx1))); bx2 <- log2((abs(bx2))^(sign(bx2)))
          bx1[!is.finite(bx1)] <- 0; bx2[!is.finite(bx2)] <- 0
          #bx1[bx1 == 0] <- 1; bx2[bx2 == 0] <- 1
          
          
          Y <- max(round(abs(range(c(bx1, bx2))) + 0.5))
          
          boxplot(bx1, bx2, outline = TRUE,
                  main = PlotHeader,
                  ylim = c(-Y, Y),
                  names = c("nascent", "parent"), las = 1, 
                  col = c(adjustcolor("#FF7F0E", alpha.f = 0.8, red.f = 1), adjustcolor("#6A0DAD", alpha.f = 0.8, green.f = 1)), 
                  ylab = "log2 of read_erichment")
          
          
          abline(h = 0, lwd = 0.5, lty = 2, col = adjustcolor('black', alpha.f = 0.75))
          
          
          if(sum(bx1) == 0 || sum(bx2) == 0){
            p_value <- "NA"
          } else {
            p_value <- signif(wilcox.test(bx1, bx2, paired = TRUE)$p.value, 2)
          }
          mtext(side = 1, line = 1, cex = 0.6, paste0(" ", "\n", p_value))
        }
        
        #
        PlotBoxes_inNorm_NPS(DataFile_ip = Lagging_Signal_NP_Data, DataFile_in = Lagging_Synthesis_NP_Data, PlotHeader = "Lagg_NP_Signal_Synthesis_data")
        PlotBoxes_dNorm_NPS(DataFile_ip = Lagging_Signal_NP_Data, SimuFile_ip = Lagging_Signal_NP_Simu, 
                            DataFile_in = Lagging_Synthesis_NP_Data, SimuFile_in = Lagging_Synthesis_NP_Simu,
                            PlotHeader = "Lagg_NP_Signal_Synthesis_data_simu")
        
        PlotBoxes_inNorm_NPS(DataFile_ip = Leading_Signal_NP_Data, DataFile_in = Leading_Synthesis_NP_Data, PlotHeader = "Lead_NP_Signal_Synthesis_data")
        PlotBoxes_dNorm_NPS(DataFile_ip = Leading_Signal_NP_Data, SimuFile_ip = Leading_Signal_NP_Simu, 
                           DataFile_in = Leading_Synthesis_NP_Data, SimuFile_in = Leading_Synthesis_NP_Simu,
                           PlotHeader = "Lead_NP_Signal_Synthesis_data_simu")
        
        #
        
        mtext("Nascent signal over parent at origins (eSPAN/eSUP)", side = 3, line = - 4, outer = TRUE, font = 2, cex = 2)
        mtext("Page 10", side = 1, line = - 2, outer = TRUE, font = 3, cex = 1)
        
        
      }
      
      dev.off()
      
      message(paste0("✅ Plot saved as ", basename(Sample_1), "_", "SUPeSPAN.pdf"))
    }
    
    for(i in 1:NumOfSamples){
      plot_results(Sample = paste0("S", i))
    }
    
    
    rm(list=ls())
    gc()
    
  }
  
  # run the analysis for a single sample with six different experiments (Input, ChIP, BrDU, eSPAN, eSUP, bSUP)
  
  Alignment_Mapping_CoverageRatio(Input_R1, Input_R2, BrDU_R1, BrDU_R2, ChIP_R1, ChIP_R2,
                                  eSPAN_R1, eSPAN_R2, bSUP_R1, bSUP_R2, eSUP_R1, eSUP_R2,
                                  ExpTitle)
  
  Simulation_Calculation_Plotting( Sample_1 = paste0("/Users/mohammed/Desktop/", ExpTitle) )
  
}
