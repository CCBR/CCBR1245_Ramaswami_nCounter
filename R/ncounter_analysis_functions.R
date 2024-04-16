# Functions for processing nCounter RNA expression data


# A function for processing RCC files into dataframes for the counts and metadata
# Assumes the sample IDs are in a field of the annotatrion named SampleID
process_RCC_files <- function(rcc.files.path, annotation.df){
  
  # Gather the names of RCC files
  files.RCC <- list.files(rcc.files.path)
  files.RCC = files.RCC[grepl('RCC',files.RCC)]
  
  # Gather number of genes and samples
  rcc.example.file <- file.path(rcc.files.path, files.RCC[1])
  genes.count = nrow(readRcc(rcc.example.file)$Code_Summary)
  ncol = length(files.RCC)
  
  # Set up the dataframe to hold all sample and gene data
  raw_expression = as.data.frame(matrix(nrow = genes.count,ncol = length(files.RCC)+2))
  colnames(raw_expression)[1:2] = c('Gene','Class')
  
  # Set up a dataframe to hold metadata and QC
  pData = as.data.frame(matrix(nrow = length(files.RCC),ncol = 11))
  colnames(pData) = c('BCAC_ID','SampleID','Owner','Comments','Date','GeneRLF','SystemAPF','imagingQC',
                      'bindingDensityQC','limitOfDetectionQC','positiveLinearityQC')
  
  # Populate gene and class columns of expression df
  raw_expression[,1:2] = readRcc(rcc.example.file)$Code_Summary[,c(2,1)]
  
  print("RCC files to be processed:")
  # Populate expression df, run QC, and record flags for each sample
  for (i in 1:length(files.RCC)){
    
    # Read in the RCC file
    print(files.RCC[i])
    rcc.path <- file.path(rcc.files.path, files.RCC[i])
    rcc = readRcc(rcc.path)
    raw = rcc$Code_Summary
    
    # Create flags for the QC thresholds
    raw_expression[,i+2] = as.numeric(raw$Count)
    colnames(raw_expression)[i+2] = strsplit(files.RCC[i],'_')[[1]][3]
    pData[i,2:7] = as.vector(rcc$Sample_Attributes)
    pData$imagingQC[i] = imagingQC(rcc)
    pData$imaging[i] <-  as.numeric(rcc$Lane_Attributes[3])
    pData$bindingDensityQC[i] = bindingDensityQC(rcc,.05,2.25)
    pData$bindingDensity[i] <- as.numeric(rcc$Lane_Attributes[6])
    pData$limitOfDetectionQC[i] = limitOfDetectionQC(rcc)
    pData$limitOfDetection[i] <- limitOfDetectionScore(rcc)
    pData$positiveLinearityQC[i] = positiveLinQC(rcc)
    pData$positiveLinearityQC[i] <- positiveLinScore(rcc)
  }
  
  # Correct the sample ID column name and sample name IDs
  annotation.df$"SampleID" <- gsub("-", " ", annotation.df$"SampleID")
  annotation.df$"SampleID" <- gsub("\\sSB", "SB", annotation.df$"SampleID")
  
  # Correct SampleID for pdata
  pData$"SampleID" <- gsub("-", " ", pData$"SampleID")
  pData$"SampleID" <- gsub("\\sSB", "SB", pData$"SampleID")
  
  # Combine all metadata into a single dataframe
  merged_pData <- merge(pData, annotation.df, by = "SampleID")
  
  
  #merged_pData[[contrast]] <- factor(merged_pData[[contrast]], levels = contrast_levels)
  
  # Create a feature df
  raw = raw_expression[,-c(1:2)]
  fData = raw_expression[,c(1:2)]
  rownames(raw) = fData$Gene
  
  # Correct the Sample ID column names for preceding spaces and hyphens
  colnames(raw) <- gsub("^\\s+", "", colnames(raw))
  colnames(raw) <- gsub("-", " ", colnames(raw))
  
  # List the information for the feature data frame
  str(fData)
  count_endog <- sum(fData$Class == "Endogenous")
  count_neg <- sum(fData$Class == "Negative")
  count_pos <- sum(fData$Class == "Positive")
  count_hk <- sum(fData$Class == "Housekeeping")
  
  print(paste0("Endogenous genes: ",count_endog))
  print(paste0("Negative genes: ",count_neg))
  print(paste0("Positive genes: ",count_pos))
  print(paste0("Housekeeping geness: ",count_hk))
  
  # List the housekeeping genes and label in the fData df
  hk.gene.list <- fData$Gene[fData$Class == "Housekeeping"]
  
  # Identify missing housekeeping genes in any samples
  merged_pData$HK_Gene_Miss = colSums(raw[hk.gene.list,] == 0)
  
  # Estbalish all of the rownames of the three dfs
  rownames(fData) = fData$Gene
  rownames(raw) = fData$Gene
  rownames(merged_pData) <- merged_pData$SampleID
  
  # Ensure the order of the SampleID columns in the expression data frame match the order of the rows in the metadata
  raw_reorder <- raw[, rownames(merged_pData)]
  
  print("Initial data frames")
  
  return(list("raw.counts" = raw_reorder, 
              "pData" = merged_pData, 
              "fData" = fData, 
              "hk.gene.list" = hk.gene.list))
  
}

run_RUV_DESEQ2 <- function(contrast, 
                           contrast.levels, 
                           annotation.data, 
                           raw.counts, 
                           feature.data, 
                           k = 1, 
                           additional.covar = NULL, 
                           hk.genes = NULL, 
                           exclude.hk.genes = NULL){
  
  # Remove samples with NA from the read counts dataframe
  NA_rows <- which(!complete.cases(annotation.data[[contrast]]))
  NA_samples <- rownames(annotation.data[NA_rows, ])
  
  if(length(NA_rows) > 0){
    raw.counts <- raw.counts[, -which(names(raw.counts) %in% NA_samples)]
  }
  
  # Remove sample with NA from the annotation
  if(length(NA_rows) > 0){
    annotation.data <- annotation.data[complete.cases(annotation.data[[contrast]]), ]
  }
  
  # Run RUVseq with k sources of variance
  set = run_RUVseq(raw.counts = raw.counts, 
                annotation.data = annotation.data, 
                feature.data = feature.data, 
                k = k, 
                hk.genes = hk.genes, 
                exclude.hk.genes = exclude.hk.genes)$set
  
  # Remove row of the contrast that contain NA values
  set <- set[complete.cases(pData(set)[[contrast]]), ]
  
  # Should date of the experiment be used as a covariate? 
  if(is.null(additional.covar)){
    dds <- DESeqDataSetFromMatrix(countData = counts(set)[1:nrow(set),],
                                  colData = pData(set),
                                  design = as.formula(paste("~ W_1 + ", 
                                                            contrast)))
  } else{
    dds <- DESeqDataSetFromMatrix(countData = counts(set)[1:nrow(set),],
                                  colData = pData(set),
                                  design = as.formula(paste("~ W_1 +", 
                                                            additional.covar, 
                                                            "+", 
                                                            contrast)))
  }
  
  
  # Run DEseq
  dds <- DESeq(dds)
  
  # Generate contrasts per group
  contrast.ruv <- as.data.frame(results(dds,contrast = c(contrast, 
                                                         contrast.levels)))
  
  # Set up final results df
  contrast.ruv$normalization <- "RUVseq"
  contrast.ruv$gene <- rownames(contrast.ruv)
  
  # Reformat the table of DEGs to be exported
  contrast.ruv.deg.list <- contrast.ruv[contrast.ruv$pvalue < 0.05, ]
  rownames(contrast.ruv.deg.list) <- NULL
  contrast.ruv.deg.list <- contrast.ruv.deg.list[, c("gene", 
                                                     setdiff(
                                                       names(contrast.ruv.deg.list), 
                                                       "gene"))]
  
  
  return(list("dds" = dds, 
              "all.de.results" = contrast.ruv, 
              "deg.list" = contrast.ruv.deg.list
              ))
  
}

run_RUVseq <- function(raw.counts, 
                       annotation.data, 
                       feature.data, 
                       k, 
                       hk.genes = NULL, 
                       exclude.hk.genes = NULL){
  
  ### INPUT: raw.counts - p x n raw expressions with p genes and n samples
  ###        annotation.data - phenotype metadata across samples
  ###        feature.data - feature metadata across genes
  ###        k - number of dimensions of unwanted variation estimated
  ###        exclude.hk.genes - vector of gene names to exclude
  

  
  # Subset the feature names and raw counts df for the same features
  feature.data = feature.data[rownames(raw.counts),]
  int = intersect(rownames(raw.counts), rownames(feature.data))
  feature.data = feature.data[int,]
  raw.counts = raw.counts[int,]
  
  # Create a Expression Set object
  set <- newSeqExpressionSet(as.matrix(round(raw.counts)),
                             phenoData=annotation.data,
                             featureData=feature.data)
  
  # Label housekeeping if they were not provided
  if (!is.null(hk.genes)){
    
    fData(set)$Class[rownames(set) %in% hk.genes] = 'Housekeeping'
    
  }
  
  # Define the housekeeping features
  hk.genes <- rownames(set)[fData(set)$Class == "Housekeeping"]
  hk.genes = hk.genes[!(hk.genes %in% exclude.hk.genes)]
  
  # Run normalization 
  set <- betweenLaneNormalization(set, which="upper")
  set <- RUVg(set, hk.genes, k=k)
  
  # Create a DESEQ2 data set
  dds <- DESeqDataSetFromMatrix(counts(set), 
                                colData=pData(set), 
                                design=~1)
  
  rowData(dds) <- feature.data
  
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds)
  cts <- counts(dds, normalized=TRUE)
  disp <- pmax((rowVars(cts, useNames = TRUE) - rowMeans(cts)),0)/rowMeans(cts)^2
  mcols(dds)$dispGeneEst <- disp
  dds <- estimateDispersionsFit(dds, fitType="mean")
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  mat <- assay(vsd)
  covars <- as.matrix(colData(dds)[,grep("W",colnames(colData(dds))),drop=FALSE])
  mat <- removeBatchEffect(mat, covariates=covars)
  assay(vsd) <- mat
  return(list(set = set,vsd = vsd))
  
  
  
}

imagingQC <- function(rcc){
  
  
  #### INPUT: rcc - input from rcc
  #### OUTPUT: flag for imaging quality
  
  fovRatio = as.numeric(rcc$Lane_Attributes[3]) / as.numeric(rcc$Lane_Attributes[2])
  if (!(fovRatio > .75)) {return('Flag')}
  if (fovRatio > .75) {return('No flag')}
  
}

bindingDensityQC <- function(rcc,low,high){
  
  
  #### INPUT: rcc - input from rcc
  ####         low, high - the lower and upper limits for binding density
  #### OUTPUT: flag for binding density
  
  bd = as.numeric(rcc$Lane_Attributes[6])
  if(!(bd < high & bd > low)) {return('Flag')}
  if (bd < high & bd > low) {return('No flag')}
  
  
}

limitOfDetectionQC <- function(rcc,numSD = 0){
  
  #### INPUT: rcc - input from rcc
  ####         numSD - number of standard deviations to calibrate the LOD
  #### OUTPUT: flag for limit of detection
  
  counts = rcc$Code_Summary
  posE = as.numeric(counts$Count[counts$Name == 'POS_E'])
  negControls = as.numeric(counts$Count[grepl('NEG',counts$Name)])
  if(!(posE > mean(negControls) + numSD*sd(negControls))) {return('Flag')}
  if (posE > mean(negControls) + numSD*sd(negControls)) {return('No flag')}
  
}

limitOfDetectionScore <- function(rcc,numSD = 0){
  
  #### INPUT: rcc - input from rcc
  ####         numSD - number of standard deviations to calibrate the LOD
  #### OUTPUT: the score for limit of detection
  
  counts = rcc$Code_Summary
  
  # Gather the counts gfor positive and neagtive probes
  posE = as.numeric(counts$Count[counts$Name == 'POS_E'])
  negControls = as.numeric(counts$Count[grepl('NEG',counts$Name)])
  
  # Calculate the negative background
  neg_background <- mean(negControls) + numSD*sd(negControls)
  
  # Return the limit of detection comparison
  return(paste0(posE, "/", neg_background))
}

positiveLinQC <- function(rcc){
  
  #### INPUT: rcc - input from rcc
  #### OUTPUT: flag for linearity for positive controls
  
  
  counts = rcc$Code_Summary
  posControls = as.numeric(counts$Count[grepl('POS_',counts$Name)])
  known = c(128,128/4,128/16,128/64,128/256,128/(256*4))
  r2 = summary(lm(sort(posControls)~sort(known)))$r.squared
  if(!(r2 > .95) | is.na(r2)) {return('Flag')}
  if(r2 > .95) {return('No flag')}
  
}

positiveLinScore <- function(rcc){
  
  #### INPUT: rcc - input from rcc
  #### OUTPUT: score for linearity for positive controls
  
  
  counts = rcc$Code_Summary
  posControls = as.numeric(counts$Count[grepl('POS_',counts$Name)])
  known = c(128,128/4,128/16,128/64,128/256,128/(256*4))
  r2 = summary(lm(sort(posControls)~sort(known)))$r.squared
  return(r2)
  
}

makeRLEplot <- function(data,metadata,id){
  
  #### INPUT: data - matrix of expressions with genes on rows and samples on columns
  ####        metadata - matrix of metadata with a column that corresponds to the colnames of data
  ####        id - colname of sample ids
  #### OUTPUT: ggplot2 RLE plot
  
  data = data - apply(data,1,median)
  stack = stack(data)
  colnames(stack)[1] = id
  stackPlot = merge(stack,metadata,by=id)
  colnames(stackPlot)[1:2] = c('Sample','values')
  rle_plots = ggplot(data = stackPlot,aes(x = Sample,y = values, color = ER_status)) +
    geom_boxplot(coef = 6) + theme_minimal() +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=24),
          plot.title = element_text(size = 30),
          legend.title=element_text(size=20),
          legend.text=element_text(size=20),
          strip.text = element_text(size=24),
          panel.spacing=unit(1, "lines"),
          panel.border = element_rect(color = "grey", 
                                      fill = NA, size = .1),
          legend.position = 'bottom',
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + xlab('Sample') +
    ylab('Median deviation of log expression') + ylim(c(-4,4))
  return(rle_plots)
  
}
