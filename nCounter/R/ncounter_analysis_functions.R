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
  
  # Establish all of the rownames of the three dfs
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
                           hk.genes){
  
  # Filter annotation for specific contrast levels
  annotation.data <- annotation.data %>% 
    filter(.data[[contrast]] %in% contrast.levels)
  
  annotation.samples <- annotation.data$SampleID
  
  rownames(annotation.data) <- annotation.samples
  
  # Filter the raw counts for the same samples
  raw.counts <- raw.counts %>% 
    select(all_of(annotation.samples))
  
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
                hk.genes = hk.genes)$set
  
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
                       hk.genes){
  
  ### INPUT: raw.counts - p x n raw expressions with p genes and n samples
  ###        annotation.data - phenotype metadata across samples
  ###        feature.data - feature metadata across genes
  ###        k - number of dimensions of unwanted variation estimated
  ###        hk.genes - vector of gene names to use as housekeeping
  

  
  # Subset the feature names and raw counts df for the same features
  feature.data = feature.data[rownames(raw.counts),]
  int = intersect(rownames(raw.counts), rownames(feature.data))
  feature.data = feature.data[int,]
  raw.counts = raw.counts[int,]
  
  # Create a Expression Set object
  set <- newSeqExpressionSet(as.matrix(round(raw.counts)),
                             phenoData=annotation.data,
                             featureData=feature.data)
  
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

create_GSEA_preranked <- function(counts, 
                                  annotation, 
                                  contrast, 
                                  contrast.levels){
  
  # Gather the signal to noise ratio for GSEA ranking
  # Default method for ranking genes from GSEA manual:
  # https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Metrics_for_Ranking
  
  # Contrast level A is the "condition" (positive when calculating fold change)
  rownames(annotation) <- annotation$SampleID
  contrast.A.annotation <- annotation %>% 
    filter(!!sym(contrast) == contrast.levels[1])
  
  contrast.A.sampleIDs <- rownames(contrast.A.annotation)
  
  contrast.A.counts <- as.data.frame(counts) %>% 
    select(all_of(contrast.A.sampleIDs))
  
  contrast.A.counts$gene <- rownames(contrast.A.counts)
  
  # Contrast level B is the "reference" (negative when calculating fold change)
  
  contrast.B.annotation <- annotation %>% 
    filter(!!sym(contrast) == contrast.levels[2])
  
  contrast.B.sampleIDs <- rownames(contrast.B.annotation)
  
  contrast.B.counts <- as.data.frame(counts) %>% 
    select(all_of(contrast.B.sampleIDs))
  
  contrast.B.counts$gene <- rownames(contrast.B.counts)
  
  # Add a column to each contrast level for the mean and standard deviation
  contrast.A.counts <- contrast.A.counts %>% 
    mutate(mean.A = rowMeans(select_if(., is.numeric))) %>%  
    mutate(stdev.A = apply(select_if(., is.numeric), 1, sd))
  
  contrast.B.counts <- contrast.B.counts %>% 
    mutate(mean.B = rowMeans(select_if(., is.numeric))) %>%  
    mutate(stdev.B = apply(select_if(., is.numeric), 1, sd))
  
  GSEA.preanked.df <- merge(contrast.A.counts, contrast.B.counts, by = "gene")
  
  GSEA.preanked.df <- GSEA.preanked.df %>% 
    mutate(signal2noise = (mean.A - mean.B)/(stdev.A + stdev.B)) %>% 
    arrange(desc(signal2noise)) %>% 
    select(c(gene, mean.A, mean.B, stdev.A, stdev.B, signal2noise))
  
  return(GSEA.preanked.df)
  
}

check_hk_for_de <- function(hk.genes, 
                            raw.counts, 
                            annotation, 
                            contrast.field){
  
  # Set up HK genes to check if they are DE
  hk.raw = raw.counts[hk.genes,]
  pval = vector(length = nrow(hk.raw))
  hk.de <- data.frame(gene = hk.genes)
  
  # Ensure the sample names are ordered the same for both read counts and annotation
  annotation <- annotation[match(colnames(hk.raw), annotation$SampleID), ]
  
  # Check HK genes for DE for the current contrast
  for (i in 1:nrow(hk.raw)){
    reg = glm.nb(as.numeric(hk.raw[i,]) ~ as.factor(annotation[[contrast.field]]))
    pval[i] = coef(summary(reg))[2,4]
  }
  
  # Create a df for the de results of HK comparisons
  hk.de <- hk.de %>% 
    mutate(pval = pval)
  
  # Track HK genes that are DEGs
  hk.degs <- c()
  for(gene in hk.de$gene){
    
    pval.gene <- hk.de$pval[hk.de$gene == gene]
    
    if(pval.gene < 0.05){
      hk.degs <- c(gene, hk.degs)
    }
  }
  
  # Remove HK that are DE from HK list 
  if(length(hk.degs > 0)){
    
    print(paste0("HK genes that are DE in ", contrast.field, ":", hk.degs))
    print("Removing HK DEGs from HK list")
    hk.genes.final <- hk.genes[!hk.genes %in% hk.degs]
    
  } else {
    
    print(paste0("No HK genes are DE for ", contrast.field))
    hk.genes.final <- hk.genes
    
  }
  
  return(list("hk.degs" = hk.degs, 
              "hk.genes.final" = hk.genes.final))
  
}

gene_counts_violin_boxplot <- function(counts, 
                                annotation.df, 
                                gene.list, 
                                annotation.field, 
                                display.summary.stat = FALSE, 
                                compare.groups = FALSE){
  
  # Set up the annotation df
  annotation.fields <- c("SampleID", annotation.field)
  subset.annotation <- as.data.frame(annotation.df)[, annotation.fields, 
                                                 drop = FALSE]
  
  # Check if goi are found in counts
  for(gene in gene.list){
    
    if(!(gene %in% counts$gene)){
      
      print(paste0(gene, " not found in counts file"))
      
      gene.list <- gene.list[-which(gene.list == gene)]
      
    }
    
  }
  
  # Convert gene counts to log2
  gene.counts <- counts %>% 
    filter(gene %in% gene.list) %>% 
    mutate(across(where(is.numeric), ~.+1)) %>% 
    mutate(across(where(is.numeric), log2))
  
  # Set up counts for merge with annotation
  gene.counts.transpose <- as.data.frame(t(gene.counts)) %>% 
    rownames_to_column(var = "SampleID")
  
  # Create master annotation/counts df
  counts.anno.df <- merge(gene.counts.transpose, 
                          subset.annotation, 
                          by = "SampleID")
  
  # Set up the annotation/counts df for ggplot2
  counts.anno.df.melt <- counts.anno.df %>% 
    pivot_longer(cols = all_of(gene.list), 
                 names_to = "gene", 
                 values_to = "log_counts")
  
  counts.anno.df.melt$log_counts <- as.numeric(counts.anno.df.melt$log_counts)
  
  # Create a combined boxplot and violin plot
  if(display.summary.stat == TRUE){
    
    field.violin.plot <- ggplot(counts.anno.df.melt, aes(x = !!sym(annotation.field), 
                                                         y = log_counts, 
                                                         fill = !!sym(annotation.field))) +
      geom_violin() + 
      geom_boxplot(width = 0.2, fill = "white") + 
      facet_wrap(~ gene) + 
      labs(x = paste(gsub("_", " ", annotation.field)), 
           y = "Log2 Counts", 
           title = paste0("Log2 Counts for ", gsub("_", " ", annotation.field))) + 
      stat_summary(
        fun = mean, 
        geom = "text", 
        aes(label = paste("mean:", round(after_stat(y), 2))), 
        vjust = -0.5, 
        color = "darkblue"
      ) + 
      stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "darkblue", 
                   fatten = 0.5) +
      stat_summary(
        fun.min = function(x) quantile(x, 0.25),
        fun = function(x) quantile(x, 0.25), 
        geom = "text", 
        aes(label = paste("Q1:", round(after_stat(y), 2))),
        vjust = 1.5, 
        color = "black"
      ) +
      stat_summary(
        fun.max = function(x) quantile(x, 0.75),
        fun = function(x) quantile(x, 0.75), 
        geom = "text", 
        aes(label = paste("Q3:", round(after_stat(y), 2))),
        vjust = -1.5, 
        color = "black"
      ) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  } else if(compare.groups == TRUE) {
    
    field.violin.plot <- ggplot(counts.anno.df.melt, aes(x = !!sym(annotation.field), 
                                                         y = log_counts, 
                                                         fill = !!sym(annotation.field))) +
      geom_violin() + 
      geom_boxplot(width = 0.2, fill = "white") + 
      stat_compare_means(method = "wilcox.test", label.y = 0.5) + 
      facet_wrap(~ gene) + 
      labs(x = paste(gsub("_", " ", annotation.field)), 
           y = "Log2 Counts", 
           title = paste0("Log2 Counts for ", gsub("_", " ", annotation.field))) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else(
    
    field.violin.plot <- ggplot(counts.anno.df.melt, aes(x = !!sym(annotation.field), 
                                                         y = log_counts, 
                                                         fill = !!sym(annotation.field))) +
      geom_violin() + 
      geom_boxplot(width = 0.2, fill = "white") + 
      facet_wrap(~ gene) + 
      labs(x = paste(gsub("_", " ", annotation.field)), 
           y = "Log2 Counts", 
           title = paste0("Log2 Counts for ", gsub("_", " ", annotation.field))) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  )
  
  
 
  
  return(field.violin.plot)
    
}

make_volcano <- function(de.results, 
                         title, 
                         legend.title, 
                         x.axis.title, 
                         fc.limit = 1, 
                         custom.gene.labels = NULL, 
                         remove.controls = TRUE, 
                         remove.genes){ 
  
  # Remove controls if applicable
  if(remove.controls == TRUE){
    
    NEG.indices <- grep("NEG_", rownames(de.results))
    POS.indices <- grep("POS_", rownames(de.results))
    control.indices <- c(NEG.indices, POS.indices)
    
    # Remove the control probes
    de.results <- de.results[-control.indices,]
    
  } 
  
  # Remove custom gene input
  if(!is.null(remove.genes)){
    
    de.results <- de.results[!rownames(de.results) %in% remove.genes, ]
    
  }
  

  
  # Create a column for direction of DEGs
  de.results$de_direction <- "NONE"
  de.results$de_direction[de.results$padj < 0.05 & 
                             de.results$log2FoldChange > fc.limit] <- "UP"
  de.results$de_direction[de.results$padj < 0.05 & 
                             de.results$log2FoldChange < -fc.limit] <- "DOWN"
  
  # Create a label for DEGs
  if(is.null(custom.gene.labels)){
    
    de.results$deglabel <- ifelse(de.results$de_direction == "NONE", 
                                   NA, 
                                   de.results$gene)
    
  } else {
    
    de.results$deglabel <- ifelse(de.results$gene %in% custom.gene.labels, 
                                   de.results$gene, 
                                   NA)
    
  }
  
  
  
  # Compute the scale for the volcano x-axis
  log2.scale <- max(abs(de.results$log2FoldChange))
  
  # Establish the color scheme for the volcano plot
  if(is.null(custom.gene.labels)){
    
    contrast.level.colors <- c("steelblue4", "grey", "violetred4")
    names(contrast.level.colors) <- c("DOWN", "NONE", "UP")
    
    volcano.plot <- ggplot(data = de.results, aes(x = log2FoldChange, 
                                                   y = -log10(padj), 
                                                   col = de_direction, 
                                                   label = deglabel)) +
      geom_vline(xintercept = c(-fc.limit, fc.limit), col = "gray", linetype = 'dashed') +
      geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
      xlim(-7.5, 7.5) + 
      labs(x = x.axis.title,
           y = "-log10 adjusted p-value", 
           title = title) + 
      geom_point(size = 2) +
      scale_color_manual(legend.title, 
                         values = contrast.level.colors) + 
      geom_text_repel(max.overlaps = Inf) + 
      xlim(-log2.scale-1, log2.scale+1) + 
      theme(plot.title = element_text(hjust = 0.5))
    
  } else {
    
    # Label the custom genes depending on significance
    de.results <- de.results %>% 
      mutate(custom.label = ifelse(!is.na(deglabel) & de_direction == "NONE", 
                                   "BLACK", 
                                   ifelse(!is.na(deglabel) & de_direction != "NONE", 
                                          de_direction, 
                                          "NONE")))
    
    contrast.level.colors <- c("steelblue4", "grey", "violetred4", "black")
    names(contrast.level.colors) <- c("DOWN", "NONE", "UP", "BLACK")
    
    de.results.labeled <- de.results %>%
      filter(custom.label != "NONE")
    
    de.results.unlabeled <- de.results %>% 
      filter(custom.label == "NONE")
    
    
    volcano.plot <- ggplot() + 
      geom_point(data = de.results.unlabeled, aes(x = log2FoldChange, 
                                                   y = -log10(padj), 
                                                   col = custom.label, 
                                                   alpha = 0.5)) + 
      geom_point(data = de.results.labeled, aes(x = log2FoldChange, 
                                                 y = -log10(padj), 
                                                 col = custom.label, 
                                                 alpha = 1)) +
      geom_vline(xintercept = c(-fc.limit, fc.limit), col = "gray", linetype = 'dashed') +
      geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
      xlim(-7.5, 7.5) + 
      labs(x = x.axis.title,
           y = "-log10 adjusted p-value", 
           title = title) + 
      geom_point(size = 2) +
      scale_color_manual(legend.title, 
                         values = contrast.level.colors, 
                         breaks = c("DOWN", "UP")) + 
      geom_text_repel(data = de.results.labeled,
                      aes(x = log2FoldChange, 
                          y = -log10(padj), 
                          label = deglabel, 
                          col = custom.label), 
                      max.overlaps = Inf, 
                      size = 6, 
                      show.legend = FALSE) + 
      xlim(-log2.scale-1, log2.scale+1) + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      scale_alpha_identity(guide = "none")
    
  }
  
  
  
  # Make the volcano plot
  
  
  return(list("volcano.plot" = volcano.plot))
  
} 

make_heatmap <- function(deseq.output, 
                         contrast, 
                         contrast.levels, 
                         annotation, 
                         remove.genes, 
                         cluster.rows, 
                         show.row.names, 
                         cluster.cols, 
                         show.col.names){
  
  # Create a heatmap for the DEGs from the DEG contrast
  
  # Get a list of the DEGS
  rownames(deseq.output$deg.list) <- deseq.output$deg.list$gene
  contrast.deg.list <-   rownames(deseq.output$deg.list[deseq.output$deg.list$padj < 0.05,     ])
  
  subtitle <- "padj < 0.05"
  
  # Loosen the cutoff if only 1 DEG found (causes an error in pheatmap)
  if(length(contrast.deg.list) < 2){ 
    
    print("No DEGs found for padj < 0.05, trying pvalue < 0.05")
    
    contrast.deg.list <-   rownames(deseq.output$deg.list[deseq.output$deg.list$pvalue <   0.05,
    ])
    subtitle <- "p-value < 0.05"
    
    # Still no DEGs after loosening threshold
    if(length(contrast.deg.list) < 2){
      
      print(paste0("No suitable DEGs to plot on heatmap for ", comparison))
      stop()
      
    }
    
  }
  
  
  # Remove the positive and negative probes
  # Identify the NEG and POS probes
  NEG.indices <- grep("NEG_", contrast.deg.list)
  POS.indices <- grep("POS_", contrast.deg.list)
  control.indices <- c(NEG.indices, POS.indices)
  
  # Remove the control probes
  if(length(control.indices) > 0){
    goi <- contrast.deg.list[-control.indices]
  } else{
    goi <- contrast.deg.list
  }
  
  
  # Gather the annotation for the contrast
  anno.columns <- c(contrast)
  annotation <- as.data.frame(annotation)
  
  # Ensure the row names are the sample names
  rownames(annotation) <- annotation$SampleID
  contrast.annotation <- annotation[ ,anno.columns, drop = FALSE]
  
  # Subset heatmap annotation for any samples taken out in counts
  # Gather the samples in both
  anno.samples.rownames <- rownames(contrast.annotation)
  counts.samples.colnames <- colnames(deseq.output$dds)
  
  # Find the samples that are present in both
  matching.names <- intersect(anno.samples.rownames, counts.samples.colnames)
  
  # Subset counts to include only overlapping samples
  heatmap.annotation <- contrast.annotation[matching.names, , drop = FALSE]
  
  # Order by contrast level
  order.rows <- factor(heatmap.annotation[[contrast]], levels = contrast.levels)
  heatmap.annotation <- heatmap.annotation[order(order.rows), , drop = FALSE]
  
  # Reorder the data columns to they match the annotation row order
  row.name.order <- rownames(heatmap.annotation)
  
  log.norm.counts <- log2(counts(deseq.output$dds, normalized = TRUE) + 1)
  log.norm.counts <- log.norm.counts[, row.name.order]
  
  # Remove custom gene input
  if(!is.null(remove.genes)){
    
    log.norm.counts <- log.norm.counts[!rownames(log.norm.counts) %in% remove.genes, ]
    
    goi <- setdiff(goi, remove.genes)
    
  }
  
  # Define the annotation colors for the contrast
  anno.colors <- setNames(list(c("cyan", "magenta")), contrast)
  
  # Assign the levels to the colors
  anno.colors[[contrast]] <- setNames(anno.colors[[contrast]], contrast.levels)
  
  
  heatmap.contrast <- pheatmap(log.norm.counts[goi, ], 
                               main = paste0(contrast, " DEGs (", subtitle, ")"), 
                               show_rownames = show.row.names, 
                               scale = "row", 
                               show_colnames = show.col.names,
                               border_color = NA, 
                               cluster_rows = cluster.rows, 
                               cluster_cols = cluster.cols, 
                               clustering_method = "average", 
                               clustering_distance_rows = "correlation", 
                               clustering_distance_cols = "correlation", 
                               color = colorRampPalette(c("blue", "white", "red"))(120),   
                               annotation_row = NA, 
                               annotation_col = heatmap.annotation, 
                               annotation_colors = anno.colors
  )
  
  
  return(list("heatmap.contrast" = heatmap.contrast))
  

  
  
}
