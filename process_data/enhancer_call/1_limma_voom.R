#!/usr/bin/env Rscript

######################################################################
# Modified TMM normalization
# This version uses user-defined negative controls instead of M/A cutoffs
# to identify non-regulatory fragments for normalization.
# Based on TMM from edgeR (Mark Robinson and Gordon Smyth)

calcFactorQuantile = function(data, lib.size, p=0.75)
#	Generalized version of upper-quartile normalization
#	Mark Robinson and Gordon Smyth
#	Created 16 Aug 2010. Last modified 12 Sep 2020.
{
	f <- rep_len(1,ncol(data))
	for (j in seq_len(ncol(data))) f[j] <- quantile(data[,j], probs=p)
	if(min(f)==0) warning("One or more quantiles are zero")
	f / lib.size
}

######################################################################
calcFactorTMM = function(obs, ref, neg_idx, libsize.obs=NULL, libsize.ref=NULL, logratioTrim=.30, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, expand)
# Compute TMM normalization factor between two libraries
# Incorporates negative control fragments instead of M/A filtering
    
{
    obs <- as.numeric(obs)
    ref <- as.numeric(ref)
    
    # Compute library sizes if not provided
    if( is.null(libsize.obs) ) nO <- sum(obs) else nO <- libsize.obs
    if( is.null(libsize.ref) ) nR <- sum(ref) else nR <- libsize.ref

    # Calculate log-ratios and average expression
    logR <- log2((obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size
    absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression
    v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref   # estimated asymptotic variance
    
    if(length(neg_idx) > 0){
        # Use user-supplied negative control indices
        logR_neg = logR[neg_idx]
        absE_neg = absE[neg_idx]
        v_neg = v[neg_idx]

        # Remaining fragments
        logR_int = logR[-neg_idx]
        absE_int = absE[-neg_idx]
        v_int = v[-neg_idx]

        # Remove infinite values and apply absolute expression cutoff
        fin_int <- is.finite(logR_int) & is.finite(absE_int) & (absE_int > Acutoff)
        fin_neg = is.finite(logR_neg) & is.finite(absE_neg) & (absE_neg > Acutoff)

        logR_neg = logR_neg[fin_neg]
        absE_neg = absE_neg[fin_neg]
        v_neg = v_neg[fin_neg]

        logR_int = logR_int[fin_int]
        absE_int = absE_int[fin_int]
        v_int = v_int[fin_int]

        if(max(abs(logR_neg)) < 1e-6) return(1)
        
        # Compute log-ratio and absolute expression trimming thresholds
        Q = quantile(logR_neg, probs = c(logratioTrim, 1-logratioTrim), na.rm=FALSE)
        loL = Q[1]
        hiL = Q[2]

        Q = quantile(absE_neg, probs=c(sumTrim, 1-sumTrim), na.rm=FALSE)
        loS = Q[1]
        hiS = Q[2]

        # Identify trimmed subset of negative controls
        logR_keep = which(logR_neg <= hiL & logR_neg >= loL)
        absE_keep = which(absE_neg <= hiS & absE_neg >= loS)
        keep_neg = intersect(logR_keep, absE_keep)

        if(expand==TRUE){
            # Optionally expand trimming to non-negative control fragments
            logR_keep = which(logR_int <= hiL & logR_int >= loL)
            absE_keep = which(absE_int <= hiS & absE_int >= loS)
            keep_int = intersect(logR_keep, absE_keep)

            # Combine negative control and additional fragments
            logR_neg = c(logR_neg[keep_neg], logR_int[keep_int])
            v_neg = c(v_neg[keep_neg], v_int[keep_int])
        }else{
            logR_neg = logR_neg[keep_neg]
            v_neg = v_neg[keep_neg]
        }
        
    }else{ # Standard TMM: no negative controls provided
        
        #	remove infinite values, cutoff based on A
        fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
        
        logR = logR[fin]
        absE = absE[fin]
        v = v[fin]
        
        if(max(abs(logR)) < 1e-6) return(1)
        
        #	taken from the original mean() function
        n <- length(logR)
        loL <- floor(n * logratioTrim) + 1
        hiL <- n + 1 - loL
        loS <- floor(n * sumTrim) + 1
        hiS <- n + 1 - loS

        #	keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
        #	a fix from leonardo ivan almonacid cardenas, since rank() can return
        #	non-integer values when there are a lot of ties
        keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)
        
        logR_neg = logR[keep]
        v_neg = v[keep]
    }
    
    if(doWeighting){
        f <- sum(logR_neg/v_neg, na.rm=TRUE) / sum(1/v_neg, na.rm=TRUE)
    }else{
        f <- mean(logR_neg, na.rm=TRUE)
    }

    #	Results will be missing if the two libraries share no features with positive counts
    #	In this case, return unity
    if(is.na(f)) f <- 0
    2^f
}
######################################################################
######################################################################
calcNormFactor_modified = function(data, idx, refColumn=NULL, expand, logratioTrim=.30, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10)
# Compute normalization factors for all samples in a count matrix
# Uses modified TMM with optional negative control fragments
{ 
    x = as.matrix(data)
    nsamples = ncol(x)
    lib.size = colSums(x)
    
    # Remove rows with all zero counts
    allzero <- .rowSums(x>0, nrow(x), nsamples) == 0L
    if(any(allzero)) x <- x[!allzero,,drop=FALSE]
    
    # Index of negative control rows (if any)
    if(idx == 0){
        neg_idx = c()
    }else{
        neg_idx = c(1:idx)
    }
    
    # Choose reference sample based on upper-quartile
    if( is.null(refColumn)){
        f75 = calcFactorQuantile(data=x, lib.size = lib.size, p = 0.75)
        if(median(f75) < 1e-20){
            refColumn = which.max(colSums(sqrt(x)))
        }else{
            refColumn = which.min(abs(f75-mean(f75)))
        } 
    }
    f <- rep_len(NA_real_,nsamples)

    # Compute TMM factor for each library
    for(i in 1:nsamples){
        f[i] = calcFactorTMM(obs=x[, i], ref=x[, refColumn], neg_idx = neg_idx, 
                             libsize.obs=lib.size[i], libsize.ref=lib.size[refColumn], 
                             expand=expand, logratioTrim, sumTrim, doWeighting=TRUE, Acutoff=-1e10)
    }
    
    # Normalize so that geometric mean is 1
    f <- f/exp(mean(log(f)))
    
    names(f) <- colnames(x)
    f # return normalization vector
}
########################################################################################


suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(edgeR)
  library(limma)
  library(DESeq2)  # kept for parity with original scripts (may not be used)
})

# Setup argument parser
parser <- ArgumentParser()

# Define input arguments
parser$add_argument("-c", "--countMatrix", help="Input UMI count")
parser$add_argument("-o", "--outDir", help="Output directory")
parser$add_argument("--num_rep_DNA", type="integer", default=3, help="Number of DNA replicates")
parser$add_argument("--num_rep_RNA", type="integer", default=3, help="Number of RNA replicates")
parser$add_argument("--corresponding", required=F, action="store_true", default=FALSE, help="If DNA/RNA are corresponding pairs")
parser$add_argument("--use_modified_tmm", action="store_true", help="Use modified TMM normalization (requires negative controls)")
parser$add_argument("--prefiltered", required=F, action="store_true", default=FALSE, help="Set if input count is already filtered")
parser$add_argument("--default_filter", required=F, action="store_true", default=FALSE, help="Use default filtering in Limma-Voom pipeline.")
parser$add_argument("--cpm_filter", required=F, action="store_true", default=FALSE, help="Filter based on nomalized read count.")
parser$add_argument("--count_cutoff", type="integer", default=10, help="Raw read count cutoff.")
parser$add_argument("--num_lib_to_filter", type="integer", default=0, help="Minimal number of libraries to pass filter")
parser$add_argument("--filter_only_DNA", action="store_true", default=FALSE, help="Filter based only on DNA libraries.")
parser$add_argument("--use_average_lib_size_in_filter", action = "store_true", default = FALSE, help = "Use average library size to compute normalized count cutoff; otherwise use the minimal library size.")
# Parse arguments
args <- parser$parse_args()

dir.create(args$outDir, showWarnings=FALSE, recursive=TRUE)

# ========================
# Load and prepare count matrix
# ========================
print(args$input)
count_mat = as.data.frame(fread(args$countMatrix, sep='\t', data.table=FALSE))
print(nrow(count_mat))

# Create DNA and RNA column labels
DNA_columns = paste0('DNA', c(1:args$num_rep_DNA))
RNA_columns = paste0('RNA', c(1:args$num_rep_RNA))

# How many metadata columns precede the DNA/RNA replicates?
meta_n <- ncol(count_mat) - (args$num_rep_DNA + args$num_rep_RNA)

# Define metadata names by count
meta_cols <- switch(
  as.character(meta_n),
  "4" = c("seqnames", "start", "end", "strand"),
  "5" = c("seqnames", "start", "end", "name", "strand"),
  "6" = c("seqnames", "start", "end", "ID", "strand", "annot"),
  # New: 7 metadata columns (penultimate = original annotation, last = annot)
  "7" = c("seqnames", "start", "end", "ID", "strand", "orig_annot", "annot"),
  stop(sprintf("Unsupported metadata column count: %d", meta_n))
)

colnames(count_mat) <- c(meta_cols, DNA_columns, RNA_columns)

# Create edgeR DGEList object
dge <- DGEList(
  counts = count_mat[, c(DNA_columns, RNA_columns)],
  genes  = count_mat[, meta_cols]
)

print(dge$samples$lib.size)

# Add metadata
dge$samples$group = as.factor(c(rep('DNA', args$num_rep_DNA), rep('RNA', args$num_rep_RNA)))
dge$samples$replicate = as.factor(c(1:args$num_rep_DNA, 1:args$num_rep_RNA))

# ========================
# Filtering
# ========================
if(!args$prefiltered){
    print('Filter')
    if(args$cpm_filter){
        print('Based on nomralized CPM count')
        cpm = cpm(dge, log=FALSE)

        # Choose how to define library size
        if(args$filter_only_DNA){
            print('Based only on DNA libraries')
            if(args$use_average_lib_size_in_filter){
                print('Use average lib size to calculate CPM threshold')
                lib_size = mean(dge$samples$lib.size[1:args$num_rep_DNA])
            }else{
                print('Use minimal lib size to calculate CPM threshold')
                lib_size = min(dge$samples$lib.size[1:args$num_rep_DNA])
            }
        }else{
            print('Filter based on all libs')
            if(args$use_average_lib_size_in_filter){
                print('Use average lib size to calculate CPM threshold')
                lib_size = mean(dge$samples$lib.size)
            }else{
                print('Use minimal lib size to calculate CPM threshold')
                lib_size = min(dge$samples$lib.size)
            }
        }
        print('Calculate CPM threshold')
        cutoff = args$count_cutoff / (lib_size / 1e6)
        print('Output CPM threshold')
        print(cutoff)

        if(args$num_lib_to_filter > 0){
            print('Output number of libraries need to pass filter')
            print(args$num_lib_to_filter)
            if(args$filter_only_DNA){
                if(args$num_rep_DNA == 1){
                    keep = cpm[, DNA_columns] >= cutoff
                }else{
                    keep = rowSums(cpm[, DNA_columns] >= cutoff) >= args$num_lib_to_filter
                }
            }else{
                keep = rowSums(cpm >= cutoff) >= args$num_lib_to_filter
            }
        }else{
            print('If num_lib_to_filter == 0, use default number determined by minimal number of replicate between DNA and RNA libraries')
            keep = rowSums(cpm >= cutoff) >= min(c(args$num_rep_DNA, args$num_rep_RNA))
        }
    }else{
        print("Filter with default edgeR's filterByExpr")
        # Default filtering method using edgeR’s filterByExpr
        keep = filterByExpr(dge, group = dge$samples$group)
    }

    print(table(keep))
    dge = dge[keep,, keep.lib.sizes=FALSE]
    # Ensure some RNA reads are retained
    check = rowSums(dge$counts[, RNA_columns] > 0) > 0
    print(table(check))
}

gc()  # Garbage collection

# Save filtered count matrix
if(!args$prefiltered){
    write.table(cbind(dge$genes, dge$counts), file=paste0(args$outDir, 'filtered_raw_count.txt'),
                sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
    print('finished saving count matrix')
}

# ========================
# Normalization (TMM)
# ========================

if(args$use_modified_tmm){# assay included negative control elements
    filtered = count_mat[count_mat$ID %in% dge$genes$ID, ]
    print(nrow(filtered))

    # Look for all negative ctrls
    neg_ctrl = filtered[filtered$annot == 'neg_ctrl', ]
    print(nrow(neg_ctrl))

    # elements of interests = all elements - neg_ctrl
    exclude_neg = filtered[filtered$annot != 'neg_ctrl', ]
    print(nrow(exclude_neg))

    filtered_count = rbind(neg_ctrl, exclude_neg)
    print(nrow(filtered_count))

    # Create edgeR DGEList object
    dge = DGEList(
        counts = filtered_count[, c(DNA_columns, RNA_columns)], 
        genes = filtered_count[, meta_cols])
    print(dge$samples$lib.size)

    # Add metadata
    dge$samples$group = as.factor(c(rep('DNA', args$num_rep_DNA), rep('RNA', args$num_rep_RNA)))
    dge$samples$replicate = as.factor(c(1:args$num_rep_DNA, 1:args$num_rep_RNA))

    dge$samples$norm.factors = calcNormFactor_modified(
        dge, nrow(neg_ctrl), 
        logratioTrim=0.3, 
        sumTrim=0.05, 
        refColumn=NULL,
        expand=TRUE)
    print(dge$samples$norm.factors)

}else{
    dge = calcNormFactors(dge, method="TMM")
    print(dge$samples$norm.factors)
}

# ========================
# Linear Modeling (limma-voom)
# ========================
# Create design matrix for DNA vs RNA
coldata = data.frame(
    RvD = c(rep(0, args$num_rep_DNA), rep(1, args$num_rep_RNA)),
    Rep = c(1:args$num_rep_DNA, 1:args$num_rep_RNA)
)
options(contrasts=c("contr.sum", "contr.poly"))
mdes = model.matrix(~RvD, data=coldata)

# If paired design (corresponding), apply duplicate correlation
if(args$corresponding){
    png(file=paste0(args$outDir, "mean_variance_trend_voom_plot.png"), width=5000, height=5000, res=300)
    vdata = voomWithQualityWeights(dge, design=mdes, plot=F)
    vcorrf = duplicateCorrelation(vdata, mdes, block=coldata$Rep)$consensus.correlation
    vdata = voomWithQualityWeights(dge, design=mdes, block=coldata$Rep, correlation=vcorrf, plot=T)
    vcorrf = duplicateCorrelation(vdata, mdes, block=coldata$Rep)$consensus.correlation
    dev.off()

    fit = eBayes(lmFit(vdata, mdes, block=coldata$Rep, correlation=vcorrf), robust=TRUE)
}else{
    png(file=paste0(args$outDir, "mean_variance_trend_voom_plot.png"), width=5000, height=5000, res=300)
    vdata = voomWithQualityWeights(dge, design=mdes, plot=T)
    dev.off()

    fit = eBayes(lmFit(vdata, mdes), robust=TRUE)
}

# Extract differential expression results
hits = topTable(fit, coef='RvD', number=Inf, sort.by="none")

# Write output
write.table(hits, file=paste0(args$outDir, 'limma_out.txt'), sep='\t', quote=FALSE)
print('finished')

