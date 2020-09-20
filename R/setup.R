#' Set up data for MR-MtRobin algorithm.
#'
#' Sets up data to run the main MR-MtRobin algorithm:
#' a Multi-Tissue transcriptome-wide
#' Mendelian Randomization method ROBust to INvalid instrumental variables.
#'
#' @param geneID character string of the gene to be tested.
#' @param snpID vector of variant identifiers be used as instruments for \code{geneID}.
#' @param eqtl_data data.frame of summary statistics from eQTL study.
#' @param gwas_data data.frame of summary statistics from GWAS study.
#' @param LD matrix of LD correlation coefficients (\eqn{r}, not \eqn{r^2}).
#' @param nTiss integer of the number of tissues analyzed by the eQTL study.
#'
#' @return A list with the following elements:
#'
#' \tabular{ll}{
#' \code{snpID} \tab vector of variant identifiers to be used as instrumental variables.\cr
#' \code{gwas_betas} \tab vector of coefficient estimates (betas) from GWAS study.\cr
#' \code{gwas_se} \tab vector of standard errors for coefficient estimates from GWAS study.\cr
#' \code{eqtl_betas} \tab matrix of coefficient estimates (betas) from eQTL study.\cr
#' \code{eqtl_se} \tab matrix of standard errors for coefficient estimates from eQTL study.\cr
#' \code{eqtl_pvals} \tab matrix of p-values from eQTL study.\cr
#' \code{LD} \tab matrix of LD correlation coefficients (\eqn{r}, not \eqn{r^2}).\cr
#' }
#'
#' The returned list elements can be passed to the main \code{\link{MR_MtRobin}} function.
#'
#' @details The following are additional details describing the input arguments.
#'
#' The data.frame \code{eqtl_data} should have the character
#' variables \code{gene_id} and \code{variant_id} (as identifiers), as well as the variables
#' \code{beta_j}, \code{SE_j} and \code{pvalue_j} for \eqn{j} in \{1,...,\code{nTiss}\},
#' corresponding to the coefficient, standard error, and p-value estimates from each respective eQTL analysis.
#'
#' The data.frame \code{gwas_data} should have the character
#' variable \code{variant_id}, as well as the variables
#' \code{beta} and \code{SE} corresponding to the coefficient and
#' standard error estimates from the GWAS analysis.
#'
#' Note that the matrix \code{LD} should hold \emph{correlation coefficients}
#' (i.e. \eqn{r}), not their squared values (\eqn{r^2}).
#'
#' @export
#'
MR_MtRobin_setup <- function(geneID, snpID, eqtl_data, gwas_data, LD, nTiss){

  ## check gene is present and subset
  if(!(geneID %in% eqtl_data$gene_id)) stop("Gene is missing in eQTL dataset")
  eqtl_data <- subset(eqtl_data, gene_id == geneID)

  ## confirm all SNPs are present in the dataset
  if(!all(snpID %in% eqtl_data$variant_id)) stop("Some SNPs missing in eQTL dataset")
  if(!all(snpID %in% gwas_data$variant_id)) stop("Some SNPs missing in GWAS dataset")
  if(!all(snpID %in% colnames(LD))) stop("Some SNPs missing in LD matrix")

  ## align data
  eqtl_data <- eqtl_data[match(snpID,eqtl_data$variant_id),]
  gwas_data <- gwas_data[match(snpID,gwas_data$variant_id),]
  LD <- LD[match(snpID,rownames(LD)),match(snpID,colnames(LD))]

  eqtl_betas <- as.matrix(eqtl_data[,paste0("beta_",1:nTiss)])
  eqtl_se <- as.matrix(eqtl_data[,paste0("SE_",1:nTiss)])
  eqtl_pvals <- as.matrix(eqtl_data[,paste0("pvalue_",1:nTiss)])
  gwas_betas <- gwas_data$beta
  gwas_se <- gwas_data$SE

  row.names(eqtl_betas) <- row.names(eqtl_se) <- row.names(eqtl_pvals) <- eqtl_data$variant_id

  return(list(snpID=snpID,gwas_betas=gwas_betas,gwas_se=gwas_se,
              eqtl_betas=eqtl_betas,eqtl_se=eqtl_se,eqtl_pvals=eqtl_pvals,LD=LD))
}


#' Select instrumental variables for the MR-MtRobin analysis using forward selection procedure.
#'
#' Selects a set of instrumental variables (IVs) based on specified criteria to be used in
#' two-sample Mendelian Randomization analysis using MR-MtRobin. The function iteratively
#' selects the SNP/variant with the smallest minimum \eqn{P}-value of association with
#' expression of \code{geneID} and having pairwise LD \eqn{r^2} less than \code{ld_thresh}
#' with each SNP already selected. Stops selections when no candidate SNPs
#' have \eqn{P}-value less than \code{pval_thresh} in at least \code{nTiss_thresh} tissues.
#'
#' @param geneID character string of the gene to be tested.
#' @param eqtl_data data.frame of summary statistics from eQTL study.
#' @param nTiss integer of the number of tissues analyzed by the eQTL study.
#' @param LD matrix of LD correlation coefficients (\eqn{r}, not \eqn{r^2}).
#' @param ld_thresh numeric of pairwise LD threshold (\eqn{r^2}) .
#' @param pval_thresh numeric of \eqn{P}-value threshold for SNP-tissue pair to be used as IV.
#' @param nTiss_thresh integer of minimum number of tissues in which a candidate IV must have
#' \eqn{P}-value less than \code{pval_thresh}.
#'
#' @return A character vector of the SNP/variant identifiers to be used as instrumental variables with \code{geneID}.
#'
#' @details The following are additional details describing the input arguments.
#'
#' The data.frame \code{eqtl_data} should have the character
#' variables \code{gene_id} and \code{variant_id} (as identifiers), as well as the variables
#' \code{pvalue_j} and \code{beta_j} for \eqn{j} in \{1,...,\code{nTiss}\}, corresponding to the
#' \eqn{P}-value testing the null hypothesis of no association between \code{gene_id} and \code{variant_id}
#' and the estimated coefficient (beta) of the marginal effect of \code{variant_id} on \code{gene_id} expression
#' in each respective eQTL analysis, respectively. Note that the names of the p-value and beta columns
#' must match the conventions \code{pvalue_j} and \code{beta_j} for \eqn{j} in \{1,...,\code{nTiss}\} exactly.
#' the convention \code{pvalue_j} for \eqn{j} in \{1,...,\code{nTiss}\} exactly.
#'
#' Note that the matrix \code{LD} should hold \emph{correlation coefficients}
#' (i.e. \eqn{r}), not their squared values (\eqn{r^2}).
#'
#' @export
#'
select_IV_fs <- function(geneID, eqtl_data, nTiss, LD, ld_thresh=0.5, pval_thresh=0.001, nTiss_thresh=3){

  ## subset eQTL dataset to gene of interest
  if(!(geneID %in% eqtl_data$gene_id)) stop("Gene is missing in eQTL dataset")
  eqtl_data <- subset(eqtl_data, gene_id == geneID)

  # pval_mat <- as.matrix(eqtl_data[,paste0("pvalue_",1:nTiss)])
  pval_mat <- as.matrix(subset(eqtl_data,select=paste0("pvalue_",1:nTiss)))
  eqtl_data$min_p <- matrixStats::rowMins(pval_mat)
  eqtl_data$nTissThresh_p <- apply(pval_mat, 1, function(x) sort(x)[nTiss_thresh])
  eqtl_data <- subset(eqtl_data, nTissThresh_p < pval_thresh)
  ## ensure same sign for tissues with p < pval_thresh
  betas <- as.matrix(subset(eqtl_data,select=paste0("beta_",1:nTiss)))
  pvals <- as.matrix(subset(eqtl_data,select=paste0("pvalue_",1:nTiss)))
  pvals_lt001 <- pvals<0.001
  betas_sig <- betas*pvals_lt001
  sig_neg_eQTL <- as.numeric(apply(betas_sig,1,function(x) any(x < -1e-15)))
  sig_pos_eQTL <- as.numeric(apply(betas_sig,1,function(x) any(x > 1e-15)))
  eqtl_data$pos_sig_eQTL <- sig_pos_eQTL
  eqtl_data$neg_sig_eQTL <- sig_neg_eQTL
  eqtl_data <- subset(eqtl_data, !(pos_sig_eQTL==1 & neg_sig_eQTL==1))

  if(nrow(eqtl_data)==0) stop(paste0("No cis-SNPs for gene with p<",pval_thresh," in at least ",nTiss_thresh," tissues."))

  ## obtain set of candidate SNPs
  SNP_pool <- eqtl_data$variant_id
  if(length(SNP_pool)==1){
    return(SNP_pool)
  }

  ## consider allowing option to pass genotype matrix and calculate LD on smaller set
  # geno_dt_sub <- subset(geno_dt,variant_id %in% SNP_pool) ## subset large genotype data.table to current SNP pool
  # geno_mat <- t(geno_dt_sub[,2:ncol(geno_dt_sub)])
  # colnames(geno_mat) <- geno_dt_sub$variant_id
  # LD_r2 <- cor(geno_mat,use="pairwise.complete")^2

  ## restrict LD matrix to SNPs in data set
  LD <- LD[eqtl_data$variant_id,eqtl_data$variant_id]
  LD_r2 <- LD^2

  selected_snps <- NULL

  ## iteratively select best SNP remaining in SNP_pool based on minimum p-value
  while(length(SNP_pool) > 1){
    ## identify top SNP in pool
    eqtl_data <- subset(eqtl_data, variant_id %in% SNP_pool)
    next_snp_idx <- which.min(eqtl_data$min_p)
    next_snp <- eqtl_data$variant_id[next_snp_idx]
    selected_snps <-c(selected_snps,next_snp)

    ## identify SNPs not in high LD with top SNP in pool
    LD_idx <- which(row.names(LD_r2)==next_snp)
    LD_keep_idx <- which(LD_r2[LD_idx,] < ld_thresh)
    LD_r2 <- LD_r2[LD_keep_idx,LD_keep_idx]

    SNP_pool <- names(LD_keep_idx)
  }

  selected_snps <- c(selected_snps, SNP_pool)

  return(selected_snps)

}

#' Select instrumental variables for the MR-MtRobin analysis using backward selection procedure.
#'
#' Selects a set of instrumental variables (IVs) based on specified criteria to be used in
#' two-sample Mendelian Randomization analysis using MR-MtRobin. After restricting the
#' candidate SNPs to those with \eqn{P}-value less than \code{pval_thresh} in at least
#' \code{nTiss_thresh} tissues, the function iteratively
#' removes one SNP/variant from the pair with the largest pairwise LD \eqn{r^2} until all
#' pairwise LD \eqn{r^2} among remaining SNPs is less than \code{ld_thresh}. The variant
#' selected to be removed is based on the ranked criteria 1) greatest number
#' of pairwise LD \eqn{r^2} > \code{ld_thresh} with other remaining SNPs;
#' 2) fewest number of tissues with \eqn{P}-value less than \code{pval_thresh}; and
#' 3) largest minimum \eqn{P}-value of association with expression of \code{geneID}.
#'
#' @param geneID character string of the gene to be tested.
#' @param eqtl_data data.frame of summary statistics from eQTL study.
#' @param nTiss integer of the number of tissues analyzed by the eQTL study.
#' @param LD matrix of LD correlation coefficients (\eqn{r}, not \eqn{r^2}).
#' @param ld_thresh numeric of pairwise LD threshold (\eqn{r^2}) .
#' @param pval_thresh numeric of \eqn{P}-value threshold for SNP-tissue pair to be used as IV.
#' @param nTiss_thresh integer of minimum number of tissues in which a candidate IV must have
#' \eqn{P}-value less than \code{pval_thresh}.
#'
#' @return A character vector of the SNP/variant identifiers to be used as instrumental variables with \code{geneID}.
#'
#' @details The following are additional details describing the input arguments.
#'
#' The data.frame \code{eqtl_data} should have the character
#' variables \code{gene_id} and \code{variant_id} (as identifiers), as well as the variables
#' \code{pvalue_j} and \code{beta_j} for \eqn{j} in \{1,...,\code{nTiss}\}, corresponding to the
#' \eqn{P}-value testing the null hypothesis of no association between \code{gene_id} and \code{variant_id}
#' and the estimated coefficient (beta) of the marginal effect of \code{variant_id} on \code{gene_id} expression
#' in each respective eQTL analysis, respectively. Note that the names of the p-value and beta columns
#' must match the conventions \code{pvalue_j} and \code{beta_j} for \eqn{j} in \{1,...,\code{nTiss}\} exactly.
#'
#' Note that the matrix \code{LD} should hold \emph{correlation coefficients}
#' (i.e. \eqn{r}), not their squared values (\eqn{r^2}).
#'
#' @export
#'
select_IV_bs <- function(geneID, eqtl_data, nTiss, LD, ld_thresh=0.5, pval_thresh=0.001, nTiss_thresh=3){

  ## subset eQTL dataset to gene of interest
  if(!(geneID %in% eqtl_data$gene_id)) stop("Gene is missing in eQTL dataset")
  eqtl_data <- subset(eqtl_data, gene_id == geneID)

  pval_mat <- as.matrix(subset(eqtl_data,select=paste0("pvalue_",1:nTiss)))
  eqtl_data$min_p <- matrixStats::rowMins(pval_mat)
  eqtl_data$nTissThresh_p <- apply(pval_mat, 1, function(x) sort(x)[nTiss_thresh])
  eqtl_data$nTiss_ltThresh <- apply(pval_mat, 1, function(x) sum(x<pval_thresh))
  eqtl_data <- subset(eqtl_data, nTissThresh_p < pval_thresh)
  ## ensure same sign for tissues with p < pval_thresh
  betas <- as.matrix(subset(eqtl_data,select=paste0("beta_",1:nTiss)))
  pvals <- as.matrix(subset(eqtl_data,select=paste0("pvalue_",1:nTiss)))
  pvals_lt001 <- pvals<0.001
  betas_sig <- betas*pvals_lt001
  sig_neg_eQTL <- as.numeric(apply(betas_sig,1,function(x) any(x < -1e-15)))
  sig_pos_eQTL <- as.numeric(apply(betas_sig,1,function(x) any(x > 1e-15)))
  eqtl_data$pos_sig_eQTL <- sig_pos_eQTL
  eqtl_data$neg_sig_eQTL <- sig_neg_eQTL
  eqtl_data <- subset(eqtl_data, !(pos_sig_eQTL==1 & neg_sig_eQTL==1))

  if(nrow(eqtl_data)==0) stop(paste0("No cis-SNPs for gene with p<",pval_thresh," in at least ",nTiss_thresh," tissues."))

  ## obtain set of candidate SNPs
  SNP_pool <- eqtl_data$variant_id
  if(length(SNP_pool)==1){
    return(SNP_pool)
  }

  ## consider allowing option to pass genotype matrix and calculate LD on smaller set
  # geno_dt_sub <- subset(geno_dt,variant_id %in% SNP_pool) ## subset large genotype data.table to current SNP pool
  # geno_mat <- t(geno_dt_sub[,2:ncol(geno_dt_sub)])
  # colnames(geno_mat) <- geno_dt_sub$variant_id
  # LD_r2 <- cor(geno_mat,use="pairwise.complete")^2

  ## restrict LD matrix to SNPs in data set
  LD <- LD[eqtl_data$variant_id,eqtl_data$variant_id]
  LD_r2 <- LD^2
  ## get pairwise r^2
  LD_r2_pw <- LD_r2
  LD_r2_pw[lower.tri(LD_r2_pw,diag=T)] <- NA
  LD_r2_pw <- reshape2::melt(LD_r2_pw)
  colnames(LD_r2_pw) <- c("SNP1","SNP2","r2")
  LD_r2_pw$SNP1 <- as.character(LD_r2_pw$SNP1)
  LD_r2_pw$SNP2 <- as.character(LD_r2_pw$SNP2)
  ## drop SNPs paired with themselves and duplicates (resulting from lower triangle, set to NA)
  LD_r2_pw <- subset(LD_r2_pw, SNP1 != SNP2 & !is.na(r2))
  ##sort by descending r^2
  LD_r2_pw <- LD_r2_pw[order(LD_r2_pw$r2,decreasing = T),]

  ## check for case all candidate SNPs have pw r^2 > thresh (return SNP with smallest p-value)
  if(min(LD_r2_pw$r2) >=ld_thresh){
    return(eqtl_data$variant_id[which.max(eqtl_data$min_p)])
  }


  while(max(LD_r2_pw$r2)>=ld_thresh){
    ## identify top SNP in pool
    eqtl_data <- subset(eqtl_data, variant_id %in% SNP_pool)
    LD_r2_pw <- subset(LD_r2_pw, SNP1 %in% SNP_pool & SNP2 %in% SNP_pool)
    LD_r2_pw_geThresh <- subset(LD_r2_pw, r2>=ld_thresh)

    ##identify SNP pair with highest remaining r^2
    next_pair <- LD_r2_pw[1,]

    ##criteria 1: drop the SNP with largest number of correlations > threshold with other SNPs
    SNP1_nSNP_geThresh <- sum(LD_r2_pw_geThresh$SNP1 == next_pair$SNP1) + sum(LD_r2_pw_geThresh$SNP2 == next_pair$SNP1)
    SNP2_nSNP_geThresh <- sum(LD_r2_pw_geThresh$SNP1 == next_pair$SNP2) + sum(LD_r2_pw_geThresh$SNP2 == next_pair$SNP2)
    if(SNP1_nSNP_geThresh > SNP2_nSNP_geThresh){
      next_SNP_drop <- next_pair$SNP1
    } else if (SNP1_nSNP_geThresh < SNP2_nSNP_geThresh){
      next_SNP_drop <- next_pair$SNP2
    } else{
      ##criteria 2: drop the SNP with smaller number of tissues with p<0.001
      next_pair_eqtl <- subset(eqtl_data, variant_id %in% c(next_pair$SNP1,next_pair$SNP2))
      SNP1_nTiss <- subset(next_pair_eqtl,variant_id==next_pair$SNP1)$nTiss_ltThresh
      SNP2_nTiss <- subset(next_pair_eqtl,variant_id==next_pair$SNP2)$nTiss_ltThresh

      if(SNP1_nTiss > SNP2_nTiss){
        next_SNP_drop <- next_pair$SNP2
      } else if(SNP1_nTiss < SNP2_nTiss){
        next_SNP_drop <- next_pair$SNP1
      } else{
        ##criteria 3: drop the SNP with larger minimum p-value
        next_SNP_drop <- next_pair_eqtl$variant_id[which.max(next_pair_eqtl$min_p)]
      }
    }

    SNP_pool <- SNP_pool[which(SNP_pool != next_SNP_drop)]

  }

  return(SNP_pool)

}


#' Select instrumental variables for the MR-MtRobin analysis.
#'
#' Selects a set of instrumental variables (IVs) based on specified criteria to be used in
#' two-sample Mendelian Randomization analysis using MR-MtRobin.
#'
#' @param geneID character string of the gene to be tested.
#' @param eqtl_data data.frame of summary statistics from eQTL study.
#' @param nTiss integer of the number of tissues analyzed by the eQTL study.
#' @param LD matrix of LD correlation coefficients (\eqn{r}, not \eqn{r^2}).
#' @param ld_thresh numeric of pairwise LD threshold (\eqn{r^2}) .
#' @param pval_thresh numeric of \eqn{P}-value threshold for SNP-tissue pair to be used as IV.
#' @param nTiss_thresh integer of minimum number of tissues in which a candidate IV must have
#' \eqn{P}-value less than \code{pval_thresh}.
#' @param back_select logical of whether to use backward selection (\code{TRUE}) or forward
#' selection (\code{FALSE}) procedure in selecting set of instrumental variables.
#'
#' @return A character vector of the SNP/variant identifiers to be used as instrumental variables with \code{geneID}.
#'
#' @details The following are additional details describing the input arguments.
#'
#' The data.frame \code{eqtl_data} should have the character
#' variables \code{gene_id} and \code{variant_id} (as identifiers), as well as the variables
#' \code{pvalue_j} and \code{beta_j} for \eqn{j} in \{1,...,\code{nTiss}\}, corresponding to the
#' \eqn{P}-value testing the null hypothesis of no association between \code{gene_id} and \code{variant_id}
#' and the estimated coefficient (beta) of the marginal effect of \code{variant_id} on \code{gene_id} expression
#' in each respective eQTL analysis, respectively. Note that the names of the p-value and beta columns
#' must match the conventions \code{pvalue_j} and \code{beta_j} for \eqn{j} in \{1,...,\code{nTiss}\} exactly.
#'
#' Note that the matrix \code{LD} should hold \emph{correlation coefficients}
#' (i.e. \eqn{r}), not their squared values (\eqn{r^2}).
#'
#' @export
#'
select_IV <- function(geneID, eqtl_data, nTiss, LD, ld_thresh=0.5, pval_thresh=0.001, nTiss_thresh=3,back_select=TRUE){

  if(back_select){
    return(select_IV_bs(geneID, eqtl_data, nTiss, LD, ld_thresh, pval_thresh, nTiss_thresh))
  } else{
    return(select_IV_fs(geneID, eqtl_data, nTiss, LD, ld_thresh, pval_thresh, nTiss_thresh))
  }

}
