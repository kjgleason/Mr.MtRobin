#' Conduct two-sample Mendelian Randomization using summary statistics
#'
#' Runs the main MR-MtRobin algorithm: a Multi-Tissue transcriptome-wide
#' Mendelian Randomization method ROBust to INvalid instrumental variables
#'
#' @param snpID vector of variant identifiers to be used as instrumental variables.
#' @param gwas_betas vector of coefficient estimates (betas) from GWAS study.
#' @param gwas_se vector of standard errors for coefficient estimates from GWAS study.
#' @param eqtl_betas matrix of coefficient estimates (betas) from eQTL study.
#' @param eqtl_se matrix of standard errors for coefficient estimates from eQTL study.
#' @param eqtl_pvals matrix of p-values from eQTL study.
#' @param LD matrix of LD correlation coefficients (\eqn{r}, not \eqn{r^2}).
#' @param pval_thresh p-value threshold for instrumental variables (IVs).
#'
#' @return A list with the following elements:
#' \tabular{ll}{
#' \code{lme_res} \tab an object of class \code{lmerMod},
#' returned from the reverse regression random slope
#' mixed model run by MR-MtRobin..\cr
#' \code{gwas_res} \tab data.table of the gwas data (snpID, gwas_beta, gwas_se).\cr
#' \code{LD} \tab matrix of LD correlation coefficients (\eqn{r}, not \eqn{r^2}).
#' }
#' The last two items are returned for use by \code{MR_MtRobin_resample()}.
#'
#' To conduct inference on the returned results, use function \code{\link{MR_MtRobin_resample}}.
#'
#' @details The following are additional details describing the input arguments.
#' For \code{eqtl_betas}, \code{eqtl_se}, and \code{eqtl_pvals}
#' each row \eqn{i} corresponds to a SNP/variant/IV
#' while each column \eqn{j} holds the summary statistics in tissue type \eqn{j}.
#' Both the rows of the eQTL data and the order of the GWAS vectors should
#' match the order of \code{snpID}.
#' Note that the matrix \code{LD} should hold \emph{correlation coefficients}
#' (i.e. \eqn{r}), not their squared values (\eqn{r^2}).
#'
#' @examples
#' ## MR_MtRobin_input created using MR_MtRobin_setup()
#' ## IV_gene1 created using select_IV()
#'
#' MR_MtRobin(snpID=IV_gene1,
#' gwas_betas=MR_MtRobin_input$gwas_betas, gwas_se=MR_MtRobin_input$gwas_se,
#' eqtl_betas=MR_MtRobin_input$eqtl_betas, eqtl_se=MR_MtRobin_input$eqtl_se,
#' eqtl_pvals=MR_MtRobin_input$eqtl_pvals, LD=MR_MtRobin_input$LD)
#'
#' @export
#'
MR_MtRobin <- function(snpID,gwas_betas,gwas_se,eqtl_betas,eqtl_se,eqtl_pvals,LD,pval_thresh=0.001){

  if(ncol(eqtl_betas) != ncol(eqtl_se) | ncol(eqtl_betas) != ncol(eqtl_pvals)){
    stop("eQTL betas, SE and p-values do not all have same number of columns")
  }

  ## determine number of tissues (or conditions/studies/etc.)
  nT <- ncol(eqtl_betas)

  ## rename columns for standardized processing
  colnames(eqtl_betas) <- paste0("beta_",1:nT)
  colnames(eqtl_se) <- paste0("SE_",1:nT)
  colnames(eqtl_pvals) <- paste0("pvalue_",1:nT)

  ## combine eQTL results into data.table
  eqtl_res <- data.table::data.table(snpID,eqtl_betas,eqtl_se,eqtl_pvals)

  # reshape to long format (each SNP-gene-tissue will be row)
  eQTL_res_melt <- data.table::melt(eqtl_res,id.vars=c("snpID"),
                                    measure.vars=patterns("^beta_","^SE_","^pvalue_"),value.name=c("beta","se","pvalue"),
                                    variable.name="tissue")

  ## subset to strong instrumental variables (based on p-value threshold)
  eQTL_res_melt_PltThresh <- subset(eQTL_res_melt, pvalue < pval_thresh)

  ## set up gwas data for merging
  gwas_res <- data.table::data.table(snpID=snpID,gwas_beta=gwas_betas,gwas_se=gwas_se)

  ## merge eQTL and GWAS results
  data.table::setkey(eQTL_res_melt_PltThresh,snpID)
  data.table::setkey(gwas_res,snpID)
  merged_res <- data.table::merge.data.table(eQTL_res_melt_PltThresh,gwas_res)

  ## set up coefficients for reverse regression
  beta_x <- matrix(merged_res$beta,ncol=1)
  beta_y <- matrix(merged_res$gwas_beta,ncol=1)

  ## standard errors (for weights)
  se_x <- matrix(merged_res$se,ncol=1)
  se_y <- matrix(merged_res$gwas_se,ncol=1) ## not used by function; used in resampling (consider returning with results)

  ##identifiers
  snpID <- merged_res$snpID

  lme_res <- lme4::lmer(beta_x~(beta_y-1)+(beta_y-1|snpID),weights=1/se_x^2)

  ## return results from weighted regression with random slopes (and correlated errors) along with GWAS SE and LD (needed for resampling)
  return(list(lme_res=lme_res,gwas_res=gwas_res,LD=LD))
}


#' Obtain p-value from MR-MtRobin results
#'
#' Uses a resampling procedure to estimate a \eqn{P}-value for a MR-MtRobin object.
#'
#' @param MR_MtRobin_res a list object returned by \code{MR_MtRobin}.
#' @param nsamp integer of the number of samples to use in estimating \eqn{P}-value
#' using a null distribution.
#' @param use_nonconverge logical of whether to use samples resulting in non-convergence
#' of random slope
#'
#' @return A list of two elements:
#' \tabular{ll}{
#' \code{pvalue} \tab numeric of the estimated \eqn{P}-value.\cr
#' \code{nsamp_used} \tab integer of the number of samples used in estimating the \eqn{P}-value.\cr
#' }
#'
#' \code{nsamp_used} is returned because not all \code{nsamp} may be used if
#' \code{use_nonconverge=FALSE} (samples will be dropped if the model does not
#' converge or results in a singular fit of the random slope).
#'
#' @examples
#' myRes <- MR_MtRobin(snpID,gwas_betas,gwas_se,eqtl_betas,eqtl_se,eqtl_pvals,LD)
#'
#' MR_MtRobin_resample(MR_MtRobin_res=myRes, nsamp=1e4, use_nonconverge=FALSE)
#'
#' @export
#'
MR_MtRobin_resample <- function(MR_MtRobin_res,nsamp=1e4,use_nonconverge=FALSE){

  lme_res <- MR_MtRobin_res$lme_res
  gwas_res <- MR_MtRobin_res$gwas_res
  LD <- MR_MtRobin_res$LD

  ## extract data from MR-MtRobin results
  eqtl_betas <- lme_res@frame$beta_x
  snpID <- lme_res@frame$snpID
  weights <- lme_res@frame$`(weights)`
  tstat_MR_MtRobin <- summary(lme_res)$coefficients[1,3]
  nT <- length(eqtl_betas)/length(unique(snpID))

  ## bootstrapped null distribution, accounting for LD correlations
  gwas_se <- gwas_res$gwas_se
  beta_gwas_nullMat <- mvtnorm::rmvnorm(nsamp,mean=rep(0,length(gwas_se)),sigma=diag(gwas_se) %*% LD %*% diag(gwas_se))
  ## name columns with corresponding variant identified
  colnames(beta_gwas_nullMat) <- gwas_res$snpID

  ## initialize return data types
  tstat_nulls <- NULL
  nsamp_used <- 0

  ## run MR-MtRobin on the null models
  if(use_nonconverge){
    for(i in 1:nsamp){
      beta_gwas_null <- beta_gwas_nullMat[i,snpID]
      lme_null <- lme4::lmer(eqtl_betas~(beta_gwas_null-1)+(beta_gwas_null-1|snpID),weights=weights)
      tstat_nulls <- c(tstat_nulls, summary(lme_null)$coefficients[1,3])
      nsamp_used <- nsamp_used + 1
    }
  } else{
    for(i in 1:nsamp){
      beta_gwas_null <- beta_gwas_nullMat[i,snpID]
      lme_null <- lme4::lmer(eqtl_betas~(beta_gwas_null-1)+(beta_gwas_null-1|snpID),weights=weights)
      if(is.null(summary(lme_null)$optinfo$conv$lme4$messages)){
        tstat_nulls <- c(tstat_nulls, summary(lme_null)$coefficients[1,3])
        nsamp_used <- nsamp_used + 1
      }
    }
  }

  ## warnings when nsamp_used is low
  if(nsamp_used==0) stop("All resampled datasets failed to converge in random slope model.")
  if(nsamp_used<=(nsamp/2)) warning("More than half of resampled datasets failed to converge in random slope model.")

  pval <- mean(abs(tstat_nulls) >= abs(tstat_MR_MtRobin))

  return(list(pvalue=pval,nsamp_used=nsamp_used))
}
