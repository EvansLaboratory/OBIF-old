#' Full factorial analysis
#'
#'  Primary function to perform full factorial analysis in an Ananlysis-Ready dataset
#'
#'  @param dataAR A numeric data matrix with the analyis-ready feature expression values where row are features and columns are samples
#'
#'  @param design A design matrix to be applied for expression analysis
#'
#'  @param contrast A contrast matrix to be applied for contrast analysis
#'
#'  @return A data frame with the individual p-values per feature for main effects, simple main effects and interaction effects.
#'
#' @export
obif.ffa <- function(dataAR, design, contrast) {
  fit <- lmFit(dataAR, design)
  fit.ebayes <- eBayes(fit)
  adj.Pvalues.lmBH <- apply(fit.ebayes$p.value, 2, p.adjust, method="BH")
  adj.Pvalues.lmBonf <- apply(fit.ebayes$p.value, 2, p.adjust, method="bonferroni")

  fit2 <- contrasts.fit(fit,contrast)
  fit2.ebayes <- eBayes(fit2)
  adj.Pvalues.contrBH <- apply(fit2.ebayes$p.value, 2, p.adjust, method="BH")
  adj.Pvalues.contrBonf <- apply(fit2.ebayes$p.value, 2, p.adjust, method="bonferroni")
  res.ffa <- data.frame(pval.exp=fit.ebayes$p.value[,2:4],
                        pval.con=fit2.ebayes$p.value,
                        adjpvBH.exp=adj.Pvalues.lmBH[,2:4],
                        adjpvBH.con=adj.Pvalues.contrBH,
                        adjpvBonf.exp=adj.Pvalues.lmBonf[,2:4],
                        adjpvBonf.con=adj.Pvalues.contrBonf)
  return(res.ffa)
}
