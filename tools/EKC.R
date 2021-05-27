EKC <- function (R, N = NA, use = c("pairwise.complete.obs", "all.obs", 
                             "complete.obs", "everything", "na.or.complete"), cor_method = c("pearson", 
                                                                                             "spearman", "kendall")) 
{

  p <- ncol(R)
  lambda <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  refs <- vector("double", p)
  for (i in seq_len(p)) {
    refs[i] <- max(((1 + sqrt(p/N))^2) * (p - sum(refs))/(p - 
                                                            i + 1), 1)
  }
  out <- list(eigenvalues = lambda, n_factors = which(lambda <= 
                                                        refs)[1] - 1, references = refs, settings = list(use = use, 
                                                                                                         cor_method = cor_method, N = N))
  class(out) <- "EKC"
  return(out)
}