# Compute summary stats of empirical null CDFs for permutations. This is a helper function for computing FDRs.

scan1perm_stats = function(genoprobs, pheno, addcovar = NULL, kinship = NULL, intcovar = NULL,
                         nperm = 100, breaks = seq(0, 1, 0.0001), scan1_out = NULL){

  # Compute mean and variance of empirical CDFs
  mn_ecdf = matrix(0, nrow = length(breaks) - 1, ncol = 1)
  sq_ecdf = matrix(0, nrow = length(breaks) - 1, ncol = 1)
  min_p = matrix(0, nrow = nperm, ncol = 1)
  for(n in 1:nperm){
    null_p = sample_perm_null(genoprobs, pheno, addcovar = addcovar, kinship = K, intcovar = intcovar)
    
    # Compute eCDF data 
    curr_hist = hist(null_p, breaks = breaks, plot = FALSE)
    curr_ecdf = as.matrix(cumsum(curr_hist$counts) / sum(curr_hist$counts))
    mn_ecdf = mn_ecdf + curr_ecdf / nperm
    sq_ecdf = sq_ecdf + curr_ecdf^2 / (nperm - 1)
    
    # Get minimum p-value
    min_p[n] = min(null_p)
  }
  
  # Compute over-dispersion parameter
  m = sum(curr_hist$counts)
  numer = sq_ecdf - (nperm / (nperm - 1)) * mn_ecdf ^ 2
  denom = mn_ecdf * (1 - mn_ecdf)
  phi = m * numer / denom 
  
  # Compute true scan
  if(is.null(scan1_out)){
    scan1_out = scan1(genoprobs, pheno, kinship = kinship, addcovar = covar, intcovar = intcovar)
  }
  nind = length(intersect(rownames(addcovar), rownames(pheno)))
  k = dim(genoprobs[[1]])[2] # Number of alleles
  true_hist = hist(lod2p(scan1_out, nind = nind, k = 8), breaks = breaks, plot = FALSE)
  true_ecdf = as.matrix(cumsum(true_hist$counts) / sum(true_hist$counts))
  
  # Compute variance of log(FDR)
  var_ind = (1 / mn_ecdf + 1 / (1 - mn_ecdf)) / (m * nperm) + (1 / true_ecdf + 1 / (1 - true_ecdf)) / m
  var_fdr = phi * var_ind
  
  fdr_hat = (mn_ecdf / true_ecdf) * ((1 - true_ecdf) / (1 - mn_ecdf))
  fdr_hat = apply(fdr_hat, 1, function(x){min(x, 1)})
  
  # Output
  null_stats = list(fdr_hat, var_fdr, min_p, breaks, k, nind)
  names(null_stats) = c("fdr_hat", "var_fdr", "min_p", "breaks", "k", "nind")
  return(null_stats)
}