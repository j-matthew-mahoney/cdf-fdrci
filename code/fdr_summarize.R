fdr_summarize = function(null_stats, fdr_cl = 0.95, fdr_thresh = c(0.05, 0.1, 0.2), geno_sig = c(0.01, 0.05), plot = TRUE, quartz = TRUE){

  alpha = 1 - fdr_cl
  
  fdr_ci_lo = exp(log(null_stats$fdr_hat) + qnorm(alpha / 2) * sqrt(null_stats$var_fdr))
  fdr_ci_hi = sapply(exp(log(null_stats$fdr_hat) - qnorm(alpha / 2) * sqrt(null_stats$var_fdr)), function(x){min(x, 1)})
  
  sl = p2lod(quantile(null_stats$min_p, sort(geno_sig)), null_stats$k, null_stats$nind)
  names(sl) = paste("Genome-wide sig.", names(sl))
  
  x_vals = p2lod(null_stats$breaks[2:length(null_stats$breaks)], null_stats$k, null_stats$nind)
  
  if(plot){
    
    if(quartz == TRUE){
      quartz()
    }
    
    xlim = c(0, 1.1 * max(sl, na.rm = TRUE))
    plot(x_vals, null_stats$fdr_hat, type = "l", ylim = c(0, 1),
         xlab = "LOD", ylab = "FDR", xlim = xlim)
    points(x_vals, fdr_ci_lo, type = "l", lty = 3)
    points(x_vals, fdr_ci_hi, type = "l", lty = 3)
    for(s in 1:length(sl)){
      abline(v = sl[s], col = "red", lty = s)
    }
    for(s in 1:length(fdr_thresh)){
      abline(h = fdr_thresh[s], col = "blue", lty = s)
    }
    
    leg = c("FDR", paste0("FDR ", fdr_cl * 100, "% CI"), paste0("Genome-wide sig. ", geno_sig), 
            paste0("FDR = ", fdr_thresh))
    leg_lty = c(1, 3, seq(1, length(geno_sig), 1), seq(1, length(fdr_thresh), 1))
    leg_col = c("black", "black", rep("red", length(geno_sig)), rep("blue", length(fdr_thresh)))
    legend("topright", legend = leg, lty = leg_lty, col = leg_col)
  }
  
  fdr_thresh = sort(fdr_thresh, decreasing = FALSE)
  fdr_thresh_fun = function(x){min(x_vals[which(fdr_ci_hi <= x)])}
  fdr_sum = sapply(fdr_thresh, fdr_thresh_fun)
  names(fdr_sum) = paste0("FDR ", fdr_thresh * 100, "%")
  
  full_sum = signif(c(sl, fdr_sum), 3)
  return(full_sum)
}


