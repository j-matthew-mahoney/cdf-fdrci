load("~/Documents/Projects/Cube/Cube Data/attie_all_qtl_viewer_v10_04.20.2020.RData")
load("~/Documents/Projects/Cube/Cube Data/dataset.DO.CUBE.multissue.RData")

library(qtl2)
source("~/Documents/Projects/Cube/FDR CI for DO/code/lod2p.R")

# Source linear adjustment function
source("/Users/mahonm/Documents/GitHub/DO-T2D-analysis/code/adjust.R")

# Adjust phenotypes
covar = dataset.invivo$covar.matrix
phys_data = dataset.invivo$data[ , c(1:6, 8, 10:14)]

# Map
ptm = proc.time()
scan1_out = scan1(genoprobs, phys_data[ , "Ins_tAUC"], kinship = K, addcovar = covar)
scan_time = proc.time() - ptm
scan_time

nind = length(intersect(rownames(covar), rownames(phys_data)))
k = 8 # Number of alleles

num_perm = 5
breaks = seq(0, 1, 0.0001)

null_p_list = vector(mode = "list", length = num_perm)

mn_ecdf = matrix(0, nrow = length(breaks) - 1, ncol = 1)
sq_ecdf = matrix(0, nrow = length(breaks) - 1, ncol = 1)

for(n in 1:num_perm){
  
  null_p_list[[n]] = sample_perm_null(genoprobs, pheno, addcovar, kinship = K, intcovar = NULL)
  
  # Compute eCDF data 
  curr_hist = hist(null_p_list[[n]], breaks = breaks, plot = FALSE)
  curr_ecdf = as.matrix(cumsum(curr_hist$counts) / sum(curr_hist$counts))
  mn_ecdf = mn_ecdf + curr_ecdf / num_perm
  sq_ecdf = sq_ecdf + curr_ecdf^2 / (num_perm - 1)
}

m = sum(curr_hist$counts)

numer = sq_ecdf - (num_perm / (num_perm - 1)) * mn_ecdf^2
denom = mn_ecdf * (1 - mn_ecdf)

phi = m * numer / denom

plot(m * numer / denom)

null_count_list = lapply(null_p_list, function(x){hist(x, breaks = seq(0, 1, 0.0001), plot = FALSE)})

find_peaks(scan1_out, map, threshold = 6, sort_by = "lod")

true_hist = hist(lod2p(scan1_out[ , "Ins_tAUC"]), breaks = breaks, plot = FALSE)
true_ecdf = as.matrix(cumsum(true_hist$counts) / sum(true_hist$counts))

var_ind = (1 / mn_ecdf + 1 / (1 - mn_ecdf)) / (m * num_perm) + (1 / true_ecdf + 1 / (1 - true_ecdf)) / m

var_fdr = phi * var_ind


fdr_hat = (mn_ecdf / true_ecdf) * ((1 - true_ecdf) / (1 - mn_ecdf))
fdr_hat = apply(fdr_hat, 1, function(x){min(x, 1)})

cl = 0.95
alpha = 1 - cl

fdr_ci_lo = exp(log(fdr_hat) + qnorm(alpha / 2) * sqrt(var_fdr))
fdr_ci_hi = apply(exp(log(fdr_hat) - qnorm(alpha / 2) * sqrt(var_fdr)), 1, function(x){min(x, 1)})

quartz()
x_vals = -log10(true_hist$mids)
plot(x_vals, fdr_hat)
points(x_vals, fdr_ci_lo)
points(x_vals, fdr_ci_hi)


## Compare to 'fdrci' package
#install.packages("fdrci")
#library(fdrci)

fdrci_out = fdrTbl(lod2p(scan1_out[ , "Ins_tAUC"]), perm.list = null_p_list, pname = "pheno1", ntests = m, cl = 0.95,
                    lowerbound = 0.1, upperbound = 8, incr = 0.1)

quartz()
x_vals = -log10(true_hist$mids)
plot(x_vals, fdr_hat, type = "l", ylim = c(0, 1))
points(x_vals, fdr_ci_lo, type = "l", lty = 2)
points(x_vals, fdr_ci_hi, type = "l", lty = 2)

x_vals = fdrci_out$threshold
points(x_vals, fdrci_out$fdr, col = "red")
points(x_vals, fdrci_out$ul, col = "red", type = "l", lty = 3)
points(x_vals, fdrci_out$ll, col = "red", type = "l", lty = 3)

# ## Is this function sensitive to changing the confidence level
# fdrci_out1 = fdrTbl(lod2p(scan1_out), perm.list = null_p_list, pname = "pheno1", ntests = m, cl = 0.50,
#                     lowerbound = 0.1, upperbound = 8, incr = 0.1)
# 
# fdrci_out2 = fdrTbl(lod2p(scan1_out), perm.list = null_p_list, pname = "pheno1", ntests = m, cl = 0.9,
#                     lowerbound = 0.1, upperbound = 8, incr = 0.1)
# 
# plot(fdrci_out1$ul, fdrci_out2$ul)
# abline(0, 1)

pheno = as.matrix(phys_data[ , "Ins_tAUC"])
ovr_names = intersect(rownames(addcovar), rownames(pheno))

genoprobs = genoprobs[ovr_names, ]

pheno = pheno[ovr_names, , drop = FALSE]
addcovar = covar[ovr_names, , drop = FALSE]


breaks = lod2p(seq(12, 0, -0.01), nind = nind, k = 8)

null_stats = scan1perm_stats(gpr, pheno, addcovar = addcovar, kinship = NULL, intcovar = NULL,
                             nperm = 2, breaks = breaks, scan1_out = scan1_out)

fdr_summarize(null_stats, fdr_cl = 0.5, geno_sig = c(0.01, 0.05), plot = TRUE)


