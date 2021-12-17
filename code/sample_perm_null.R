sample_perm_null = function(genoprobs, pheno, addcovar = NULL, intcovar = NULL, kinship = NULL){
    
  nind = length(pheno)
  k = dim(genoprobs[[1]])[2] # Number of alleles
  
  perm = sample(1:length(pheno), replace = FALSE)
  null_names = rownames(genoprobs[[1]])[perm]
  
  null_geno <- genoprobs
  for(i in 1:length(null_geno)){
    rownames(null_geno[[i]])  <- null_names
  }
  
  scan1_null = scan1(null_geno, pheno, kinship = kinship, addcovar = addcovar, intcovar = intcovar)
  
  null_p = lod2p(scan1_null, nind, k)
  return(null_p)
}


