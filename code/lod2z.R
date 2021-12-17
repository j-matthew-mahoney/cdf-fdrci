# Converts LOD scores to z-scores

lod2z = function(l, nind, k = 8){

  # Helper function to convert F-score
  lod2f = function(l){
    f = ((nind - k)/(k-1)) * (10^(2 * l / nind) - 1)
    return(f)
  }
  
  # Helper function to convert to p-value
  lod2p = function(l){
    f = lod2f(l)
    p = 1 - pf(f, df1 = k - 1, df2 = nind - k)
    return(p)
  }
  
    p = lod2p(l)
    z = qnorm(p)
    return(z)
}