# Convert LOD score to p-value

lod2p = function(l, nind, k = 8){
  
  # Helper function to convert F-score
  lod2f = function(l){
    f = ((nind - k)/(k-1)) * (10^(2 * l / nind) - 1)
    return(f)
  }
  
  f = lod2f(l)
  p = 1 - pf(f, df1 = k - 1, df2 = nind - k)
  return(p)
}