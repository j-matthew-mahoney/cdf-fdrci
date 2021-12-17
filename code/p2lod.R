p2lod = function(p, k = 8, nind){
  f = qf(1 - p, df1 = k - 1, df2 = nind - k)
  l = (nind / 2) * log10(1+f * ((k - 1) / (nind - k)))
  return(l)
}