### HOSVDs
library(rTensor)

renumber <- function(z) {
  z_re <- rep(0, length(z))
  uniq <- unique(z)
  
  for (i in 1:length(uniq)) {
    z_re[z == uniq[i]] <- i
  }
  return(z_re)
}

HOSVD_kmeans = function(Y, r){
  decomp = hosvd(as.tensor(Y), rep(r,3))
  U1 = decomp$U[[1]]
  
  z = kmeans(U1, r)$cluster
  z = renumber(z)
  
  return(z)
}


HOSVD_score = function(Y,r){
  z = rep(0,dim(Y)[1])
  sc = 1:dim(Y)[1]
  
  decomp = hosvd(as.tensor(Y), rep(r,3))
  U1 = decomp$U[[1]]
  
  l2 <- apply(U1, 1, function(x) sum(x^2))
  if (length(which(l2 == 0) > 0)) {
    sc <- which(l2 != 0)
    U1 <- U1[sc, ]
    l2 = l2[sc]
  }
  
  Us = diag(sqrt(l2)^(-1)) %*% U1
  
  z[sc] =  kmeans(Us, r)$cluster
  z[-sc] = sample(unique(z[sc]), length(z[-sc]), replace = T)
  
  z = renumber(z)
  
  
  return(z)
}
