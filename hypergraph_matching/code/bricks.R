### bricks  for Gaussian tensor matching
### Feb 26 Jiaxin

# calcualte the dp distance 
dist_p = function(Ai,Bi,p){
  
  if(p == 1){
    # 1-Wasserstein 
    # Ai, Bi are two vectors of dim d^2
    dp = sum(abs(sort(Ai) - sort(Bi)))/length(Ai)
  }
  if(p == 2){
    cdf1 = ecdf(Ai)
    cdf2 = ecdf(Bi)
    
    knots = c(Ai,Bi)
    
    dp = sqrt(sum( (cdf1(knots) - cdf2(knots))^2 ))
  }
  
  return(dp)
}

# accuracy measurement
match_error = function(m1, m2){
  # m1, m2 two-column matrices of the same dimension
  # m2 true
  
  d = dim(m2)[1]
  
  if(!perfect_matching(m2)){
    warning("m2 should be the true one-to-one matching!" )
    return()
  }
  
  err = 0
  
  m2 = m2[order(m2[,2]),]
  for (i in 1:d) {
    
    tm = m2[i,1]
    m1_cut = m1[m1[,2] == i, 1]
    
    err = err + sum(m1_cut!= tm)
  }
  
  err = err/d
  
  # m1 = m1[order(m1[,2]),]
  # m2 = m2[order(m2[,2]),]
  
  return(err)
}


######## functions for clean_seed
# judge whether it is a perfect matching
perfect_matching <- function(pi) {
  # pi should be 2 column matrix
  
  if(is.null(dim(pi))){
    return(TRUE)
  }
  
  n = dim(pi)[1]
  
  if(length(unique(pi[,1]))!= n |length(unique(pi[,2]))!= n){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

extract_unique_pair = function(pi, margin){
  
  good_name = sort(unique(pi[,margin]))[table(pi[,margin]) == 1]
  good_pairs = pi[pi[,margin] %in% good_name,]
  
  good_pairs = matrix(good_pairs, ncol = 2)
  
  return(good_pairs)
}

# determin_margin = function(pi){
#   if( length(unique(pi[,1])) >= length(unique(pi[,2]))){
#     margin = 2
#   }else{
#     margin = 1
#   }
#   
#   return(margin)
# }


# clean the repeated pairs in the seed
clean_seed = function(seed_set, d_collect){
  # preserve as much pairs as we can
  
  s = dim(seed_set)[1]
  
  if(perfect_matching(seed_set)){
    return(seed_set)
  }
  
  # choose the main column
  # margin = 1, 3 - margin = 2; 
  # margin = 2, 3 - margin = 1
  if( length(unique(seed_set[,1])) >= length(unique(seed_set[,2]))){
    margin = 2
  }else{
    margin = 1
  }
  
  keep = extract_unique_pair(seed_set,margin)
  remain = seed_set[! seed_set[,margin] %in% keep[,margin],]
  # count = 0
  while (dim(matrix(remain, ncol = 2))[1] >= 1) { # until we deal with all the pairs in remain
    # count = count +1
    # cat(count)
    # if(count == 10){
    #   print(remain)
    #   print(margin)
    #   print(keep)
    # }
    
    # delete the pairs in remain that contradict to keep
    if(! is.null(dim(remain))){
      remain = remain[! remain[,(3 - margin)] %in% keep[,(3 - margin)],]
    }
  
    if(perfect_matching(remain)){
      keep = rbind(keep, remain)
      break
    }
    
    # extract the unique pairs in the current remain 
    keep = rbind(keep,extract_unique_pair(remain,margin))
    remain = remain[! remain[,margin] %in% keep[,margin],]
    
    if(perfect_matching(remain)){
      keep = rbind(keep, remain)
      break
    }
    
    # remove the duplication for the first pair in remain

    remain = matrix(remain, ncol = 2)
    ind = remain[1,margin]
    cut_ind = which(remain[,margin] == ind )
    choose_ind = cut_ind[which.min(d_collect[remain[cut_ind,]])]
    
    keep = rbind(keep, remain[choose_ind,])
    remain = remain[! remain[,margin] == ind,]

  }
  
  rownames(keep) = NULL
  
  return(keep)
  
}

