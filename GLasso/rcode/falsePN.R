falsePN = function(Omega,Omega_est){
  
  # function calculates the false positive and false negative rates
  # based on the true network Omega and the estimated network Omega_est
  
  # avoid the elements < 1e-16 are counted as nonzero
  Omega_est[abs(Omega_est) <= 1e-16] = 0
  
  falsepositive = Omega_est != 0 & Omega == 0
  num_FP = sum( falsepositive  - diag(diag(falsepositive)) )/2
  # in matlab, sum(sum(mat)) = sumvector(colsum()) = sum all ele,  
  # in r, sum(mat) = sum all ele 
  # diag(diag(fp)) = become_diag_mat( extract_diag_ele(fp) )
  negative = Omega == 0
  num_N = sum( negative - diag(diag(negative))  )/2
  
  if(num_N == 0){
    cat("All elements in Omega are nonzero.\n FP is equal to 0.\n")
    FP = 2
  }else{
    FP = num_FP/num_N # false negative rate
  }
  
  
  falsenegative = Omega_est == 0 & Omega != 0
  num_FP = sum( falsenegative - diag(diag(falsenegative)) )/2
  positive = Omega != 0
  num_P = sum( positive - diag(diag(positive)) )/2
  
  if(num_P == 0){
    cat("Only diagnol elements in Omega are nonzero.\n FN is equal to 0.\n")
    FN = 2
  }else{
    FN = num_FN/num_P # false negative rate
  }
  
  return(list(FP = FP, FN = FN))
  
}

