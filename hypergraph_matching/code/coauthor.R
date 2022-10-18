# Coauthorship data analysis for hypergraph matching
# Oct 18 another analysis idea

# dependencies 

library(data.table)
library(tidyr)
library(gtools) 

setwd("/Volumes/GoogleDrive/My Drive/projects/tensor_matching/code/")
source("functions_mm.R")
setwd("/Volumes/GoogleDrive/My Drive/projects/tensor_matching/code/real_data/")

##### load coauthor data

# # node - author name
# node_label = fread("coauth-DBLP/coauth-DBLP-node-labels.txt", sep = ",", header = F)
# colnames(node_label) = "node_label"
# node_label = separate(node_label, "node_label", into = c("node", "label"), sep = "^\\S*\\K\\s+")

# number of authors in the pub
nverts = fread("coauth-DBLP/coauth-DBLP-nverts.txt", header = F)
nverts = nverts$V1

# timestamps
timestamps = fread("coauth-DBLP/coauth-DBLP-times.txt", header = F)
timestamps = timestamps$V1

# simplex - authors in one pub
simplices = fread("coauth-DBLP/coauth-DBLP-simplices.txt", header = F)
simplices = simplices$V1
simplex_list = split(simplices, as.factor( rep(1:length(nverts), times = nverts) ) )
remove(simplices)


###### subsample data

pub_year = c(2010,2011,2012,2013, 2014,2015,2016, 2017)

pub_ind1 = timestamps %in% pub_year[1]
pub_ind2 = timestamps %in% pub_year[2]
pub_ind3 = timestamps %in% pub_year[3]
pub_ind4 = timestamps %in% pub_year[4]
pub_ind5 = timestamps %in% pub_year[5]
pub_ind6 = timestamps %in% pub_year[6]
pub_ind7 = timestamps %in% pub_year[7]
pub_ind8 = timestamps %in% pub_year[8]

# 2010 is the baseline, correlation between 201x and 2010 will decrease as x goes larger

##### select pub with >=3 authors
# ind1 = pub_ind1 & nverts >= 3
# ind2 = pub_ind2 & nverts >= 3

simplex1 = simplex_list[pub_ind1 & nverts >= 3]
simplex2 = simplex_list[pub_ind2 & nverts >= 3]
simplex3 = simplex_list[pub_ind3 & nverts >= 3]
simplex4 = simplex_list[pub_ind4 & nverts >= 3]
simplex5 = simplex_list[pub_ind5 & nverts >= 3]
simplex6 = simplex_list[pub_ind6 & nverts >= 3]
simplex7 = simplex_list[pub_ind7 & nverts >= 3]
simplex8 = simplex_list[pub_ind8 & nverts >= 3]

freq_author1 = as.data.frame(table(unlist(simplex1)))
freq_author2 = as.data.frame(table(unlist(simplex2)))
freq_author3 = as.data.frame(table(unlist(simplex3)))
freq_author4 = as.data.frame(table(unlist(simplex4)))
freq_author5 = as.data.frame(table(unlist(simplex5)))
freq_author6 = as.data.frame(table(unlist(simplex6)))
freq_author7 = as.data.frame(table(unlist(simplex7)))
freq_author8 = as.data.frame(table(unlist(simplex8)))

##### select highly active authors
author1 = unique(as.numeric(freq_author1[freq_author1$Freq >= 2, ]$Var1))
author2 = unique(as.numeric(freq_author2[freq_author2$Freq >= 2, ]$Var1))
author3 = unique(as.numeric(freq_author3[freq_author3$Freq >= 2, ]$Var1))
author4 = unique(as.numeric(freq_author4[freq_author4$Freq >= 2, ]$Var1))
author5 = unique(as.numeric(freq_author5[freq_author5$Freq >= 2, ]$Var1))
author6 = unique(as.numeric(freq_author6[freq_author6$Freq >= 2, ]$Var1))
author7 = unique(as.numeric(freq_author7[freq_author7$Freq >= 2, ]$Var1))
author8 = unique(as.numeric(freq_author8[freq_author8$Freq >= 2, ]$Var1))

# highly active authors
# active_author = Reduce(intersect, list(author1, author2))
# active_author = Reduce(intersect, list(author1, author2, author3, author4))
active_author = Reduce(intersect, list(author1, author2, author3, author4,author5, author6))
# active_author = Reduce(intersect, list(author1, author2, author3, author4,author5, author6, author7, author8))
#active_author = intersect(intersect(intersect(author1, author2), author3), author4)

# selected authors
# more selected authors with pubs that have at least two selected authors
rel_pub1_ind = unlist(lapply(simplex1, function(x) sum(x %in% active_author) >= 2))
sum(rel_pub1_ind == 1)

rel_pub2_ind = unlist(lapply(simplex2, function(x) sum(x %in% active_author) >= 2))
sum(rel_pub2_ind == 1)

sel_author1 = intersect(unique(unlist(simplex1[rel_pub1_ind])), active_author) 
sel_author2 = intersect(unique(unlist(simplex2[rel_pub2_ind])), active_author) 

sel_author = intersect(sel_author1, sel_author2)

# random choose 50
set.seed(200)
# sel_author = active_author
sel_author = sample(sel_author, 50)



# random selection
# # sel_author = active_author
# set.seed(052)
# sel_author = sample(active_author, 50)

##### find the pubs related to the selected authors 
sel_pub1_ind = unlist(lapply(simplex_list[pub_ind1], function(x) sum(x %in% sel_author) >= 1))
sel_pub2_ind = unlist(lapply(simplex_list[pub_ind2], function(x) sum(x %in% sel_author) >= 1))
sel_pub3_ind = unlist(lapply(simplex_list[pub_ind3], function(x) sum(x %in% sel_author) >= 1))
sel_pub4_ind = unlist(lapply(simplex_list[pub_ind4], function(x) sum(x %in% sel_author) >= 1))
sel_pub5_ind = unlist(lapply(simplex_list[pub_ind5], function(x) sum(x %in% sel_author) >= 1))
sel_pub6_ind = unlist(lapply(simplex_list[pub_ind6], function(x) sum(x %in% sel_author) >= 1))
sel_pub7_ind = unlist(lapply(simplex_list[pub_ind7], function(x) sum(x %in% sel_author) >= 1))
sel_pub8_ind = unlist(lapply(simplex_list[pub_ind8], function(x) sum(x %in% sel_author) >= 1))

sel_simplex1 = simplex_list[pub_ind1][sel_pub1_ind]
sel_simplex2 = simplex_list[pub_ind2][sel_pub2_ind]
sel_simplex3 = simplex_list[pub_ind3][sel_pub3_ind]
sel_simplex4 = simplex_list[pub_ind4][sel_pub4_ind]
sel_simplex5 = simplex_list[pub_ind1][sel_pub5_ind]
sel_simplex6 = simplex_list[pub_ind2][sel_pub6_ind]
sel_simplex7 = simplex_list[pub_ind3][sel_pub7_ind]
sel_simplex8 = simplex_list[pub_ind4][sel_pub8_ind]

sel_simplex_list = list(sel_simplex1, sel_simplex2,sel_simplex3,sel_simplex4,
                        sel_simplex5, sel_simplex6,sel_simplex7,sel_simplex8)

##### construct co-authorship tensor

co_tensor_list = list()

for (j in 1:8) {

  co_tensor = array(0, dim = rep(length(sel_author),3))
  
  cat("co_tensor ", j, "\n")
  simplex_j = sel_simplex_list[[j]]

  # coauthorship tensor 1
  for (i in 1:length(simplex_j)) {
    
    aut_ind = which(sel_author %in% simplex_j[[i]])
    
    if(length(aut_ind) == 1){
      ind_arr = t(as.matrix(rep(aut_ind, 3)))
    }else if(length(aut_ind) == 2){
      ind_arr = permutations(n = 2, r = 3, v = aut_ind, repeats.allowed = T)[2:7,]
    }else if(length(aut_ind) > 2){
      ind_arr = permutations(n = length(aut_ind), r = 3, v = aut_ind, repeats.allowed = F)
    }
    
    # for (j in 1:dim(ind_arr)[1]) {
    #   co_tensor1[ind_arr[j,][1], ind_arr[j,][2], ind_arr[j,][3]] = 1
    # }
    # 
    co_tensor[ind_arr] = co_tensor[ind_arr] + 1
    
  }
  co_tensor_list[[j]] = co_tensor
}

##### test 


######## normalized co-authorship tensor
norm_co_tensor_list = lapply(co_tensor_list, function(x) array( scale(as.vector(x)), dim = dim(x)))

# norm_co_tensor1 = array( scale(as.vector(co_tensor1)), dim = dim(co_tensor1))
# norm_co_tensor2 = array( scale(as.vector(co_tensor2)), dim = dim(co_tensor2))
# norm_co_tensor3 = array( scale(as.vector(co_tensor3)), dim = dim(co_tensor3))
# norm_co_tensor4 = array( scale(as.vector(co_tensor4)), dim = dim(co_tensor4))

pi_star = list(1:length(sel_author), 1:length(sel_author), 1:length(sel_author))


# normalized
co_match12_ini = W1_initial_mm(norm_co_tensor_list[[1]],norm_co_tensor_list[[2]])
matching_error(co_match12_ini, pi_star)

co_match12_re = alter_cleanup_mm(norm_co_tensor_list[[1]],norm_co_tensor_list[[2]],co_match12_ini, 20)
matching_error(co_match12_re, pi_star)

# 
co_match13_ini = W1_initial_mm(norm_co_tensor_list[[1]],norm_co_tensor_list[[8]])
matching_error(co_match13_ini, pi_star)
co_match13_re = alter_cleanup_mm(norm_co_tensor1,norm_co_tensor3,co_match13_ini, 20)
matching_error(co_match13_re, pi_star)



# 
co_match14_ini = W1_initial_mm(norm_co_tensor1,norm_co_tensor4)
matching_error(co_match14_ini, pi_star)
co_match14_re = alter_cleanup_mm(norm_co_tensor1,norm_co_tensor4,co_match14_ini, 20)
matching_error(co_match14_re, pi_star)


cor_vec = c()
tensor1 = norm_co_tensor_list[[1]]
tensor2 = norm_co_tensor_list[[8]]
for (i in 1:dim(norm_co_tensor_list[[1]])[3]) {
  cor_vec = c(cor_vec, cor(as.vector(tensor1[,,i]), as.vector(tensor2[,,i])))
}

###################################



co_tensor1 = array(0, dim = rep(length(sel_author),3))

# coauthorship tensor 1
for (i in 1:length(sel_simplex1)) {
  
  aut_ind = which(sel_author %in% sel_simplex1[[i]])
  
  if(length(aut_ind) == 1){
    ind_arr = t(as.matrix(rep(aut_ind, 3)))
  }else if(length(aut_ind) == 2){
    ind_arr = permutations(n = 2, r = 3, v = aut_ind, repeats.allowed = T)[2:7,]
  }else if(length(aut_ind) > 2){
    ind_arr = permutations(n = length(aut_ind), r = 3, v = aut_ind, repeats.allowed = F)
  }
  
  # for (j in 1:dim(ind_arr)[1]) {
  #   co_tensor1[ind_arr[j,][1], ind_arr[j,][2], ind_arr[j,][3]] = 1
  # }
  # 
  co_tensor1[ind_arr] = co_tensor1[ind_arr] + 1
  
}

# coauthorship tensor 2
co_tensor2 = array(0, dim = rep(length(sel_author),3))
for (i in 1:length(sel_simplex2)) {
  
  aut_ind = which(sel_author %in% sel_simplex2[[i]])
  
  if(length(aut_ind) == 1){
    ind_arr = t(as.matrix(rep(aut_ind, 3)))
  }else if(length(aut_ind) == 2){
    ind_arr = permutations(n = 2, r = 3, v = aut_ind, repeats.allowed = T)[2:7,]
  }else if(length(aut_ind) > 2){
    ind_arr = permutations(n = length(aut_ind), r = 3, v = aut_ind, repeats.allowed = F)
  }
  
  # for (j in 1:dim(ind_arr)[1]) {
  #   co_tensor1[ind_arr[j,][1], ind_arr[j,][2], ind_arr[j,][3]] = 1
  # }
  # 
  co_tensor2[ind_arr] = co_tensor2[ind_arr] + 1
}

# coauthorship tensor 3
co_tensor3 = array(0, dim = rep(length(sel_author),3))
for (i in 1:length(sel_simplex3)) {
  
  aut_ind = which(sel_author %in% sel_simplex3[[i]])
  
  if(length(aut_ind) == 1){
    ind_arr = t(as.matrix(rep(aut_ind, 3)))
  }else if(length(aut_ind) == 2){
    ind_arr = permutations(n = 2, r = 3, v = aut_ind, repeats.allowed = T)[2:7,]
  }else if(length(aut_ind) > 2){
    ind_arr = permutations(n = length(aut_ind), r = 3, v = aut_ind, repeats.allowed = F)
  }
  
  # for (j in 1:dim(ind_arr)[1]) {
  #   co_tensor1[ind_arr[j,][1], ind_arr[j,][2], ind_arr[j,][3]] = 1
  # }
  # 
  co_tensor3[ind_arr] = co_tensor3[ind_arr] + 1
}

# coauthorship tensor 4
co_tensor4 = array(0, dim = rep(length(sel_author),3))
for (i in 1:length(sel_simplex4)) {
  
  aut_ind = which(sel_author %in% sel_simplex4[[i]])
  
  if(length(aut_ind) == 1){
    ind_arr = t(as.matrix(rep(aut_ind, 3)))
  }else if(length(aut_ind) == 2){
    ind_arr = permutations(n = 2, r = 3, v = aut_ind, repeats.allowed = T)[2:7,]
  }else if(length(aut_ind) > 2){
    ind_arr = permutations(n = length(aut_ind), r = 3, v = aut_ind, repeats.allowed = F)
  }
  
  # for (j in 1:dim(ind_arr)[1]) {
  #   co_tensor1[ind_arr[j,][1], ind_arr[j,][2], ind_arr[j,][3]] = 1
  # }
  # 
  co_tensor4[ind_arr] = co_tensor4[ind_arr] + 1
}

