#### merge the output file in form paste0("s_",s,"_r_",r,"_i_",i,"_n_",n,"_fig5_sample.RData") 
final_sample=final_sd_sample=array(0,dim=c(2,2,5,5)) # signal, rank, info, method
final_sample_cor=final_sd_sample_cor=array(0,dim=c(2,2,5,5)) 

dup = 30
for (s in 1:2) {
    for (r in 1:2) {
        for (i in 1:5){

            #s = 1; r = 1; i = 3
            err_table = cor_table = c() # std, env, mreg, rrr
            for (n in 1:dup){
               # s = 1; r = 1; i = 3; n = 26
                file = paste0("home/jhu267/jcgs-code/Figure5/sample_output_v1/","s_",s,"_r_",r,"_i_",i,"_n_",n,"_fig5_sample.RData")
                load(file)

                err_table = rbind(err_table, result$err)
                cor_table = rbind(cor_table, result$cor)
            }

            final_sample[s,r,i,1:4] = apply(err_table, 2, mean)
            final_sd_sample[s,r,i,1:4] = apply(err_table, 2, sd)

            final_sample_cor[s,r,i,1:4] = apply(cor_table, 2, mean)
            final_sd_sample_cor[s,r,i,1:4] = apply(cor_table, 2, sd)
        }
    }
}

save(final_sample, final_sd_sample, final_sample_cor, final_sd_sample_cor, file = "sample_fig5_v1.RData")


####### merge .mat file in form into one big file paste0("s_",s,"_r_",r,"_i_",i,"_n_",n,"_fig5_sample.mat")

if(!require(rmatio)){ # write .mat file for SupCP
  install.packages("rmatio")
  library(rmatio)
}

tsr_list_s = list(); x_list_s = list(); u_list_s = list();para_list_s = list()
for (s in 1:2) {
    tsr_list_r = list(); x_list_r = list(); u_list_r = list();para_list_r = list()
    for (r in 1:r) {
        tsr_list_i = list(); x_list_i = list(); u_list_i = list();  para_list_i = list()
        for (i in 1:5) {
           tsr_list = list(); x_list = list(); u_list = list();  para_list = list()
           for (n in 1:dup) {
              file = paste0("home/jhu267/jcgs-code/Figure5/sample_output_v1/","s_",s,"_r_",r,"_i_",i,"_n_",n,"_fig5_sample.mat")
              mat_data = read.mat(file)

              tsr_list[[n]] = mat_data$tsr
              x_list[[n]] = mat_data$y
              u_list[[n]] = mat_data$x_true
              para_list[[n]] = mat_data$para
           }
           tsr_list_i[[i]] = tsr_list; x_list_i[[i]] = x_list; u_list_i[[i]] = u_list; para_list_i[[i]] = para_list 
        }
       tsr_list_r[[r]] = tsr_list_i; x_list_r[[r]] = x_list_i; u_list_r[[r]] = u_list_i; para_list_r[[r]] = para_list_i
    }
    tsr_list_s[[s]] = tsr_list_r; x_list_s[[s]] = x_list_r; u_list_s[[s]] = u_list_r; para_list_s[[s]] = para_list_r
}

sample_list = list(tsr = tsr_list_s, y = x_list_s, x_true = u_list_s, para = para_list_s)
write.mat(sample_list,"sample_fig5_v1.mat")


