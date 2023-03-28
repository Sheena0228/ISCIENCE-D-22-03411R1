setwd("/Users/yfl_ex3/Desktop/KTL_RNA seq")
data_full = read.csv("FPKM.csv")
dim(data_full)

sample_var = colnames(data_full[, 5:16])
DAC_dummy = c(0,0, 0, 1, 1, 1, 0,0,0,1,1,1)
DEX_dummy = c(0,0,0,0,0,0,1,1,1,1,1,1)
sample_data =t(rbind(DAC_dummy, DEX_dummy)) 
rownames(sample_data) = sample_var
sample_data = as.data.frame((sample_data))

data = as.matrix(data_full[, 5:16])
mod = model.matrix(~as.factor(sample_data$DAC_dummy)*as.factor(sample_data$DEX_dummy))
fit = lm.fit(mod, t(data))
cof = fit$coefficients
cof[, 1]
dim(cof)
max(cof[4, ])
min(cof[4, ])

interact = as.array(cof[4, ])
data_inter = cbind(data_full, interact)
data_inter[1, ]

write.csv(data_inter, file = "KTL_seq.csv")

p_val = vector()

for (i in 1:dim(data)[1]){
  
  fit_single = lm(data[i, ] ~sample_data$DAC_dummy*sample_data$DEX_dummy)
  p_val = rbind(p_val, tidy(fit_single)[4, 5])
}


write.csv(p_val, file = "KTL_pval2.csv")
