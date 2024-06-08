#文章复现
#安装"TwoSampleMR"包，若已安装，则不需要跑这两行代码
install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR")



#加载包
library(TwoSampleMR)
exp_dat <- read_exposure_data(
  filename = "smkinit.txt",
  sep = "\t",
  snp_col = "rsID",
  beta_col = "Variant Effect Size",
  se_col = "Standard Error",
  effect_allele_col = "Alternate Allele",
  other_allele_col = "Reference Allele",
  eaf_col = "Alternate Allele Frequency",
  pval_col = "P-value",
  samplesize_col = "N"
)

View(exp_dat)


out_dat <- extract_outcome_data(
  snps = exp_dat$SNP,
  outcomes = 'ebi-a-GCST90020053'
)

write.csv(out_dat,file="frailty.csv")

#协调数据
dat <- harmonise_data(
  exposure_dat =  exp_dat, 
  outcome_dat = out_dat,
  action = 2
  )

###计算效应量
res<-mr(dat)
res


#异质性
mr_heterogeneity(dat)
#跑离群值，这行代码时间比较久（迭代次数越大，时间越久）
run_mr_presso(dat,NbDistribution = 1000)
#####水平多效性检验
mr_pleiotropy_test(dat)

##散点图
res <- mr(dat)
p1 <- mr_scatter_plot(res, dat)
p1[[1]]
library(ggplot2)
ggsave(p1[[1]], file="smoking-frailty.tiff", width=7, height=7)



