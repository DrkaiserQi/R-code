install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)


###读取暴露数据（下面输入文件中具体参数需要根据文件信息修改）
exp_dat <- read_exposure_data(
  filename = "All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "Freq1.Hapmap",
  pval_col = "p",
  samplesize_col = "N"
)

#提取与暴露强相关（P< 5e-8）的SNP，阈值可以适当调整
exp_dat <- exp_dat[exp_dat$pval.exposure < 5e-8,]

#剔除连锁不平衡LD（默认kb = 10000，r2 = 0.001），阈值可以适当调整
exp_dat <- clump_data(exp_dat)

#输出暴露-SNP文件
write.csv(exp_dat, file="exp_dat.csv")

###提取结局数据---第一种方法（自动做“代理”SNP处理）
out_dat <- extract_outcome_data(
  snps = exp_dat$SNP,
  outcomes = 'ieu-a-7'
)
#提取结局数据---第二种方法（本地数据）（下面输入文件中具体参数需要根据文件信息修改）
out_dat <- read_outcome_data(
  snps = exp_dat$SNP,
  filename = "cad.add.160614.website.txt",
  sep = "\t",
  snp_col = "markername",
  beta_col = "beta",
  se_col = "se_dgc",
  effect_allele_col = "effect_allele",
  other_allele_col = "noneffect_allele",
  eaf_col = "effect_allele_freq",
  pval_col = "p_dgc")

###协调暴露-结局效应，合并数据
#为什么要协调：1、暴露与SNP的效应方向与结局与SNP的效应方向相反；2、剔除不兼容的SNP（暴露：A/G，结局A/C）
dat <- harmonise_data(
  exposure_dat =  exp_dat, 
  outcome_dat = out_dat
)
write.csv(dat, file="dat.csv")#输出最后用于分析的暴露-结局数据

###进行MR分析（5种方法，以️IVW为主）
res <- mr(dat)
res
generate_odds_ratios(res)

#查看MR总共有多少方法，可用于后续挑选
mr_method_list()
#比如用IVW，mr-egger两种方法
res1 <- mr(dat,method_list = c("mr_ivw","mr_egger_regression"))
res1

#敏感性分析
#异质性检验
mr_heterogeneity(dat)

#跑离群值，这行代码时间比较久（迭代次数越大，时间越久）
run_mr_presso(dat,NbDistribution = 1000)

#####水平多效性检验
mr_pleiotropy_test(dat)

##获取单独的每个SNP的效应值
res_single <- mr_singlesnp(dat)
res_single

##留一法分析
res_loo <- mr_leaveoneout(dat)
res_loo 

######做图######
##散点图 （res1为为上述两种MR方法的结果散点图，具体可以调整）
p1 <- mr_scatter_plot(res1, dat)
p1[[1]]
library(ggplot2)
ggsave(p1[[1]], file="res.pdf", width=7, height=7)


#森林图
p2 <- mr_forest_plot(res_single)
p2[[1]]
##保存图片
ggsave(p2[[1]], file="p2.pdf", width=7, height=7)


##留一法图（本质也是森林图） 
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]
##保存图片
ggsave(p3[[1]], file="p3.pdf", width=7, height=7)


##漏斗图 
res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)
p4[[1]]
##保存图片
ggsave(p4[[1]], file="p4.pdf", width=7, height=7)

