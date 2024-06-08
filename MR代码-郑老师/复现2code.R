library(TwoSampleMR)
exp_dat <- extract_instruments(outcomes="ieu-a-835")
out_dat <- extract_outcome_data(
  snps = exp_dat$SNP,
  outcomes = 'ukb-a-105'
)

dat <- harmonise_data(
  exposure_dat =  exp_dat, 
  outcome_dat = out_dat,
  action = 1
)

res<-mr(dat)
res

mr_heterogeneity(dat)
mr_pleiotropy_test(dat)


p1 <- mr_scatter_plot(res, dat)
p1[[1]]

library(ggplot2)
ggsave(p1[[1]], file="BMI-RA.tiff", width=7, height=7)

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]
ggsave(p2[[1]], file="BMI-RA-forestplot.tiff", width=7, height=7)


##漏斗图 
res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)
p4[[1]]
##保存图片
ggsave(p4[[1]], file="BMI-RA-funnelplot.tiff", width=7, height=7)


