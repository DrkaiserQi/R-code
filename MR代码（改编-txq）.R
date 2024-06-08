###MR Analysis###

#1.运行两样本孟德尔分析R包（TwoSampleMR）
library(TwoSampleMR)
library(data.table)



#2.导入本地（SNPs-暴露）原始数据
raw_dat<- fread("file",sep = "\t",header = T) #file为文件具体名称
#如果所有的数据都在一列中，则需要进行分列并予注释(注意行名，重复就白跑了),*****相关代码如下：
library(tidyverse)
raw_dat1<-raw_dat %>%
       separate(col=c("ukb-d-30880_irnt"),into=c('ES','SE','LP','AF','SS','ID'),sep=":")
           
        
#3.提取P值小于5e-8的SNPs，并进行文件保存 
#将raw_dat矩阵中的p_value转化为数值型数据
exp_dat<- subset(raw_dat, p<5e-8)
write.csv(exp_dat,file = "exposure.csv") ###有些原始数据并无SNP注释，只有Chr-bp数据，需要进行转化处理

#4.重新用适合于MR分析格式的方式读取已经保存的"exposure.csv" SNPs矩阵,命名为exposure_dat
#下方的黄色字符都可以根据原数据修改,***其中样本量数samplesize是非必需的*** 
exposure_dat<- read_exposure_data(
  filename = "exposure.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "Freq1.Hapmap",
  pval_col = "p",
  samplesize_col = "N"
)
#5.剔除连锁不平衡（LD）的SNPs，默认:r2 = 0.001，kb = 10000.
exposure_dat_clumped<- clump_data(exposure_dat)

#6.提取结局数据中相对应暴露中的SNPs数据。***导入本地文件可以用fread函数（data.table R包）
#本地数据（方法一，建议用这种方法）（下方的黄色字符都可以根据原数据修改）
outcome_dat<- fread ("file",sep = "\t",header = T) #file为文件具体名称
out_dat <- read_outcome_data(
  snps = exposure_dat_clumped$SNP,
  filename = "cad.add.160614.website.txt",
  sep = "\t",
  snp_col = "markername",
  beta_col = "beta",
  se_col = "se_dgc",
  effect_allele_col = "effect_allele",
  other_allele_col = "noneffect_allele",
  eaf_col = "effect_allele_freq",
  pval_col = "p_dgc")
#在线数据（方法二，这种方法会自动做“代理”SNP处理）
out_dat <- extract_outcome_data(
  snps = exposure_dat_clumped$SNP,
  outcomes = 'ieu-a-7'
)

#7.协调暴露-结局效应，合并数据
#为什么要协调：1、暴露与SNP的效应方向与结局与SNP的效应方向相反；2、剔除不兼容的SNP（暴露：A/G，结局A/C）
dat <- harmonise_data(
  exposure_dat =  exposure_dat_clumped, 
  outcome_dat = out_dat
)

#8.输出最后用于分析的暴露-结局数据，###并在Excel中计算F值（10）以排除弱工具变量
write.csv(dat, file="dat.csv")

#9.进行MR分析（默认5种方法，以️IVW为主）
res <- mr(dat)
res
#若结局为二分类变量，应计算OR值（OR和b可以进行转化）
generate_odds_ratios(res)

#9*.查看MR总共有多少方法（共18种），可用于后续方法挑选
mr_method_list()
#*比如用IVW，mr-egger两种方法
res1 <- mr(dat,method_list = c("mr_ivw","mr_egger_regression"))
res1

#10.1.散点图 （res为默认5种MR方法的结果散点图,res1为上述两种MR方法的结果散点图，具体根据要求选择）
p1 <- mr_scatter_plot(res, dat)
p1[[1]]
library(ggplot2)
ggsave(p1[[1]], file="res.pdf", width=7, height=7)

#10.2.1.异质性检验(在meta分析中存在异质性需要谨慎解释，而在MR中存在异质性十分常见，不影响结果可靠性，可以用许多原因进行解释)
mr_heterogeneity(dat)
#10.2.2.跑离群值（也可以探究异质性来源），这行代码时间比较久（迭代次数越大，时间越久一般SNP越多，迭代次数要求越多，否则报错）
run_mr_presso(dat,NbDistribution = 1000)
#10.2.3.异质性可视化---漏斗图 (有时候异质性检测P>0.05,但漏斗图显示不对称，写作时应酌情考虑放图)
res_single <- mr_singlesnp(dat)
p2 <- mr_funnel_plot(res_single)
p2[[1]]
#保存漏斗图图片
ggsave(p2[[1]], file="p2.pdf", width=7, height=7)

#10.3.水平多效性检验,若检测出多效性（p<0.05)，建议停止分析（******已经违背MR分析前提******），找其他主题进行分析
mr_pleiotropy_test(dat)

#10.4.留一法分析（leaveoneout）
res_loo <- mr_leaveoneout(dat)
res_loo 
#留一法作图（本质也是森林图） 
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]
#保存图片
ggsave(p3[[1]], file="p3.pdf", width=7, height=7)

#10.5.森林图
p4 <- mr_forest_plot(res_single)
p4[[1]]
#保存图片
ggsave(p4[[1]], file="p4.pdf", width=7, height=7)


