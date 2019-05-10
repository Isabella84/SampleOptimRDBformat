#######################################################################################################
#######################################################################################################
#######################################################################################################
##
##
##   Biological sampling optimization (Script SampleOptim)
##   Developed by: Patricia Goncalves (IPMA)
##   Version for: Regional DataBase (RDB) exchange format
##
##   Reference: 
##   Goncalves, Patricia 2019. SampleOptim a data analysis R-tool to optimize fish sampling for 
##   biological parameters as input on fish stock assessment.
##
##
##   github link: https://github.com/gonpatricia/SampleOptimRDBformat/3_Simulations_results_data_analysis_RDB.R
## 
##
#######################################################################################################
#######################################################################################################
#######################################################################################################
##Packages: 
library(FSA)
library(FSAdata)
library(nlstools)
library(reshape)
library(ggplot2)
library(ggthemes)
library(cvTools)
library(dplyr)
library("robustbase")
library(MASS)
library(psyphy)
library(boot)
library(RCurl)


###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
############# 3. Simulation results (Plots preparation)
#############   Agregate all the data available from the different simulation runs
#############   For all the matrixes generated from the code

#### Data of age, length by year (from the individuals selected in each simulation run)
dados_lt_age_1<-read.table("dados_lt_age_1.csv",sep=",", header=T)
dados_lt_age_2<-read.table("dados_lt_age_2.csv",sep=",", header=T)
dados_lt_age_3<-read.table("dados_lt_age_3.csv",sep=",", header=T)
dados_lt_age_4<-read.table("dados_lt_age_4.csv",sep=",", header=T)
dados_lt_age_5<-read.table("dados_lt_age_5.csv",sep=",", header=T)
dados_lt_age_6<-read.table("dados_lt_age_6.csv",sep=",", header=T)
dados_lt_age_7<-read.table("dados_lt_age_7.csv",sep=",", header=T)
dados_lt_age_8<-read.table("dados_lt_age_8.csv",sep=",", header=T)
dados_lt_age_9<-read.table("dados_lt_age_9.csv",sep=",", header=T)
dados_lt_age_10<-read.table("dados_lt_age_10.csv",sep=",", header=T)
dados_lt_age_15<-read.table("dados_lt_age_15.csv",sep=",", header=T)
dados_lt_age_20<-read.table("dados_lt_age_20.csv",sep=",", header=T)

data_selected_lt_age<-rbind(dados_lt_age_1,dados_lt_age_2,dados_lt_age_3,dados_lt_age_4,dados_lt_age_5,dados_lt_age_6,dados_lt_age_7,dados_lt_age_8,dados_lt_age_9,dados_lt_age_10,dados_lt_age_15,dados_lt_age_20)
summary(data_selected_lt_age)


############################################################################
####   von Bertalanffy growth model parameters by year and simulation run
############################################################################
#############################################################################
results_simulvbgm_1<-read.table("results_simulvbgm_1.csv",sep=",", header=T)
results_simulvbgm_2<-read.table("results_simulvbgm_2.csv",sep=",", header=T)
results_simulvbgm_3<-read.table("results_simulvbgm_3.csv",sep=",", header=T)
results_simulvbgm_4<-read.table("results_simulvbgm_4.csv",sep=",", header=T)
results_simulvbgm_5<-read.table("results_simulvbgm_5.csv",sep=",", header=T)
results_simulvbgm_6<-read.table("results_simulvbgm_6.csv",sep=",", header=T)
results_simulvbgm_7<-read.table("results_simulvbgm_7.csv",sep=",", header=T)
results_simulvbgm_8<-read.table("results_simulvbgm_8.csv",sep=",", header=T)
results_simulvbgm_9<-read.table("results_simulvbgm_9.csv",sep=",", header=T)
results_simulvbgm_10<-read.table("results_simulvbgm_10.csv",sep=",", header=T)
results_simulvbgm_15<-read.table("results_simulvbgm_15.csv",sep=",", header=T)
results_simulvbgm_20<-read.table("results_simulvbgm_20.csv",sep=",", header=T)

simulvbgm_param_year<-rbind(results_simulvbgm_1,results_simulvbgm_2,results_simulvbgm_3,results_simulvbgm_4,results_simulvbgm_5,results_simulvbgm_6,results_simulvbgm_7,results_simulvbgm_8,results_simulvbgm_9,results_simulvbgm_10,results_simulvbgm_15,results_simulvbgm_20)

####For all the years, to compare the VGBGM parameters between years
####Figure 7 - Summary of parameters by year for the full set of simulations (n=100), for j's (number of selected otoliths)
Fig7_K_VBGM <- file.path(paste("Fig7_K_VBGM_", tm, ".png", sep = ""))
png(file=Fig7_K_VBGM)
fig7_K<-ggplot(simulvbgm_param_year, aes(x=factor(type), y=K, fill=factor(type))) + 
  geom_boxplot(outlier.shape=NA)+xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
fig7_K
dev.off()

Fig7_t0_VBGM <- file.path(paste("Fig7_t0_VBGM_", tm, ".png", sep = ""))
png(file=Fig7_t0_VBGM)
fig7_t0<-ggplot(simulvbgm_param_year, aes(x=factor(type), y=t0,  fill=factor(type))) +
  geom_boxplot(outlier.shape=NA)+xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
fig7_t0
dev.off()

Fig7_Linf_VBGM <- file.path(paste("Fig7_Linf_VBGM_", tm, ".png", sep = ""))
png(file=Fig7_Linf_VBGM)
fig7_Linf<-ggplot(simulvbgm_param_year, aes(x=factor(type), y=Linf, fill=factor(type))) +
  geom_boxplot(outlier.shape=NA)+xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+#ylim(c(20,100))+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
fig7_Linf
dev.off()


####################################################################################
#####  Mean length at age data original and by simulation run
###################################################################################
####################################################################################
table_original_simulsub_mla_1<-read.table("table_original_simulsub_mla_1.csv",sep=",", header=T)
table_original_simulsub_mla_2<-read.table("table_original_simulsub_mla_2.csv",sep=",", header=T)
table_original_simulsub_mla_3<-read.table("table_original_simulsub_mla_3.csv",sep=",", header=T)
table_original_simulsub_mla_4<-read.table("table_original_simulsub_mla_4.csv",sep=",", header=T)
table_original_simulsub_mla_5<-read.table("table_original_simulsub_mla_5.csv",sep=",", header=T)
table_original_simulsub_mla_6<-read.table("table_original_simulsub_mla_6.csv",sep=",", header=T)
table_original_simulsub_mla_7<-read.table("table_original_simulsub_mla_7.csv",sep=",", header=T)
table_original_simulsub_mla_8<-read.table("table_original_simulsub_mla_8.csv",sep=",", header=T)
table_original_simulsub_mla_9<-read.table("table_original_simulsub_mla_9.csv",sep=",", header=T)
table_original_simulsub_mla_10<-read.table("table_original_simulsub_mla_10.csv",sep=",", header=T)
table_original_simulsub_mla_15<-read.table("table_original_simulsub_mla_15.csv",sep=",", header=T)
table_original_simulsub_mla_20<-read.table("table_original_simulsub_mla_20.csv",sep=",", header=T)
table_mla_sims<-rbind(table_original_simulsub_mla_1,table_original_simulsub_mla_2,table_original_simulsub_mla_3,table_original_simulsub_mla_4,table_original_simulsub_mla_5,table_original_simulsub_mla_6,table_original_simulsub_mla_7,table_original_simulsub_mla_8,table_original_simulsub_mla_9,table_original_simulsub_mla_10,table_original_simulsub_mla_15,table_original_simulsub_mla_20)

###############################################################################################################
###### Figure 8 - compare mean length at age (oringinal versus simulation) by year and simulation run
#######
tm<- "T"
years<-unique(table_mla_sims$year)

for(i in years)
{
  fig8_mla<- ggplot(data=subset(table_mla_sims, year==i),aes(x=factor(age), y=m_lt,fill=factor(type))) +
    geom_boxplot(outlier.shape = NA)+xlab("age")+ylab("mean length at age (cm)")+scale_color_brewer(palette = "Paired")+ theme_classic()+
    facet_wrap(~year)+ggtitle(i)+
    theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
          axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
  ggsave(fig8_mla, file=paste0("Fig8_mla_",i,tm,".png"),width=14, height = 10, units="cm")
}


################################

###########################################################################################
### Standard deviation from length at age by simualtion run
#####################################################################
#############################################################################################
table_original_simulsub_sd_1<-read.table("table_original_simulsub_sd_1.csv",sep=",", header=T)
table_original_simulsub_sd_2<-read.table("table_original_simulsub_sd_2.csv",sep=",", header=T)
table_original_simulsub_sd_3<-read.table("table_original_simulsub_sd_3.csv",sep=",", header=T)
table_original_simulsub_sd_4<-read.table("table_original_simulsub_sd_4.csv",sep=",", header=T)
table_original_simulsub_sd_5<-read.table("table_original_simulsub_sd_5.csv",sep=",", header=T)
table_original_simulsub_sd_6<-read.table("table_original_simulsub_sd_6.csv",sep=",", header=T)
table_original_simulsub_sd_7<-read.table("table_original_simulsub_sd_7.csv",sep=",", header=T)
table_original_simulsub_sd_8<-read.table("table_original_simulsub_sd_8.csv",sep=",", header=T)
table_original_simulsub_sd_9<-read.table("table_original_simulsub_sd_9.csv",sep=",", header=T)
table_original_simulsub_sd_10<-read.table("table_original_simulsub_sd_10.csv",sep=",", header=T)
table_original_simulsub_sd_15<-read.table("table_original_simulsub_sd_15.csv",sep=",", header=T)
table_original_simulsub_sd_20<-read.table("table_original_simulsub_sd_20.csv",sep=",", header=T)
table_sdla_sims<-rbind(table_original_simulsub_sd_1,table_original_simulsub_sd_2,table_original_simulsub_sd_3,table_original_simulsub_sd_4,table_original_simulsub_sd_5,table_original_simulsub_sd_6,table_original_simulsub_sd_7,table_original_simulsub_sd_8,table_original_simulsub_sd_9,table_original_simulsub_sd_10,table_original_simulsub_sd_15,table_original_simulsub_sd_20)

###############################################################################################################
###### Figure 9 - compare the sd length at age (oringinal versus simulation) by year and simulation run
#######
#table_sdla_sims<-table_sdla_sims[!is.na(table_sdla_sims$sd_lt),]##remove NAs on data[table_sdla_sims$year==years[[nb]],]

years<-unique(table_sdla_sims$year)

for(i in years)
{
  fig9_sdla<- ggplot(data=subset(table_sdla_sims, year==i),aes(x=factor(age), y=sd_lt,fill=factor(type))) +
    geom_boxplot(outlier.shape =NA)+xlab("age")+ylab("sd length at age (cm)")+scale_color_brewer(palette = "Paired")+ theme_classic()+
    facet_wrap(~year)+ggtitle(i)+
    theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
          axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
  ggsave(fig9_sdla, file=paste0("Fig9_sdla_",i,tm,".png"),width=14, height = 10, units="cm")
}


############################################################################################
###
### Data combine (table_mla_sims, table_sdla_sims)
##################################################################################################

data_mla_sd<- cbind(table_mla_sims, table_sdla_sims)
data_mla_sd<- data_mla_sd[,c(1,2,3,9,10,11,12)]
data_mla_sd<-data_mla_sd[complete.cases(data_mla_sd$sd_lt),]
data_mla_sd$coefv<-data_mla_sd$m_lt/data_mla_sd$sd_lt


years<-unique(data_mla_sd$year)

for(i in years)
{
  fig9_cvla<- ggplot(data=subset(data_mla_sd, year==i),aes(x=factor(age), y=coefv,fill=factor(type))) +
    geom_boxplot(outlier.shape =NA)+xlab("age")+ylab("CV (lenght)")+scale_color_brewer(palette = "Paired")+ theme_classic()+
    facet_wrap(~year)+ggtitle(i)+
    theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
          axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
  ggsave(fig9_cvla, file=paste0("Fig9a_cvla_",i,tm,".png"),width=14, height = 10, units="cm")
}





####################################################################################################
### Stats (mape, rmspe, mspe) from each simualations run and year
######################################################################
######################################################################################################
table_res_stat_1<-read.table("table_res_stat_1.csv",sep=",", header=T)
table_res_stat_2<-read.table("table_res_stat_2.csv",sep=",", header=T)
table_res_stat_3<-read.table("table_res_stat_3.csv",sep=",", header=T)
table_res_stat_4<-read.table("table_res_stat_4.csv",sep=",", header=T)
table_res_stat_5<-read.table("table_res_stat_5.csv",sep=",", header=T)
table_res_stat_6<-read.table("table_res_stat_6.csv",sep=",", header=T)
table_res_stat_7<-read.table("table_res_stat_7.csv",sep=",", header=T)
table_res_stat_8<-read.table("table_res_stat_8.csv",sep=",", header=T)
table_res_stat_9<-read.table("table_res_stat_9.csv",sep=",", header=T)
table_res_stat_10<-read.table("table_res_stat_10.csv",sep=",", header=T)
table_res_stat_15<-read.table("table_res_stat_15.csv",sep=",", header=T)
table_res_stat_20<-read.table("table_res_stat_20.csv",sep=",", header=T)
stats_simul<-rbind(table_res_stat_1,table_res_stat_2,table_res_stat_3,table_res_stat_4,table_res_stat_5,table_res_stat_6,table_res_stat_7,table_res_stat_8,table_res_stat_9,table_res_stat_10,table_res_stat_15,table_res_stat_20)

summary(stats_simul)

###############################################################################################################
###### Figure 10 - compare the stats (mape, mspe, rtmspe) by year and by simulation type (number of otoliths/length class)
#######

Fig10_mspe <- file.path(paste("Fig10_mspe_", tm, ".png", sep = ""))
png(file=Fig10_mspe)
fig10_mspe<-ggplot(stats_simul, aes(x=factor(type), y=mspe, group=1)) + geom_step()+
  #geom_point(size=2)+
  xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
fig10_mspe
dev.off()

Fig10_mape <- file.path(paste("Fig10_mape_", tm, ".png", sep = ""))
png(file=Fig10_mape)
fig10_mape<-ggplot(stats_simul, aes(x=factor(type), y=mape, group=1)) + geom_step()+
  xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
fig10_mape
dev.off()

Fig10_rtmspe <- file.path(paste("Fig10_rtmspe_", tm, ".png", sep = ""))
png(file=Fig10_rtmspe)
fig10_rtmspe<-ggplot(stats_simul, aes(x=factor(type), y=rtmspe, group=1)) + geom_step()+
  xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
fig10_rtmspe
dev.off()


###### Normalized data stats
library(dplyr)

##Mean by stats by year
mean_mspe <- stats_simul %>% group_by(year) %>% summarize(Mean = mean(mspe, na.rm=TRUE))
mean_mape <- stats_simul %>% group_by(year) %>% summarize(Mean = mean(mape, na.rm=TRUE))
mean_rtmspe <- stats_simul %>% group_by(year) %>% summarize(Mean = mean(rtmspe, na.rm=TRUE))

### Normalization of data by stats and year (max)
nor_mspe<- stats_simul %>% group_by(year) %>% mutate(norm_mspe = mspe/max(mspe))
nor_mape<- stats_simul %>% group_by(year) %>% mutate(norm_mape = mape/max(mape))
nor_rtmspe<- stats_simul %>% group_by(year) %>% mutate(norm_rtmspe = rtmspe/max(rtmspe))


Fig10a_norm_mspe <- file.path(paste("Fig10_norm_mspe_", tm, ".png", sep = ""))
png(file=Fig10a_norm_mspe)
fig10a_mspe<-ggplot(nor_mspe, aes(x=factor(type), y=norm_mspe, group=1)) + geom_step()+
  #geom_point(size=2)+
  xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
fig10a_mspe
dev.off()

Fig10a_norm_mape <- file.path(paste("Fig10_norm_mape_", tm, ".png", sep = ""))
png(file=Fig10a_norm_mape)
fig10a_mape<-ggplot(nor_mape, aes(x=factor(type), y=norm_mape, group=1)) + geom_step(stat="summary", fun.y = mean)+
  xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
fig10a_mape
dev.off()

Fig10a_norm_rtmspe <- file.path(paste("Fig10_norm_rtmspe_", tm, ".png", sep = ""))
png(file=Fig10a_norm_rtmspe)
fig10a_rtmspe<-ggplot(nor_rtmspe, aes(x=factor(type), y=norm_rtmspe, group=1)) + geom_step(stat="summary", fun.y = mean)+
  xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
fig10a_rtmspe
dev.off()

### Normalization of data by stats and year (mean)
nor_mn_mspe<- stats_simul %>% group_by(year) %>% mutate(norm_mn_mspe = mspe/mean(mspe))
nor_mn_mape<- stats_simul %>% group_by(year) %>% mutate(norm_mn_mape = mape/mean(mape))
nor_mn_rtmspe<- stats_simul %>% group_by(year) %>% mutate(norm_mn_rtmspe = rtmspe/mean(rtmspe))


Fig10a_norm_mn_mspe <- file.path(paste("Fig10_norm_mn_mspe_", tm, ".png", sep = ""))
png(file=Fig10a_norm_mn_mspe)
fig10a_mn_mspe<-ggplot(nor_mn_mspe, aes(x=factor(type), y=norm_mn_mspe, group=1)) + geom_step()+
  #geom_point(size=2)+
  xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
fig10a_mn_mspe
dev.off()

Fig10a_norm_mn_mape <- file.path(paste("Fig10_norm_mn_mape_", tm, ".png", sep = ""))
png(file=Fig10a_norm_mn_mape)
fig10a_mn_mape<-ggplot(nor_mn_mape, aes(x=factor(type), y=norm_mn_mape, group=1)) + geom_step(stat="summary", fun.y = mean)+
  xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
fig10a_mn_mape
dev.off()

Fig10a_norm_mn_rtmspe <- file.path(paste("Fig10_norm_mn_rtmspe_", tm, ".png", sep = ""))
png(file=Fig10a_norm_mn_rtmspe)
fig10a_mn_rtmspe<-ggplot(nor_mn_rtmspe, aes(x=factor(type), y=norm_mn_rtmspe, group=1)) + geom_step(stat="summary", fun.y = mean)+
  xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
fig10a_mn_rtmspe
dev.off()



###################################################################################################################
####### Comparison of von Bertalanffy growth model parameters estimated by year
###### and for each j (number of otoliths by length class) the result of the 100 simulations aggregated
###################################################################################################################
Fig11_K_VBGM <- file.path(paste("Fig11_K_VBGM_", tm, ".png", sep = ""))
png(file=Fig11_K_VBGM)
fig11_K<-ggplot(stats_simul, aes(x=factor(type), y=k)) + 
  geom_point(size=2,colour="green")+xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
fig11_K
dev.off()

Fig11_t0_VBGM <- file.path(paste("Fig11_t0_VBGM_", tm, ".png", sep = ""))
png(file=Fig11_t0_VBGM)
fig11_t0<-ggplot(stats_simul, aes(x=factor(type), y=t0)) +
  geom_point(size=2,colour="green")+xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
fig11_t0
dev.off()

Fig11_Linf_VBGM <- file.path(paste("Fig11_Linf_VBGM_", tm, ".png", sep = ""))
png(file=Fig11_Linf_VBGM)
fig11_Linf<-ggplot(stats_simul, aes(x=factor(type), y=Linf)) +
  geom_point(size=2,colour="green")+xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+#ylim(c(20,100))+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
fig11_Linf
dev.off()



################################################################################################################
###    Based on the VBGM by sim and year, the length for each selected fish was predicted based on age data
###########################################################################################################
#################################################################################################################
vb_predict_melt_1<-read.table("vb_predict_melt_1.csv",sep=",", header=T)
vb_predict_melt_2<-read.table("vb_predict_melt_2.csv",sep=",", header=T)
vb_predict_melt_3<-read.table("vb_predict_melt_3.csv",sep=",", header=T)
vb_predict_melt_4<-read.table("vb_predict_melt_4.csv",sep=",", header=T)
vb_predict_melt_5<-read.table("vb_predict_melt_5.csv",sep=",", header=T)
vb_predict_melt_6<-read.table("vb_predict_melt_6.csv",sep=",", header=T)
vb_predict_melt_7<-read.table("vb_predict_melt_7.csv",sep=",", header=T)
vb_predict_melt_8<-read.table("vb_predict_melt_8.csv",sep=",", header=T)
vb_predict_melt_9<-read.table("vb_predict_melt_9.csv",sep=",", header=T)
vb_predict_melt_10<-read.table("vb_predict_melt_10.csv",sep=",", header=T)
vb_predict_melt_15<-read.table("vb_predict_melt_15.csv",sep=",", header=T)
vb_predict_melt_20<-read.table("vb_predict_melt_20.csv",sep=",", header=T)
lpredict_vbsim<-rbind(vb_predict_melt_1,vb_predict_melt_2,vb_predict_melt_3,vb_predict_melt_4,vb_predict_melt_5,vb_predict_melt_6,vb_predict_melt_7,vb_predict_melt_8,vb_predict_melt_9,vb_predict_melt_10,vb_predict_melt_15,vb_predict_melt_20)


###############################################################################################################
###### Figure 12 - compare the stats (mape, mspe, rtmspe) by year and by simulation type (number of otoliths/length class)
#######
Fig12_predictLt <- file.path(paste("Fig12_predictLt_", tm, ".png", sep = ""))
png(file=Fig12_predictLt)
fig12_predictLt<-ggplot(data=lpredict_vbsim, aes(x=ID_ind, y=pred_lt,colour=type)) + 
  geom_line()+xlab("ID_ind")+ylab("Predicted length")+theme_classic()+
  facet_wrap(~ano)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
fig12_predictLt
dev.off()


Fig12a_predictLt <- file.path(paste("Fig12a_predictLt_", tm, ".png", sep = ""))
png(file=Fig12a_predictLt)
fig12a_predictLt<-ggplot(data=lpredict_vbsim, aes(x=ID_ind, y=pred_lt,colour=factor(ano))) + 
  geom_point()+xlab("ID_ind")+ylab("Predicted length")+theme_classic()+
  facet_wrap(~factor(type))+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
fig12a_predictLt
dev.off()



################################################################################################################################
###############################################################################################################################
### Data from biological samples 
###############################################################################################################################
###############################################################################################################################
dados_bio_1<-read.table("dados_bio_1.csv",sep=",", header=T)
dados_bio_2<-read.table("dados_bio_2.csv",sep=",", header=T)
dados_bio_3<-read.table("dados_bio_3.csv",sep=",", header=T)
dados_bio_4<-read.table("dados_bio_4.csv",sep=",", header=T)
dados_bio_5<-read.table("dados_bio_5.csv",sep=",", header=T)
dados_bio_6<-read.table("dados_bio_6.csv",sep=",", header=T)
dados_bio_7<-read.table("dados_bio_7.csv",sep=",", header=T)
dados_bio_8<-read.table("dados_bio_8.csv",sep=",", header=T)
dados_bio_9<-read.table("dados_bio_9.csv",sep=",", header=T)
dados_bio_10<-read.table("dados_bio_10.csv",sep=",", header=T)
dados_bio_15<-read.table("dados_bio_15.csv",sep=",", header=T)
dados_bio_20<-read.table("dados_bio_20.csv",sep=",", header=T)
dados_bio_simul<-rbind(dados_bio_1,dados_bio_2,dados_bio_3,dados_bio_4,dados_bio_5,dados_bio_6,dados_bio_7,dados_bio_8,dados_bio_9,dados_bio_10,dados_bio_15,dados_bio_20)


####################################################################################################################
####################################################################################################################
##### Data from the maturity ogive model adjustment
###################################################################################################################
####################################################################################################################
dados_mo_1<-read.table("table_res_mo_1.csv",sep=",", header=T)
dados_mo_2<-read.table("table_res_mo_2.csv",sep=",", header=T)
dados_mo_3<-read.table("table_res_mo_3.csv",sep=",", header=T)
dados_mo_4<-read.table("table_res_mo_4.csv",sep=",", header=T)
dados_mo_5<-read.table("table_res_mo_5.csv",sep=",", header=T)
dados_mo_6<-read.table("table_res_mo_6.csv",sep=",", header=T)
dados_mo_7<-read.table("table_res_mo_7.csv",sep=",", header=T)
dados_mo_8<-read.table("table_res_mo_8.csv",sep=",", header=T)
dados_mo_9<-read.table("table_res_mo_9.csv",sep=",", header=T)
dados_mo_10<-read.table("table_res_mo_10.csv",sep=",", header=T)
dados_mo_15<-read.table("table_res_mo_15.csv",sep=",", header=T)
dados_mo_20<-read.table("table_res_mo_20.csv",sep=",", header=T)
dados_mo_simul<-rbind(dados_mo_1,dados_mo_2,dados_mo_3,dados_mo_4,dados_mo_5,dados_mo_6,dados_mo_7,dados_mo_8,dados_mo_9,dados_mo_10,dados_mo_15,dados_mo_20)
library(tidyr)
dados_mo_re <- dados_mo_simul %>% gather(type_mat, L_mean, L25:L75)


years<- sort(unique(dados_mo_re$year)) ##list of years on the samples data

for(bb in 1:length(years))
{
  Figmo <- file.path(paste("FigmoLt_", years[bb], ".png", sep = ""))
  png(file=Figmo)
  plotmo<-
    ggplot(dados_mo_re[dados_mo_re$year==years[bb],], aes(x=factor(type),y=L_mean, fill=type_mat))+xlab("Number of otoliths selected by length class (cm)")+ ylab("Length")+
    geom_boxplot()+ theme_classic() +
    theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
          axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
  print(plotmo)
  dev.off()
}

###Sem outliers

for(bb in 1:length(years))
{
  Figmo_sna <- file.path(paste("Figmo_snaLt_", years[bb], ".png", sep = ""))
  png(file=Figmo_sna)
  plotmo_sna<-
    ggplot(dados_mo_re[dados_mo_re$year==years[bb],], aes(x=factor(type),y=L_mean, fill=type_mat))+xlab("Number of otoliths selected by length class (cm)")+ ylab("Length")+
    geom_boxplot(outlier.shape = NA)+ theme_classic() +
    theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
          axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
  print(plotmo_sna)
  dev.off()
}


##########################################################################################################
###########################################################################################################
#############################     END CODE ;)    ###########################################################
###########################################################################################################
###########################################################################################################
