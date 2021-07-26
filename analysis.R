#
library(data.table)
library(haven)
library(gee)
library(gridExtra)
library(clusrank)
library(ggpubr)
library(perm)

source("geesmv/R/GEE.var.md.R")

# ---- read in the data 
# treatment dates/changes for participants: 
treatment_dates = fread("Tx_dates.txt",sep=",",data.table = F)
treatment_dates$Name = treatment_dates$num
# vitamin C breaking strength for participants:
vitC_breakingstrength = fread("BreakingStrength.txt",sep=",",data.table = F)
vitC_breakingstrength$name = vitC_breakingstrength$num
# ---------------
treatment_dates_reduced = NULL
for(person in unique(treatment_dates$Name)){
  df_subset = treatment_dates[treatment_dates$Name == person,]
  # find start of deprivation
  start_Row = which(df_subset$Deprivation == 1)[1]
  df_subset = df_subset[start_Row:nrow(df_subset),]
  # find where Vit C amount changed
  df_subset$shifted_AA_dose = c(0,df_subset$AA_dose_numeric[1:(nrow(df_subset)-1)])
  df_subset$date = paste(df_subset$day,df_subset$month,df_subset$year,sep="/")
  diff_dose = which(df_subset$shifted_AA_dose != df_subset$AA_dose_numeric)
  diff_dose = c(1,diff_dose,nrow(df_subset))
  differential_dose = df_subset[diff_dose,]
  start_date = as.Date(differential_dose$date[1],format = "%d/%m/%Y")
  differential_dose$days  = c(0,
                              as.numeric(as.Date(differential_dose$date[2:nrow(differential_dose)],
                                                 format = "%d/%m/%Y") - start_date))

  treatment_dates_reduced = rbind(treatment_dates_reduced,
                                  data.frame(Name=person,
                                           AAdose = differential_dose$AA_dose_numeric,
                                           Days = differential_dose$days))
  }

# treatment_dates_reduced: data frame in which for each individual -- day of vitC dose changes 
# (since) start of experiment  + AA dose

expandedDose=NULL
for(person in unique(treatment_dates_reduced$Name)){
  df_subset = treatment_dates_reduced[treatment_dates_reduced$Name == person,]
  vec_days = seq(1:df_subset$Days[nrow(df_subset)])
  vec_doses = NULL
  for(i in 1:(nrow(df_subset)-1)){
    vec_doses = c(vec_doses,
                  rep(df_subset$AAdose[i],df_subset$Days[i+1]-df_subset$Days[i]))
  }
  expandedDose = rbind(expandedDose,
                     data.frame(Name=person,Dose=vec_doses,Days=vec_days))
}

# ------------------------------------------#
average_prior = NULL
days_before = c(10,30,60,90,120,150,180,240)
for(person in unique(vitC_breakingstrength$name)){
  vitC_subset = vitC_breakingstrength[vitC_breakingstrength$name == person,]
  df_Date = treatment_dates[treatment_dates$Name == person,]
  start_Row = which(df_Date$Deprivation == 1)[1]
  df_Date = df_Date[start_Row:nrow(df_Date),]
  df_Date$date = as.Date(paste(df_Date$day,df_Date$month,df_Date$year,sep="/"),
                           format = "%d/%m/%Y")
  start_Date = min(df_Date$date)
  vitC_subset$time_Since_Entry = as.numeric(as.Date(vitC_subset$datescar_R) - start_Date)+1
  df_Sub = expandedDose[expandedDose$Name == person,]
  
  vitC_subset$sum_Dose = rep(0,nrow(vitC_subset)) 
  sum_Dose_Prev = matrix(0,nrow=nrow(vitC_subset),ncol=length(days_before))
  sum_Dose_Time = matrix(0,nrow=nrow(vitC_subset),ncol=length(days_before))
  
  for(times in 1:nrow(vitC_subset)){
    time_val = vitC_subset$time_Since_Entry[times]
    vitC_subset$sum_Dose[times] = sum(df_Sub$Dose[1:(time_val)])
    for(vals in 1:length(days_before)){
      val_int = max(1,time_val-days_before[vals])
      sum_Dose_Prev[times,vals] = sum(df_Sub$Dose[(val_int+1):(time_val)])
      sum_Dose_Time[times,vals] = min(days_before[vals],vitC_subset$time_Since_Entry[times])
    }

  }
  vitC_subset$average  = vitC_subset$sum_Dose/vitC_subset$time_Since_Entry
  sum_Dose_Prev_Avg = sum_Dose_Prev/sum_Dose_Time
  colnames(sum_Dose_Prev_Avg) =  c("10day_avg","30day_avg",
                                   "60day_avg","90day_avg","120day_avg",
                                   "150day_avg","180day_avg","240day_avg")
  vitC_subset = cbind(vitC_subset,sum_Dose_Prev_Avg)
  average_prior = rbind(average_prior,data.frame(vitC_subset))
}

# ------------- #
uniquenames = unique(average_prior$name)
name_id = NULL
for(i in 1:nrow(vitC_breakingstrength)){
  name_id = c(name_id,which(uniquenames == average_prior$name[i]))
}
average_prior$name_id = name_id
#---------change scar strength measurement to units of newton:
average_prior$oldstrength = average_prior$strength
average_prior$strength = average_prior$oldstrength/1000 * 9.80665
# ------------#
library(geesmv)
Rx_group=c(1,2,3)
avg_alltime = gee(strength~average,id = factor(name),data=average_prior[average_prior$treatment%in% Rx_group,],
                      corstr = "independence")
avg_alltimevcov = GEE.var.md.special(strength~average,id="name_id",
                     data=average_prior[average_prior$treatment%in% Rx_group,],
                     corstr="independence")
  
coefficinets = coef(summary(avg_alltime))
vccov_matrix = avg_alltimevcov$cov.beta # avg_alltime$robust.variance
v_a = vccov_matrix[1,1]
v_b = vccov_matrix[2,2]
cov_ab = vccov_matrix[1,2]
 
# ---- compare to when subset to scars with average daily vit C >= 10----- #

avg_prior_remove10 = average_prior[average_prior$treatment%in% Rx_group & average_prior$average >= 10,]
avg_alltime_remove10 = gee(strength~average,id = factor(name),data=avg_prior_remove10,
                           corstr = "independence")
avg_alltimecov_remove10 = GEE.var.md.special(strength~average,id="name_id", data=avg_prior_remove10,
                                             corstr="independence")
coefficinets_remove10 = coef(summary(avg_alltime_remove10))
vccov_matrix_remove10 = avg_alltimecov_remove10$cov.beta # avg_alltime$robust.variance
v_a_remove10 = vccov_matrix_remove10[1,1]
v_b_remove10 = vccov_matrix_remove10[2,2]
cov_ab_remove10 = vccov_matrix_remove10[1,2]

# ---- compare results----- #

print("All average daily Vit C intakes:")
print(paste(signif(coefficinets["(Intercept)","Estimate"],3)," + ", signif(coefficinets["average","Estimate"],3)," * average Vit C intake",sep=""))
print(paste("Estimated slope: ",signif(coefficinets["average","Estimate"],3),
            " (p: ", signif( 2 * pnorm(abs(coefficinets["average","Estimate"]/sqrt(v_b)), lower.tail = FALSE),3),
            ")",sep=""))

print("Restricting to average daily Vit C >= 10:")
print(paste(signif(coefficinets_remove10["(Intercept)","Estimate"],3),
                     " + ", signif(coefficinets_remove10["average","Estimate"],3),
                     " * average Vit C intake",sep=""))
print(paste("Estimated slope: ",signif(coefficinets_remove10["average","Estimate"],3),
                     " (p: ", signif( 2 * pnorm(abs(coefficinets_remove10["average","Estimate"]/sqrt(v_b_remove10)), lower.tail = FALSE),3),
                     ")",sep=""))
# ---- look at association of scar strength with prior vit C intakes averaged over different time periods ----- #
 
  md.ind.10 <- GEE.var.md(strength~X10day_avg,id="name_id",data=average_prior[average_prior$treatment%in% Rx_group,],corstr="independence")
  md.ind.30 <- GEE.var.md(strength~X30day_avg,id="name_id",data=average_prior[average_prior$treatment%in% Rx_group,],corstr="independence")
  md.ind.60 <- GEE.var.md(strength~X60day_avg,id="name_id",data=average_prior[average_prior$treatment%in% Rx_group,],corstr="independence")
  md.ind.90 <- GEE.var.md(strength~X90day_avg, id="name_id",data=average_prior[average_prior$treatment%in% Rx_group,],corstr="independence")
  md.ind.120 <- GEE.var.md(strength~X120day_avg, id="name_id",data=average_prior[average_prior$treatment%in% Rx_group,],corstr="independence")
  md.ind.150 <- GEE.var.md(strength~X150day_avg, id="name_id",data=average_prior[average_prior$treatment%in% Rx_group,],corstr="independence")
  md.ind.180 <- GEE.var.md(strength~X180day_avg, id="name_id",data=average_prior[average_prior$treatment%in% Rx_group,],corstr="independence")
  md.ind.240 <- GEE.var.md(strength~X240day_avg,id="name_id",data=average_prior[average_prior$treatment%in% Rx_group,],corstr="independence")
  md.ind.all <- GEE.var.md(strength~average,id="name_id",data=average_prior[average_prior$treatment%in% Rx_group,],corstr="independence")
  
  avg_10day = gee(strength~X10day_avg,id = factor(name),data=average_prior[average_prior$treatment%in% Rx_group,],corstr = "independence")
  avg_30day = gee(strength~X30day_avg,id = factor(name),data=average_prior[average_prior$treatment%in% Rx_group,],corstr = "independence")
  avg_60day = gee(strength~X60day_avg,id = factor(name),data=average_prior[average_prior$treatment %in% Rx_group,],corstr = "independence")
  avg_90day = gee(strength~X90day_avg,id = factor(name),data=average_prior[average_prior$treatment%in% Rx_group,],corstr = "independence")
  avg_120day = gee(strength~X120day_avg,id = factor(name),data=average_prior[average_prior$treatment %in% Rx_group,],corstr = "independence")
  avg_150day = gee(strength~X150day_avg,id = factor(name),data=average_prior[average_prior$treatment %in% Rx_group,],corstr = "independence")
  avg_180day = gee(strength~X180day_avg,id = factor(name),data=average_prior[average_prior$treatment %in% Rx_group,],corstr = "independence")
  avg_240day = gee(strength~X240day_avg,id = factor(name),data=average_prior[average_prior$treatment%in% Rx_group,],corstr = "independence")
  

  zscore = data.frame(zcore_robust =   c((summary(avg_10day)$coefficients[,"Estimate"]/sqrt(md.ind.10$cov.beta))[2],
                                        (summary(avg_30day)$coefficients[,"Estimate"]/sqrt(md.ind.30$cov.beta))[2],
                                          (summary(avg_60day)$coefficients[,"Estimate"]/sqrt(md.ind.60$cov.beta))[2],
                                         (summary(avg_90day)$coefficients[,"Estimate"]/sqrt(md.ind.90$cov.beta))[2],
                                         (summary(avg_120day)$coefficients[,"Estimate"]/sqrt(md.ind.120$cov.beta))[2],
                                         (summary(avg_150day)$coefficients[,"Estimate"]/sqrt(md.ind.150$cov.beta))[2],
                                         (summary(avg_180day)$coefficients[,"Estimate"]/sqrt(md.ind.180$cov.beta))[2],
                                         (summary(avg_240day)$coefficients[,"Estimate"]/sqrt(md.ind.240$cov.beta))[2],
                                         (summary(avg_alltime)$coefficients[,"Estimate"]/sqrt(md.ind.all$cov.beta))[2]),
                      val = c(coef(summary(avg_10day))[2,1],
                              coef(summary(avg_30day))[2,1],
                              coef(summary(avg_60day))[2,1],
                              coef(summary(avg_90day))[2,1],
                              coef(summary(avg_120day))[2,1],
                              coef(summary(avg_150day))[2,1],
                              coef(summary(avg_180day))[2,1],
                              coef(summary(avg_240day))[2,1],
                              coef(summary(avg_alltime))[2,1]),
                      var = c("10 day prior", "30 day prior", "60 day prior",
                              "90 day prior", "120 day prior",  "150 day prior","180 day prior", "240 day prior",
                              "Entire study period"))
  zscore$var = factor(zscore$var,levels = c("10 day prior","30 day prior", "60 day prior",
                                            "90 day prior", "120 day prior",
                                            "150 day prior",  "180 day prior","240 day prior",
                                            "Entire study period"))
  zscore$pval = signif(2 * pnorm(abs(zscore$zcore_robust), lower.tail = FALSE),digits = 3)
  print(zscore,row.names = F)
  
  # ------ permutation tests of post-saturation vit C scar strengths
  # two-sample permutation tests (permTS)
  post_sat_data = vitC_breakingstrength[which(vitC_breakingstrength$postsaturation == 1),]
  
  # 1 : 70 mg // 2 : 10 mg // 3 : 0 mg
  #  p_val70v10
  permTS(post_sat_data[post_sat_data$treatment==1,"strength"],
         post_sat_data[post_sat_data$treatment==2,"strength"], 
         alternative = "two.sided", method = "exact.ce", 
         control = permControl(tsmethod = "abs"))
  # p_val70v0
  permTS(post_sat_data[post_sat_data$treatment==1,"strength"],
         post_sat_data[post_sat_data$treatment==3,"strength"], 
         alternative = "two.sided", method = "exact.ce", 
         control = permControl(tsmethod = "abs"))
  # p_val10v0
  permTS(post_sat_data[post_sat_data$treatment==2,"strength"],
         post_sat_data[post_sat_data$treatment==3,"strength"], 
         alternative = "two.sided", method = "exact.ce", 
         control = permControl(tsmethod = "abs"))