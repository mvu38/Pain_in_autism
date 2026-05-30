# this code takes extracts from the EMA study of pain experiences in autism conducted 
# by a team of Dr. Katelynn Boerner at BC Children's Hospital Research Institute.
# Ecological momentary assessment (EMA), demographic, and medical data relevant for the
# paper titled
# "Ecological momentary assessment of pain in autistic children, adolescents, and young adults"
# by Dudarev et al.

# loading libraries and functions
library(dplyr)
library(tidyr)
library(psych)
library(lmerTest)
library(sjPlot)
library(ggplot2)
library(effectsize)

# read the data
setwd(getwd())
ema_prompt <- read.csv("ema_byprompt.csv")
demog_med <- read.csv('demog_med.csv')

# treat gender in demog-med as described in the paper
# notice that in raw data sex was coded as 1=male, 2 = female, but
# gender was coded as 1 = girl/woman, 2 = boy/man, 3 + - other.
demog_med$gender <- ifelse(demog_med$demo_sex==1 & demog_med$demo_gender==2,"boy",
                           ifelse(demog_med$demo_sex==2 & demog_med$demo_gender==1,"girl","TGD"))
demog_med$gender <- ifelse(is.na(demog_med$gender),"TGD",demog_med$gender)

###########
# Aggregate EMA data into days per participant
###########
ema_summary <- ema_prompt %>%
  group_by(record_id,redcap_event_name) %>%
  summarise_at(vars(pain_stomach,headache,low_back_pain,limb_pain,pain_other,no_pain,n_of_pains_prompt,
                    pain_intensity_max,pain_interference_max,
                    faint_dizzy,tachycardia,nausea,weakness,migraine,fatigue,other_symptoms,no_symptoms),
               list(max = ~max(.,na.rm=TRUE),
                    mean = ~mean(.,na.rm=TRUE)))

# code pain and symptom presence as two binary variables
ema_summary$pain_bin <- ifelse(ema_summary$no_pain_mean ==1,0,1)
ema_summary$symptoms_bin <- ifelse(ema_summary$no_symptoms_mean ==1,0,1)

# compute n of pains across the day
ema_summary$n_of_pains <- rowSums(ema_summary[,c("pain_stomach_max","headache_max",
                                                 "low_back_pain_max","limb_pain_max",
                                                 "pain_other_max")],na.rm=TRUE)

# compute n of symptoms across the day
ema_summary$n_of_symptoms <- rowSums(ema_summary[,c("faint_dizzy_max","tachycardia_max",
                                                    "nausea_max","weakness_max",
                                                    "migraine_max","fatigue_max","other_symptoms_max")],na.rm=TRUE)
###########
# Aggregate EMA data per participant across days
###########
ema_summary_pp <- ema_summary %>%
  group_by(record_id) %>%
  summarize_at(vars(pain_bin,symptoms_bin,pain_intensity_max_mean,pain_interference_max_mean,
                    n_of_pains_prompt_mean),
               list(sum = ~sum(.,na.rm=TRUE),
                    mean = ~mean(.,na.rm=TRUE),
                    n = ~sum(!is.na(.x))))

ema_summary_pp$proportion_days_with_pain <- ema_summary_pp$pain_bin_sum/ema_summary_pp$pain_bin_n
ema_summary_pp$proportion_days_with_symp <- ema_summary_pp$symptoms_bin_sum/ema_summary_pp$symptoms_bin_n

# add demographics
ema_summary_pp <- merge(ema_summary_pp,demog_med[,c("record_id","meddx_CP","meddx_sum",
                                                    "gender","age")],by="record_id")

# add chronic pain/illness/none
ema_summary_pp$CP_nonCP_healthy <- ifelse(ema_summary_pp$meddx_CP==0,
                                          ifelse(ema_summary_pp$meddx_sum==0,0,1),2)

###########
# COMPLIANCE AND EXCLUSIONS
###########
# summarize compliance (compliance defined as at least half of EMA questions answered)
compliance <- ema_summary %>%
  group_by(record_id) %>%
  summarize_at(vars(redcap_event_name),
               list(n_days_completed = ~sum(!is.na(.x))))

hist(compliance$n_days_completed)

# exclude those with less than 7 days of answers
include <- compliance[compliance$n_days_completed > 6, "record_id"]

ema_summary <- ema_summary[ema_summary$record_id %in% include$record_id, ]
ema_summary_pp <- ema_summary_pp[ema_summary_pp$record_id %in% include$record_id, ]
ema_prompt <- ema_prompt[ema_prompt$record_id %in% include$record_id, ]
demog_med <- demog_med[demog_med$record_id %in% include$record_id, ]

##############
# compliance for report
##############
ema_prompt <- ema_prompt[is.na(ema_prompt$timestamp)==FALSE, ]

# compute compliance
nrow(ema_prompt) / (14 * 3 * length(include$record_id)) * 100 

#-----------------------------------------------------------------------------------------------------
#                               Age and gender on chronic and momentary pain 
#-----------------------------------------------------------------------------------------------------
# code gender
ema_summary_pp$cisboyvsgirl <- ifelse(ema_summary_pp$gender=="boy",1,ifelse(ema_summary_pp$gender=="girl",-1,0))
ema_summary_pp$genderothervsboy <- ifelse(ema_summary_pp$gender =="TGD",1,ifelse(ema_summary_pp$gender=="boy",-1,0))
# centralize age
ema_summary_pp$age_cent <- ema_summary_pp$age - mean(ema_summary_pp$age,na.rm=TRUE)
# impute missing age with sample average
ema_summary_pp$age_cent <- ifelse(is.na(ema_summary_pp$age_cent),0,ema_summary_pp$age_cent)

# test CP by age by gender
chronicpain2 <- glm(meddx_CP ~ age_cent * cisboyvsgirl * genderothervsboy, family=binomial,data=ema_summary_pp)
summary(chronicpain2)
standardize_parameters(chronicpain2,method="basic")

# daily pain regardless of CP
painagegender <- lm(proportion_days_with_pain ~ age_cent * cisboyvsgirl * genderothervsboy, 
                    data=ema_summary_pp)
summary(painagegender)
standardize_parameters(painagegender,method="basic")

# daily pain within no-illness group
painagegender <- lm(proportion_days_with_pain ~ age_cent * cisboyvsgirl * genderothervsboy, 
                    data=ema_summary_pp[ema_summary_pp$CP_nonCP_healthy==0, ])
summary(painagegender)
standardize_parameters(painagegender,method="basic")

# pain intensity
painagegender2 <- lm(pain_intensity_max_mean_mean ~ age_cent * cisboyvsgirl * genderothervsboy, 
                     data=ema_summary_pp[ema_summary_pp$CP_nonCP_healthy==0, ])
summary(painagegender2)

#-----------------------------------------------------------------------------------------------------
#                               Chronic pain vs. chronic illness vs. no illness 
#-----------------------------------------------------------------------------------------------------
# ANOVAS comparing chronic pain, chronic illness, and no illness, were performed in SPSS, using 
# datafile exported here:
write.csv(ema_summary_pp,"pain_acute_chronic.csv")


#-----------------------------------------------------------------------------------------------------
#                                           Control group
#-----------------------------------------------------------------------------------------------------
# for ASD: from ema_summary compute first and second week separately
# read a file with control group data
# create a file with pain frequency and intensity, age, group (ASD, non-ASD), gender.
# 2 measures: pain+symptoms; only pain (careful with 18 yo)

# 18 yo
control_ema <- read.csv('18_yo_daily_ema.csv')
control_demog <- read.csv('18_yo_baselines_youth.csv')

# select and parse data from 18yo
control_ema$physical_symp_bin <- ifelse(control_ema$physical.symp_mean > 0,1,0)

control_ema <- control_ema %>%
  rowwise() %>%
  mutate(pain_clean = max(stomach.ache_mean,headache_mean,lbp_mean,armslegspain_mean,other_mean)) 

# clean out empty rows
control_ema <- control_ema[is.na(control_ema$day_from_start_redcap)==FALSE, ]
control_ema$pain_clean_bin <- ifelse(control_ema$pain_clean > 0,1,0)

control_ema$pain_intensity <- ifelse(control_ema$pain_clean_bin==1,control_ema$symptom_intensity_mean,NA)

control_ema_summary <- control_ema %>%
  group_by(family_id) %>%
  summarize_at(vars(physical_symp_bin,pain_clean_bin,symptom_intensity_mean,pain_intensity,day_from_start_redcap),
               list(mean = ~mean(.,na.rm=TRUE),
                    sum = ~sum(.,na.rm=TRUE),
                    n = length))

control_ema_summary <- merge(control_ema_summary[is.na(control_ema_summary$family_id)==FALSE, 
                                                 c("family_id","physical_symp_bin_sum","pain_clean_bin_sum",
                                                   "symptom_intensity_mean_mean","pain_intensity_mean",
                                                   "day_from_start_redcap_n")],
                             control_demog[,c("family_id","age","gender")],by="family_id") %>%
  rename(record_id = family_id)

control_ema_summary$group <- "control"

hist(control_ema_summary$age)

# RECODE GENDER (in the control sample, coding was different): 1 = woman, 2 = man, 3 = gender-diverse
control_ema_summary$gender <- ifelse(control_ema_summary$gender==1,2,
                                           ifelse(control_ema_summary$gender==2,1,3))
control_ema_summary$gender <- ifelse(is.na(control_ema_summary$gender),3,control_ema_summary$gender)

# in ASD data, parse pains to match the control sample
ema_summary <- ema_summary %>%
  rowwise() %>%
  mutate(
    physical.symp_mean = max(pain_stomach_mean,headache_mean,low_back_pain_mean,limb_pain_mean,pain_other_mean,
                             faint_dizzy_mean,tachycardia_mean,nausea_mean,weakness_mean,other_symptoms_mean),
    pain_clean = max(pain_stomach_mean,headache_mean,low_back_pain_mean,limb_pain_mean,pain_other_mean)
  )

ema_summary$physical_symp_bin <- ifelse(ema_summary$physical.symp_mean > 0,1,0)
ema_summary$pain_clean_bin <- ifelse(ema_summary$pain_clean > 0,1,0)

# in ASD sample, split the data into week 1 vs. week 2
ema_summary$day_from_start_redcap <- as.numeric(
  gsub("ema_day","",gsub("_arm_1","",ema_summary$redcap_event_name)))

ema_summary$week <- ifelse(ema_summary$day_from_start_redcap < 8,1,2)

asd_ema_summary <- ema_summary %>%
  group_by(record_id,week) %>%
  summarize_at(vars(physical_symp_bin,pain_clean_bin,pain_intensity_max_mean,day_from_start_redcap),
               list(mean = ~mean(.,na.rm=TRUE),
                    sum = ~sum(.,na.rm=TRUE),
                    n = ~sum(!is.na(.x))))

# figure out meddx for ASD
asd_ema_summary <- merge(asd_ema_summary[is.na(asd_ema_summary$record_id)==FALSE, 
                                         c("record_id","week","physical_symp_bin_sum","pain_clean_bin_sum",
                                           "pain_intensity_max_mean_mean",
                                           "day_from_start_redcap_n")],
                         demog_med[,c("record_id","age","gender","meddx_sum","meddx_CP")],
                         by="record_id")

asd_ema_summary$group <- "autism"

# reshape ASD for rep.measures
asd_ema_wide <- reshape(asd_ema_summary,direction="wide",idvar = "record_id",
                        timevar = "week")

# pretty up
asd_ema_wide <- rename(asd_ema_wide[,1:14],
                       "age"="age.1","gender"="gender.1","group"="group.1",
                       "pain_intensity_mean.1" = "pain_intensity_max_mean_mean.1",
                       "pain_intensity_mean.2" = "pain_intensity_max_mean_mean.2")

# for controls, duplicate variables to match asd's week 1 and 2
control_ema_summary$pain_clean_bin_sum.2 <- control_ema_summary$pain_clean_bin_sum
control_ema_summary$physical_symp_bin_sum.2 <- control_ema_summary$physical_symp_bin_sum
control_ema_summary$pain_intensity_mean.2 <- control_ema_summary$pain_intensity_mean
control_ema_summary$day_from_start_redcap_n.2 <- control_ema_summary$day_from_start_redcap_n

control_ema_summary <- rename(control_ema_summary,
                              "pain_clean_bin_sum.1"="pain_clean_bin_sum",
                              "physical_symp_bin_sum.1"="physical_symp_bin_sum",
                              "pain_intensity_mean.1"="pain_intensity_mean",
                              "day_from_start_redcap_n.1"="day_from_start_redcap_n")

# merge
asd_ema_wide$symptom_intensity_mean_mean <- NA
control_ema_summary$meddx_sum.1 <- NA
control_ema_summary$meddx_CP.1 <- NA
all_ema_summary <- rbind(asd_ema_wide,control_ema_summary)

# save file for analyses in SPSS
write.csv(all_ema_summary,"asd to control comparison.csv")


##################################################################################################
#                                     f i g u r e s
##################################################################################################
# total pain prevalence by age and gender
ggplot(ema_summary_pp,aes(x=age,y=proportion_days_with_pain,color=gender))+
  geom_jitter(size=2)+
  theme_minimal()

# prepare for barl plot of pain prevalence per group
ema_summary_pp$days_with_pain_s <- ifelse(ema_summary_pp$pain_bin_sum ==0, "no pain",
                                      ifelse(ema_summary_pp$pain_bin_sum < 5, "1-4",
                                             ifelse(ema_summary_pp$pain_bin_sum < 9, "5-8",
                                                    ifelse(ema_summary_pp$pain_bin_sum < 13,"9-12",
                                                           ifelse(ema_summary_pp$pain_bin_sum < 15,"all the time",ema_summary_pp$sum)))))

ema_summary_pp$days_with_symp_s <- ifelse(ema_summary_pp$symptoms_bin_sum ==0, "no symptoms",
                                          ifelse(ema_summary_pp$symptoms_bin_sum < 5, "1-4",
                                                 ifelse(ema_summary_pp$symptoms_bin_sum < 9, "5-8",
                                                        ifelse(ema_summary_pp$symptoms_bin_sum < 13,"9-12",
                                                               ifelse(ema_summary_pp$symptoms_bin_sum < 15,"all the time",ema_summary_pp$sum)))))

ema_summary_pp$days_with_pain_s <- factor(ema_summary_pp$days_with_pain_s,levels=c("no pain","1-4","5-8","9-12","all the time"))
ema_summary_pp$days_with_symp_s <- factor(ema_summary_pp$days_with_symp_s,levels=c("no symptoms","1-4","5-8","9-12","all the time"))


# bar plot of pain prevalence per group
library(ggpattern)
ema_summary_groups <- ema_summary_pp %>%
  group_by(CP_nonCP_healthy,days_with_pain_s) %>%
  summarize_at(vars(record_id),list(n=length))

ggplot(ema_summary_groups,aes(x=days_with_pain_s,y=n,
                          group=as.factor(CP_nonCP_healthy),
                          pattern=as.factor(CP_nonCP_healthy),
                          pattern_angle=as.factor(CP_nonCP_healthy)))+
  geom_bar_pattern(
    stat = "identity",
    fill = "white",              # base fill (usually white for clarity),
    color = "firebrick1",
    pattern_fill = "firebrick1",      # pattern color
    pattern_size = 0.01,
    pattern_density = 0.1,
    pattern_spacing = 0.02
  )+
  scale_pattern_manual(values = c("none","circle","stripe"))+
  scale_pattern_angle_manual(values = c(0, 135, 45))+
  theme_minimal() +
  labs(x = "Days with Pain (out of 14)", y = "Number of Participants") +
  ylim(0,30)+
  theme(
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),  # top margin
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),  # right margin
    axis.text = element_text(size = 14)
  ) 

# same for the symptom
ema_summary_groups2 <- ema_summary_pp %>%
  group_by(CP_nonCP_healthy,days_with_symp_s) %>%
  summarize_at(vars(record_id),list(n=length))

ggplot(ema_summary_groups2,aes(x=days_with_symp_s,y=n,
                              group=as.factor(CP_nonCP_healthy),
                              pattern=as.factor(CP_nonCP_healthy),
                              pattern_angle=as.factor(CP_nonCP_healthy)))+
  geom_bar_pattern(
    stat = "identity",
    fill = "white",              # base fill (usually white for clarity),
    color = "steelblue",
    pattern_fill = "steelblue",      # pattern color
    pattern_size = 0.01,
    pattern_density = 0.1,
    pattern_spacing = 0.02
  )+
  scale_pattern_manual(values = c("none","circle","stripe"))+
  scale_pattern_angle_manual(values = c(0, 135, 45))+
  theme_minimal() +
  labs(x = "Days with Symptoms (out of 14)", y = "Number of Participants") +
  #ylim(0,30)+
  theme(
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),  # top margin
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),  # right margin
    axis.text = element_text(size = 14)
  ) 

#----------------------------------------------------------------------------------
#                               types of pain
#----------------------------------------------------------------------------------
# types of pain
pains <- ema_summary %>%
  group_by(record_id) %>%
  summarize_at(vars(pain_stomach_max,headache_max,
                    low_back_pain_max,limb_pain_max,
                    pain_other_max,pain_bin,n_of_pains,n_of_pains_prompt_mean),
               list(mean = ~mean(.,na.rm=TRUE),
                    max = ~max(.,na.rm=TRUE),
                    n = ~sum(!is.na(.x))))

# preparing to plot
pains <- merge(pains,ema_summary_pp[,c("record_id","CP_nonCP_healthy")],by="record_id")
pains <- arrange(pains,desc(CP_nonCP_healthy),desc(pain_bin_mean))
for (i in 1:nrow(pains)) {
  if (i==1) {pains$new_id[i] <- 1} else {
    pains$new_id[i] <- ifelse(pains$record_id[i]==pains$record_id[i-1],pains$new_id[i-1],pains$new_id[i-1]+1)  
  }
}

# normalize number of pains per day and prompt by dividing them by 5 
# (total number of pain categories prompted)
pains$percent_pains_day = pains$n_of_pains_mean/5
pains$percent_pains_prompt = pains$n_of_pains_prompt_mean_mean/5

pains_long <- rename(pains,"stomach pain"="pain_stomach_max_mean",
                     "headache"="headache_max_mean",
                     "low back pain"="low_back_pain_max_mean",
                     "limb pain"="limb_pain_max_mean",
                     "other pain"= "pain_other_max_mean",
                     "max number of areas per day" = "percent_pains_day",
                     "mean number of areas simultaneously" = "percent_pains_prompt",
                     "overall (any pain)" = "pain_bin_mean") %>%
  select("record_id","new_id","stomach pain","headache","low back pain","limb pain", "other pain",
         "mean number of areas simultaneously","max number of areas per day","overall (any pain)") %>%
  pivot_longer(cols = c("stomach pain","headache","low back pain","limb pain", "other pain",
                        "mean number of areas simultaneously","max number of areas per day","overall (any pain)"),
               names_to = "variable",
               values_to = "proportion")
pains_long$variable <- factor(pains_long$variable,
                              levels=c("overall (any pain)","stomach pain","headache","low back pain","limb pain", "other pain",
                                       "mean number of areas simultaneously","max number of areas per day"))
pains_long$proportion <- ifelse(is.na(pains_long$proportion) |
                                  is.infinite(pains_long$proportion),0,pains_long$proportion)
ggplot(pains_long, 
       aes(x = factor(new_id), y = variable, fill = proportion)) +
  geom_tile(color = "white") + #don't change this, this is just background
  geom_tile(fill = NA, color = "grey70", linewidth = 0.04) +
  scale_fill_gradient(low="white",high="firebrick") +   # This can be adjusted. This thing produces the green to yellow to orange as in my graph.
  theme_minimal() +
  #scale_x_continuous(breaks = seq(1, 7, by = 1)) + #if either of your scales is continuous, this will parse it into tiles. If all categorical, no need.
  xlab("participant number (not record_id)") +
  ylab("")+
  ggtitle("Proportion of days with:")

#----------------------------------------------------------------------------------
#                               types of symptom
#----------------------------------------------------------------------------------
# normalize number of symptoms by total n of categories = 7:
ema_summary$percent_symptoms_day <- ema_summary$n_of_symptoms/7

# prepare to plot
symptoms <- ema_summary %>%
  group_by(record_id) %>%
  summarize_at(vars(faint_dizzy_max,tachycardia_max,
                    nausea_max,weakness_max,
                    migraine_max,fatigue_max,
                    other_symptoms_max,symptoms_bin,percent_symptoms_day),
               list(mean = ~mean(.,na.rm=TRUE),
                    max = ~max(.,na.rm=TRUE),
                    n = ~sum(!is.na(.x))))

# match IDs to pain
symptoms <- merge(symptoms,pains[,c("record_id","new_id")], by="record_id")
symptoms <- arrange(symptoms,new_id)

symptoms_long <- symptoms %>%
  select("record_id","new_id","faint_dizzy_max_mean","tachycardia_max_mean",
         "nausea_max_mean","weakness_max_mean","migraine_max_mean","fatigue_max_mean",
         "other_symptoms_max_mean","percent_symptoms_day_max","symptoms_bin_mean") %>%
  pivot_longer(cols = c("faint_dizzy_max_mean","tachycardia_max_mean",
                        "nausea_max_mean","weakness_max_mean","migraine_max_mean","fatigue_max_mean",
                        "other_symptoms_max_mean","percent_symptoms_day_max","symptoms_bin_mean"),
               names_to = "variable",
               values_to = "proportion")
symptoms_long$variable <- sub("_"," ",sub("_max_mean","",symptoms_long$variable))
symptoms_long$variable <- factor(symptoms_long$variable,
                                 levels=c("symptoms bin_mean","faint dizzy","tachycardia",
                                          "nausea","weakness","migraine","fatigue",
                                          "other symptoms","percent symptoms_day_max"))
symptoms_long$proportion <- ifelse(is.na(symptoms_long$proportion) |
                                  is.infinite(symptoms_long$proportion),0,symptoms_long$proportion)

ggplot(symptoms_long, 
       aes(x = factor(new_id), y = variable, fill = proportion)) +
  geom_tile(color = "white") + #don't change this, this is just background
  geom_tile(fill = NA, color = "grey70", linewidth = 0.04) +
  scale_fill_gradient(low="white",high="steelblue") +   # This can be adjusted. This thing produces the green to yellow to orange as in my graph.
  theme_minimal() +
  #scale_x_continuous(breaks = seq(1, 7, by = 1)) + #if either of your scales is continuous, this will parse it into tiles. If all categorical, no need.
  xlab("participant number (not record_id)") +
  ylab("")+
  ggtitle("Proportion of days with:")


#-----------------------------------------------------------------------------------------------------------
#                           pain intensity and interference
#-----------------------------------------------------------------------------------------------------------
# plot intensities
ggplot(ema_summary_pp,aes(x = pain_intensity_max_mean_mean,group=as.factor(CP_nonCP_healthy),
                          pattern=as.factor(CP_nonCP_healthy),
                          pattern_angle=as.factor(CP_nonCP_healthy)))+
  #scale_fill_manual(values=c("white","grey","grey43"))+
  ylim(0,0.03)+
  geom_density_pattern(
    pattern_size = 0.1,
    pattern_density = 0.1,
    pattern_spacing = 0.05)+
  scale_pattern_manual(values = c("none","stripe","stripe"))+
  scale_pattern_angle_manual(values = c(0, 135, 45))+
  labs(x = "Mean Pain Intensity per participant")+
  theme_minimal()

# plot interference
ggplot(ema_summary_pp,aes(x = pain_interference_max_mean_mean,group=as.factor(CP_nonCP_healthy),
                          pattern=as.factor(CP_nonCP_healthy),
                          pattern_angle=as.factor(CP_nonCP_healthy)))+
  #scale_fill_manual(values=c("white","grey","grey43"))+
  ylim(0,0.03)+
  geom_density_pattern(
    pattern_size = 0.1,
    pattern_density = 0.1,
    pattern_spacing = 0.05)+
  scale_pattern_manual(values = c("none","stripe","stripe"))+
  scale_pattern_angle_manual(values = c(0, 135, 45))+
  labs(x = "Mean Pain Interference per participant")+
  theme_minimal()

#---------------------------------------------------------------------------------------------------------------------------------------------
#                               # plot control pains
#---------------------------------------------------------------------------------------------------------------------------------------------

control_ema$stomach.ache_bin <- ifelse(control_ema$stomach.ache_mean>0,1,0)
control_ema$headache_bin <- ifelse(control_ema$headache_mean>0,1,0)
control_ema$lbp_bin <- ifelse(control_ema$lbp_mean>0,1,0)
control_ema$armslegspain_bin <- ifelse(control_ema$armslegspain_mean>0,1,0)
control_ema$pain_other_bin <- ifelse(control_ema$other_mean>0,1,0)
control_ema$dizziness_bin <- ifelse(control_ema$dizziness_mean > 0,1,0)
control_ema$nausea_bin <- ifelse(control_ema$nausea_mean > 0,1,0)
control_ema$weakness_bin <- ifelse(control_ema$weakness_mean > 0,1,0)
control_ema$tachycardia_bin <- ifelse(control_ema$tachycardia_mean > 0,1,0)

control_ema$n_of_pains <- rowSums(control_ema[,c("stomach.ache_bin","headache_bin",
                                                 "lbp_bin","armslegspain_bin","pain_other_bin")])
control_ema$pain_present <- ifelse(control_ema$n_of_pains >0,1,0)
control_ema_summary2 <- control_ema %>%
  group_by(family_id) %>%
  summarize_at(vars(stomach.ache_bin,headache_bin,lbp_bin,armslegspain_bin,pain_other_bin,n_of_pains,pain_present,
                    dizziness_bin,nausea_bin,weakness_bin,tachycardia_bin),
               list(mean = ~mean(.,na.rm=TRUE),
                    max = ~max(.,na.rm=TRUE),
                    n = length))

control_ema_summary2 <- arrange(control_ema_summary2,desc(pain_present_mean))
for (i in 1:nrow(control_ema_summary2)) {
  if (i==1) {control_ema_summary2$new_id[i] <- 1} else {
    control_ema_summary2$new_id[i] <- ifelse(control_ema_summary2$family_id[i]==control_ema_summary2$family_id[i-1],
                                             control_ema_summary2$new_id[i-1],control_ema_summary2$new_id[i-1]+1)  
  }
}

pains_long_control <- rename(control_ema_summary2,"stomach pain"="stomach.ache_bin_mean",
                             "headache"="headache_bin_mean",
                             "low back pain"="lbp_bin_mean",
                             "limb pain"="armslegspain_bin_mean",
                             "other pain"= "pain_other_bin_mean",
                             "number of areas per day" = "n_of_pains_max",
                             "overall (any pain)" = "pain_present_mean",
                             "weakness" = "weakness_bin_mean",
                             "nausea" = "nausea_bin_mean",
                             "tachycardia" = "tachycardia_bin_mean",
                             "dizziness" = "dizziness_bin_mean") %>%
  select("family_id","new_id","stomach pain","headache","low back pain","limb pain", "other pain",
         "number of areas per day","overall (any pain)") %>%
  pivot_longer(cols = c("stomach pain","headache","low back pain","limb pain", "other pain",
                        "number of areas per day","overall (any pain)"),
               names_to = "variable",
               values_to = "proportion")
pains_long_control$proportion <- ifelse(pains_long_control$variable=="number of areas per day",
                                        pains_long_control$proportion/5,pains_long_control$proportion)
pains_long_control$variable <- factor(pains_long_control$variable,
                                      levels=c("overall (any pain)","stomach pain","headache","low back pain","limb pain", "other pain",
                                               "number of areas per day"))

ggplot(pains_long_control, 
       aes(x = factor(new_id), y = variable, fill = proportion)) +
  geom_tile(color = "white") +
  geom_tile(fill = NA, color = "grey70", linewidth = 0.3) +
  scale_fill_gradient(low="white",high="#995050") +
  theme_minimal() +
  xlab("participant number (not record_id)") +
  ylab("") +
  ggtitle("Proportion of days with:") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )


# plot symptoms
pains_long_control2 <- rename(control_ema_summary2,"stomach pain"="stomach.ache_bin_mean",
                              "headache"="headache_bin_mean",
                              "low back pain"="lbp_bin_mean",
                              "limb pain"="armslegspain_bin_mean",
                              "other pain"= "pain_other_bin_mean",
                              "number of areas per day" = "n_of_pains_max",
                              "overall (any pain)" = "pain_present_mean",
                              "weakness" = "weakness_bin_mean",
                              "nausea" = "nausea_bin_mean",
                              "tachycardia" = "tachycardia_bin_mean",
                              "dizziness" = "dizziness_bin_mean") %>%
  select("family_id","new_id","overall (any pain)","weakness","nausea","tachycardia","dizziness") %>%
  pivot_longer(cols = c("overall (any pain)","weakness","nausea","tachycardia","dizziness"),
               names_to = "variable",
               values_to = "proportion")
pains_long_control2$variable <- factor(pains_long_control2$variable,
                                       levels=c("dizziness","tachycardia","nausea","weakness","overall (any pain)"))

ggplot(pains_long_control2, 
       aes(x = factor(new_id), y = variable, fill = proportion)) +
  geom_tile(color = "white") +
  geom_tile(fill = NA, color = "grey70", linewidth = 0.3) +
  scale_fill_gradient(low="white",high="#105090") +
  theme_minimal() +
  xlab("participant number (not record_id)") +
  ylab("") +
  ggtitle("Proportion of days with:") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )




