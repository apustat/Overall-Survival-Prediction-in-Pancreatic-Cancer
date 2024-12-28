
######################################
library(glmnet)
library(stabs)
library("rpart")
library("survival")
library(dplyr)


dat<-read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Pancreatic cancer/logpc_new.csv", header=TRUE)
dat2<-as.data.frame(na.omit(dat)) #129 observations

#many_na=(colMeans(is.na(dat2)))*100
#y=Surv(dat2$time_to_mortality_exit/365.25,dat2$is_dead)
y=Surv(dat2$oasyrs,dat2$is_dead)
biodat<-dat2[,c(45,47:77)]

biodat=as.matrix(biodat)

set.seed(3355660)
cvfit.bio=cv.glmnet(biodat, y, family = "cox",alpha=1,standardize =TRUE)
(bestlam.bio=cvfit.bio$lambda.min)
plot(cvfit.bio)
fit.bio = glmnet(biodat,y , family = "cox",alpha=1,standardize =TRUE) #lasso regression
plot(fit.bio)
coef(fit.bio,s =bestlam.bio)

#stab.rconcave=stabsel(x=biodat, y=y, fitfun=glmnet.lasso, cutoff=0.75, PFER = 1, sampling.type = "SS", assumption="r-concave")
##################################################################################
#dat2$panc_stage<- recode(dat2$panc_pcs_summary_stage , "1"=1, "2"=2, "7"=7, "9"=9, "N"=10)

clindat=dat2[, c("Age","GENDER", "bmi_curr","panc_pcs_summary_stage","cig_stat","diabetes_f","panc_fh","arthrit_f","bronchit_f","colitis_f",
                 "divertic_f","emphys_f","gallblad_f","hearta_f","hepatit_f","hyperten_f","osteopor_f","polyps_f","stroke_f","race7")]
clindat$Age=as.numeric(clindat$Age)
clindat$GENDER=as.factor(ifelse(clindat$GENDER=="F",0,1))
clindat$panc_pcs_summary_stage=as.factor(clindat$panc_pcs_summary_stage)

clindat$cig_stat=as.factor(clindat$cig_stat)
clindat$diabetes_f=as.factor(clindat$diabetes_f)
clindat$panc_fh=as.factor(clindat$panc_fh)
clindat$arthrit_f=as.factor(clindat$arthrit_f)
clindat$bronchit_f=as.factor(clindat$bronchit_f)
#clindat$cirrhos_f=as.factor(clindat$cirrhos_f) ##Has one category only
clindat$colitis_f=as.factor(clindat$colitis_f)
#clindat$crohn_f=as.factor(clindat$crohn_f) ##Has one category only
clindat$divertic_f=as.factor(clindat$divertic_f)
clindat$emphys_f=as.factor(clindat$emphys_f)
clindat$gallblad_f=as.factor(clindat$gallblad_f)
#clindat$gardner_f=as.factor(clindat$gardner_f) ##Has one category only
clindat$hearta_f=as.factor(clindat$hearta_f)
clindat$hepatit_f=as.factor(clindat$hepatit_f)
clindat$hyperten_f=as.factor(clindat$hyperten_f)
clindat$osteopor_f=as.factor(clindat$osteopor_f)
#clindat$polypos_f=as.factor(clindat$polypos_f) ##Has one category only after replacing the NA by 0
clindat$polyps_f=as.factor(clindat$polyps_f)
clindat$stroke_f=as.factor(clindat$stroke_f)
clindat$race7=as.factor(clindat$race7)

summary(clindat)

clindat1=as.matrix(clindat)
clindat2=matrix(as.numeric(clindat1), ncol=ncol(clindat1))

set.seed(3355660)
cvfit.clin=cv.glmnet(clindat2, y , family = "cox",alpha=1,standardize =TRUE)
(bestlam.clin=cvfit.clin$lambda.min)
plot(cvfit.clin)
fit.clin = glmnet(clindat, y , family = "cox",alpha=1,standardize =TRUE) #lasso regression
plot(fit.clin)
coef(fit.clin,s =bestlam.clin)
############################################

bothdat=cbind(biodat, clindat)
bothdat1=as.matrix(bothdat)
bothdat2=matrix(as.numeric(bothdat1), ncol=ncol(bothdat1))

set.seed(3355660)
cvfit.both=cv.glmnet(bothdat2,y , family = "cox",alpha=1,standardize =TRUE)
(bestlam.both=cvfit.both$lambda.min)
plot(cvfit.both)
fit.both = glmnet(bothdat,y , family = "cox",alpha=1,standardize =TRUE) #lasso regression
plot(fit.both)
coef(fit.both,s =bestlam.both)

############################
#fulldata=cbind(as.data.frame(biodat), as.data.frame(clindat), time_to_mortality_exit=dat2$time_to_mortality_exit, is_dead=dat2$is_dead)
#write.csv(dat2, "C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Pancreatic cancer/fulldat.csv", row.names=FALSE)

###################################################################################################################################
dat2$oasyrs_cat=ifelse(dat2$oasyrs<1, 0,1)

y_cat=dat2$oasyrs_cat
biodat<-dat2[,c(45,47:77)]

biodat=as.matrix(biodat)

set.seed(3)
cvfit.bio=cv.glmnet(biodat, y_cat, family = "binomial",alpha=1, standardize =TRUE)
(bestlam.bio=cvfit.bio$lambda.min)
fit.bio = glmnet(biodat,y_cat, family = "binomial",alpha=1,standardize =TRUE) #lasso regression
coef(fit.bio,s =bestlam.bio)



###################
clindat=dat2[, c("Age","GENDER", "bmi_curr","panc_pcs_summary_stage","cig_stat","diabetes_f","panc_fh","arthrit_f","bronchit_f","colitis_f",
                 "divertic_f","emphys_f","gallblad_f","hearta_f","hepatit_f","hyperten_f","osteopor_f","polyps_f","stroke_f","race7")]
clindat$Age=as.numeric(clindat$Age)
clindat$GENDER=as.factor(ifelse(clindat$GENDER=="F",0,1))
clindat$panc_pcs_summary_stage=as.factor(clindat$panc_pcs_summary_stage)

clindat$cig_stat=as.factor(clindat$cig_stat)
clindat$diabetes_f=as.factor(clindat$diabetes_f)
clindat$panc_fh=as.factor(clindat$panc_fh)
clindat$arthrit_f=as.factor(clindat$arthrit_f)
clindat$bronchit_f=as.factor(clindat$bronchit_f)
#clindat$cirrhos_f=as.factor(clindat$cirrhos_f) ##Has one category only
clindat$colitis_f=as.factor(clindat$colitis_f)
#clindat$crohn_f=as.factor(clindat$crohn_f) ##Has one category only
clindat$divertic_f=as.factor(clindat$divertic_f)
clindat$emphys_f=as.factor(clindat$emphys_f)
clindat$gallblad_f=as.factor(clindat$gallblad_f)
#clindat$gardner_f=as.factor(clindat$gardner_f) ##Has one category only
clindat$hearta_f=as.factor(clindat$hearta_f)
clindat$hepatit_f=as.factor(clindat$hepatit_f)
clindat$hyperten_f=as.factor(clindat$hyperten_f)
clindat$osteopor_f=as.factor(clindat$osteopor_f)
#clindat$polypos_f=as.factor(clindat$polypos_f) ##Has one category only after replacing the NA by 0
clindat$polyps_f=as.factor(clindat$polyps_f)
clindat$stroke_f=as.factor(clindat$stroke_f)
clindat$race7=as.factor(clindat$race7)

summary(clindat)

clindat1=as.matrix(clindat)
clindat2=matrix(as.numeric(clindat1), ncol=ncol(clindat1))

set.seed(3)
cvfit.clin=cv.glmnet(clindat2, y_cat , family = "binomial",alpha=1,standardize =TRUE)
(bestlam.clin=cvfit.clin$lambda.min)
fit.clin = glmnet(clindat, y_cat, family = "binomial",alpha=1,standardize =TRUE) #lasso regression
coef(fit.clin,s =bestlam.clin)


############################################
bothdat=cbind(biodat, clindat)
bothdat1=as.matrix(bothdat)
bothdat2=matrix(as.numeric(bothdat1), ncol=ncol(bothdat1))

set.seed(300)
cvfit.both=cv.glmnet(bothdat2,y_cat, family = "binomial",alpha=1,standardize =TRUE)
(bestlam.both=cvfit.both$lambda.min)
fit.both = glmnet(bothdat,y_cat, family = "binomial",alpha=1,standardize =TRUE) #lasso regression
coef(fit.both,s =bestlam.both)


###################Heatmap#############################
biomarkers<-dat2[,47:77]

library(reshape2)
library(ggplot2)
# creating correlation matrix
corr_mat <- round(cor(biomarkers),2)

# reduce the size of correlation matrix
melted_corr_mat <- melt(corr_mat)


ggplot(data = melted_corr_mat, aes(x=Var1, y=Var2, fill=value)) + 
geom_tile()+
geom_text(aes(Var2, Var1, label = value), color = "black", size = 2)+
scale_fill_gradient2(low = "royalblue",mid="snow3",high = "red", guide = "colorbar",breaks = c(-1.0,-0.5,0.0,0.5,1.0),limits = c(-1, 1))+
  theme(axis.text.x = element_text(angle = 50,vjust = 1,hjust = 1))+
  xlab("Biomarkers (log scale)") + ylab("Biomarkers (log scale)")


###########Stability selection####################################
library(stabs)
dat<-read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Pancreatic cancer/fulldat1.csv", header=TRUE)

biomarkers<-dat[,c(24,47:77)]

stab.rconcave.bio=stabsel(x=biomarkers, y=dat$oasyrs_cat,fitfun=glmnet.lasso, cutoff=0.55, PFER = 3, sampling.type = "SS", assumption="r-concave")
stab.rconcave.bio

clindat=dat[, c("Age","GENDER", "bmi_curr","panc_pcs_summary_stage","cig_stat","diabetes_f","panc_fh","arthrit_f","bronchit_f","colitis_f",
                 "divertic_f","emphys_f","gallblad_f","hearta_f","hepatit_f","hyperten_f","osteopor_f","polyps_f","stroke_f","race")]
clindat$Age=as.numeric(clindat$Age)
clindat$GENDER=as.factor(ifelse(clindat$GENDER=="F",0,1))
clindat$panc_pcs_summary_stage=as.factor(clindat$panc_pcs_summary_stage)

clindat$cig_stat=as.factor(clindat$cig_stat)
clindat$diabetes_f=as.factor(clindat$diabetes_f)
clindat$panc_fh=as.factor(clindat$panc_fh)
clindat$arthrit_f=as.factor(clindat$arthrit_f)
clindat$bronchit_f=as.factor(clindat$bronchit_f)
#clindat$cirrhos_f=as.factor(clindat$cirrhos_f) ##Has one category only
clindat$colitis_f=as.factor(clindat$colitis_f)
#clindat$crohn_f=as.factor(clindat$crohn_f) ##Has one category only
clindat$divertic_f=as.factor(clindat$divertic_f)
clindat$emphys_f=as.factor(clindat$emphys_f)
clindat$gallblad_f=as.factor(clindat$gallblad_f)
#clindat$gardner_f=as.factor(clindat$gardner_f) ##Has one category only
clindat$hearta_f=as.factor(clindat$hearta_f)
clindat$hepatit_f=as.factor(clindat$hepatit_f)
clindat$hyperten_f=as.factor(clindat$hyperten_f)
clindat$osteopor_f=as.factor(clindat$osteopor_f)
#clindat$polypos_f=as.factor(clindat$polypos_f) ##Has one category only after replacing the NA by 0
clindat$polyps_f=as.factor(clindat$polyps_f)
clindat$stroke_f=as.factor(clindat$stroke_f)
clindat$race=as.factor(clindat$race)

bothdat=cbind(biomarkers, clindat)

stab.rconcave.both=stabsel(x=bothdat, y=dat$oasyrs_cat,fitfun=glmnet.lasso, cutoff=0.55, PFER = 3, sampling.type = "SS", assumption="r-concave")
stab.rconcave.both


###########################################################
library(haven)
data=read_sas("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Pancreatic cancer/panc_biomarker_new.sas7bdat")
data2=as.data.frame(na.omit(data))
