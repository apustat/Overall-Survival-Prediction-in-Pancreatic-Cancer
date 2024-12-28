dat<-read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/Pancreatic cancer/fulldat.csv", header=TRUE)
dat2<-as.data.frame(na.omit(dat))
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
