METHODS 

Data source: This study used the data collected as part of the Prostate, Lung, Colorectal, and Ovarian (PLCO) trial, a randomized multicenter trial conducted in the United States, 
to see if specific screening exams reduced mortality from prostate, lung, colorectal, and ovarian cancer (Andriole et al., 2012). All subjects in this study were 
enrolled with written consent as a requirement for participation in the PLCO trial. Recruitment and randomization for the study began in November 1993 and ended in 
July 2001. Subjects with a prior personal history of PLCO cancers, ongoing cancer treatment (excluding basal-cell and squamous-cell skin cancer), participation in 
another cancer screening or cancer primary prevention trial, and a recent screening test for prostate or colorectal cancer were excluded from the study. At the time 
of baseline entry, the cohort consisted of approximately 155,000 men and women aged 55 to 74 years. At the start of the study, participants completed a baseline 
questionnaire that included demographic, personal, and medical information, including diabetes status. Cases of pancreatic cancer were discovered through self-report 
by yearly mail-in surveys, state cancer registries, death certificates, physician referrals, and reports from next of relatives for dead individuals. PLCO staff obtained
and verified all medical and pathologic records related to pancreatic cancer diagnosis and supporting documentation. During the study's 13-year follow-up period, 767 
people developed pancreatic cancer. Pancreatic cancer cases with information on gender, age, race, smoking status, co-morbid conditions, family history of pancreatic 
cancer, and cancer characteristics at diagnosis (morphology, stage) were extracted from PLCO datasets for this study. 

Biomarker measurement: Biomarkers were measured on pre-diagnostic samples collected at regular intervals through the PLCO screening trial. The time from biomarker sample
collection to diagnosis was calculated in months. Thirty-one targeted biomarkers were measured on the Luminx platform.

Statistical Methods: Subject demographic and clinical variables are summarized as means and standard deviations (SD) for continuous variables and counts and percentages 
for categorical variables. Biomarker data were evaluated for the skewness and kurtosis by generating histograms. Biomarkers were log2-transformed to normalize the distributions. 
Pearson correlation was calculated for the biomarkers to assess the extent of multicollinearity and was displayed using a heatmap. 
Overall survival is defined as the time from diagnosis to last follow-up. The Kaplan–Meier method was used to estimate overall survival and the log-rank test was used to 
compare survival distributions between groups. Cox regression was used to assess the relationship between continuous biomarkers with overall survival, adjusting for multiple 
comparisons with the Benjamini-Hochberg false discovery rate (FDR) method. Three multivariate Cox models were run to determine independent predictors of overall survival among the
clinical variables only, biomarkers only and then the combined clinical and biomarkers. Lasso was used to select the important variables and then models were refit outside of the Lasso 
framework in order to estimate hazard ratios and confidence intervals. Uno’s c-statistic with 95% confidence interval (CI) is used to assess the models predictive performance, which is 
similar to the area under the receiver operator curve (AUC). Additionally, time-dependent receiver operating characteristic (ROC) curves were generated to evaluate the discriminative 
ability of the survival models at different time points. The time-dependent AUC was calculated to assess the predictive accuracy of the models over time, accounting for censoring in 
survival data. These curves provide insight into how well the selected biomarkers and clinical variables predict survival at various time horizons. Statistical analysis was performed
using R software and SAS Software version 9.4 (SAS Institute Inc., Cary, NC). 
