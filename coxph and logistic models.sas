proc import file="/home/u58303998/P_cancer/fulldat.csv"
    out=fulldat dbms=csv replace; /*No missing value in fulldat; */
    getnames=yes;
run;

proc means data=fulldat median;
var oasyrs;
run;

data fulldat; set fulldat;
If race7=1 then race="White,Non-Hispanic";
Else race="Others";
run; 

proc freq data=fulldat;
table race; run;

ods rtf file="/home/u58303998/P_cancer/model_summary_final.rtf" style=Statistical;
************************Cox PH model*******************;
*Model 1: Cox PH model with biomarkers only;
Proc phreg data=fulldat;
MODEL oasyrs*is_dead(0)= time_to_dx_mon logHE4_CD_ logCA_72_4_CD_ logCytokeratin_19_CD_ logCA_125_Fujio_CD_ logCEA_Fujio_CD_ logErbB2_CD_
            logErbB2_CD_ logEGFR_CD_ logCA_19_9_CD_ logOPG_M_ logOC_M_ logPTH_M_ logFSH__M_ logBDNF__M_ logGH__M_ logtPAI_1__M_ 
            logNCAM__M_ logMIF_M_ logNSE__M_ logOsteonectin__M_ logPeriostin__M_ logTRAP5__M_  logYKL40__M_ logsIL_1R1__M_ logsIL_4R__M_ 
            logsVEGFR1__M_ logsVEGFR2__M_ logsVEGFR3__M_  logHIF_1a__CD_ logALCAM__CD_ logCEACAM_1__CD_ logCEACAM_6__CD_/ 
            rl selection=backward slstay=0.05;
run;

*time_to_dx_month is not selected by backward selection,  so included forcefully;
Proc phreg data=fulldat;
MODEL  oasyrs*is_dead(0)= time_to_dx_mon logCA_19_9_CD_ logMIF_M_ logOsteonectin__M_ logCEACAM_6__CD_ / rl;
run;

*model 2: Cox PH model with clinical variables only; 
Proc phreg data=fulldat;
class GENDER panc_pcs_summary_stage (param=ref ref="1") cig_stat diabetes_f(param=ref ref="0") panc_fh arthrit_f bronchit_f colitis_f divertic_f
                emphys_f gallblad_f hearta_f hepatit_f hyperten_f osteopor_f polyps_f stroke_f race;
MODEL oasyrs*is_dead(0)= Age GENDER bmi_curr panc_pcs_summary_stage cig_stat diabetes_f panc_fh 
                 hyperten_f osteopor_f polyps_f stroke_f race/ rl selection=backward slstay=0.05;
run;

*Model 3: Cox PH model with selected variables from model 1 and 2;
proc phreg data=fulldat;
class panc_pcs_summary_stage (param=ref ref="1") diabetes_f(param=ref ref="0");
MODEL oasyrs*is_dead(0)= time_to_dx_mon logCA_19_9_CD_ logMIF_M_ logOsteonectin__M_ logCEACAM_6__CD_
	Age panc_pcs_summary_stage diabetes_f/ rl;
run;

******************Logistic regression************************;
*************************************************************;
data fulldat;
set fulldat;
if oasyrs < 1 then  oasyrs_cat="0";
else oasyrs_cat=1;
run;

proc contents data=fulldat; run;

*************Chisquare**********************;
/* PROC FREQ data=fulldat; */
/*     TABLE oasyrs_cat*GENDER oasyrs_cat*panc_pcs_summary_stage oasyrs_cat*cig_stat oasyrs_cat*diabetes_f */
/*     oasyrs_cat*panc_fh oasyrs_cat*arthrit_f oasyrs_cat*bronchit_f oasyrs_cat*colitis_f oasyrs_cat*divertic_f */
/*     oasyrs_cat*emphys_f oasyrs_cat*gallblad_f oasyrs_cat*hearta_f oasyrs_cat*hepatit_f oasyrs_cat*hyperten_f */
/*     oasyrs_cat*osteopor_f oasyrs_cat*polyps_f oasyrs_cat*stroke_f oasyrs_cat*race/ CHISQ; */
/* RUN; */

*chisquare tests show warning becasue of the small cell counts;

*************Logistic model*****************;

*Model 4: Logistic model with Biomarkers only: Used Stability selection with LASSO, then logistic model with selected biomarkers;

Proc logistic data=fulldat descending;
class oasyrs_cat;
model oasyrs_cat = logMIF_M_ logCEACAM_6__CD_ /rl;
run;

*Model 5: Logistic model with clinical variables only;

Proc logistic data=fulldat descending;
class oasyrs_cat GENDER panc_pcs_summary_stage (param=ref ref="1") cig_stat diabetes_f (param=ref ref="0")  panc_fh arthrit_f bronchit_f colitis_f divertic_f
                emphys_f gallblad_f hearta_f hepatit_f hyperten_f osteopor_f polyps_f stroke_f race;
MODEL oasyrs_cat= Age panc_pcs_summary_stage diabetes_f /  firth;
run;



*Model 6: Logistic model biomarkers and clinical variables: Used Stability selection with LASSO, then logistic model with selected variables;

Proc logistic data=fulldat descending;
class oasyrs_cat panc_pcs_summary_stage (param=ref ref="1") diabetes_f (param=ref ref="0") polyps_f(param=ref ref="0");
MODEL oasyrs_cat= logMIF_M_ panc_pcs_summary_stage diabetes_f polyps_f/rl firth;
run;

******************************************************************;
*Kaplan-Meier for significant biomarkers;
proc univariate data=fulldat;
var logHE4_CD_ logCA_19_9_CD_ logMIF_M_ logOsteonectin__M_ logCEACAM_6__CD_ logOPG_M_ logGH__M_ logTRAP5__M_;
run;

data fulldat1;
set fulldat;
if logHE4_CD_ < 7.31288 then  HE4_cat="1";
else if logHE4_CD_ >= 7.31288 and logHE4_CD_< 8.38802 then  HE4_cat="2";
else if logHE4_CD_ >= 8.38802 and logHE4_CD_< 9.58496 then  HE4_cat="3";
else HE4_cat="4";

if logCA_19_9_CD_ < 6.83289 then  CA19_cat="1";
else if logCA_19_9_CD_ >= 6.83289 and logCA_19_9_CD_< 7.57743 then  CA19_cat="2";
else if logCA_19_9_CD_ >= 7.57743 and logCA_19_9_CD_< 9.30834 then  CA19_cat="3";
else CA19_cat="4";

if logMIF_M_ < 3.62059 then  MIF_cat="1";
else if logMIF_M_ >= 3.62059 and logMIF_M_< 4.06178 then  MIF_cat="2";
else if logMIF_M_ >= 4.06178 and logMIF_M_< 4.74416 then  MIF_cat="3";
else MIF_cat="4";

if logOsteonectin__M_ < 4.59096 then  Osteonectin_cat="1";
else if logOsteonectin__M_ >= 4.59096 and logOsteonectin__M_< 4.80735 then  Osteonectin_cat="2";
else if logOsteonectin__M_ >= 4.80735 and logOsteonectin__M_< 5.05311 then  Osteonectin_cat="3";
else Osteonectin_cat="4";

if logCEACAM_6__CD_ < 5.08746 then  CEACAM_cat="1";
else if logCEACAM_6__CD_ >= 5.08746 and logCEACAM_6__CD_< 5.14975 then  CEACAM_cat="2";
else if logCEACAM_6__CD_ >= 5.14975 and logCEACAM_6__CD_< 5.23266 then  CEACAM_cat="3";
else CEACAM_cat="4";

if logOPG_M_ < 7.33985 then  OPG_cat="1";
else if logOPG_M_ >= 7.33985 and logOPG_M_< 7.77479 then  OPG_cat="2";
else if logOPG_M_ >= 7.77479 and logOPG_M_< 8.09803 then  OPG_cat="3";
else OPG_cat="4";

if logGH__M_ < 7.90087 then  GH_cat="1";
else if logGH__M_ >= 7.90087 and logGH__M_< 9.08481 then  GH_cat="2";
else if logGH__M_ >= 9.08481 and logGH__M_< 10.37178 then  GH_cat="3";
else GH_cat="4";

if logTRAP5__M_ < 9.45533 then  TRAP5_cat="1";
else if logTRAP5__M_ >= 9.45533 and logTRAP5__M_< 9.89178 then  TRAP5_cat="2";
else if logTRAP5__M_ >= 9.89178 and logTRAP5__M_< 10.16616 then  TRAP5_cat="3";
else TRAP5_cat="4";
run;

PROC EXPORT DATA=fulldat1 OUTFILE="/home/u58303998/P_cancer/fulldat1.csv"; run;

proc freq data=fulldat1;
table HE4_cat CA19_cat MIF_cat Osteonectin_cat CEACAM_cat OPG_cat GH_cat TRAP5_cat;
run;

*******************************************************;
proc lifetest data=fulldat1 plots=Survival (TEST);
time oasyrs*is_dead(0);
strata HE4_cat;
run;

proc lifetest data=fulldat1 plots=Survival (TEST);
time oasyrs*is_dead(0);
strata CA19_cat;
run;

proc lifetest data=fulldat1 plots=Survival (TEST);
time oasyrs*is_dead(0);
strata MIF_cat;
run;

proc lifetest data=fulldat1 plots=Survival (TEST);
time oasyrs*is_dead(0);
strata Osteonectin_cat;
run;

proc lifetest data=fulldat1 plots=Survival (TEST);
time oasyrs*is_dead(0);
strata CEACAM_cat;
run;

proc lifetest data=fulldat1 plots=Survival (TEST);
time oasyrs*is_dead(0);
strata OPG_cat;
run;

proc lifetest data=fulldat1 plots=Survival (TEST);
time oasyrs*is_dead(0);
strata GH_cat;
run;

proc lifetest data=fulldat1 plots=Survival (TEST);
time oasyrs*is_dead(0);
strata TRAP5_cat;
run;

ods rtf close;