library(data.table)
library(ggplot2)
set.seed(3.14159)
datafolder <- "InputData"
### POD from BBMD
bmd_m.df <- fread(file.path(datafolder,"iverson-f-et-al-1995-male-bmds.csv"))
bmd_f.df <- fread(file.path(datafolder,"iverson-f-et-al-1995-female-bmds.csv"))
bmd_m <- bmd_m.df$model_average
bmd_f <- bmd_f.df$model_average
bmd <- 1000*c(bmd_m,bmd_f)[sample.int(10000)]
cat("bmd",quantile(bmd,prob=c(0.5,0.05,0.95)),"\n")

### Inter-species TK
# Mouse TK from literature search [NEED TO CHECK THIS]
# 1.64 mg-hr/L per mg/kg from Clark et al. 2015 in Adult mice
# 4.51 mg-hr/L per mg/kg from Clark et al. 2015 in Aged mice
# Assume GM = 1.64, upper 95% Conf bound = 4.51
# 90% CI is 0.60 - 4.51, and includes other values from Pestka et al. 2017
# of 0.66 - 0.91 mg-hr/L per mg/kg
AUC_dose_m <- rlnorm(10000,meanlog = log(1.64), sdlog=log((4.51/1.64)^(1/qnorm(0.95))))
cat("AUC_dose_m",quantile(AUC_dose_m,prob=c(0.5,0.05,0.95)),"\n")

# Mouse TK from Faeste et al.
# CL in L/hr = 0.4082 * BW^0.8737
# AUC/dose = BW/CL = Fabs*BW^0.1263 / 0.4082 (kg*hr/L)
# BW in Iverson study is (0.04385,0.04154) for males/females
# = Fabs*1.64 

# Human from TK model
# ktot = 0.4 hr^-1
# Fabs = 0.6
# Vd = 1.24 L/kg
# AUC/dose = Fabs/(Vd * ktot)
AUC_dose_h.df <- fread(file.path(datafolder,"AUC_dose_samples.csv"))
AUC_dose_h.GM <- rep(AUC_dose_h.df$GM_AUC_dose,20)[sample.int(10000)]
cat("AUC_dose_h.GM",quantile(AUC_dose_h.GM,prob=c(0.5,0.05,0.95)),"\n")

# Human TK from Faeste et al.
# CL in L/hr = 0.4082 * BW^0.8737
# AUC/dose = Fabs*BW/CL = Fabs*BW^0.1263 / 0.4082 (kg*hr/L)
# BW of 70 kg 
# = Fabs*4.19

# Maybe should just do ratio of clearances -- assume Fabs is the same for diet
# CL_m = 0.4082*0.04^(0.8737-1) = 0.613 L/(hr kg)
# CL_h = 1.24 * 0.4 = 0.5 L/(hr kg)
# CL_m/CL_h = 1.23 [instead of 0.86]
# From allometric scaling from Faeste (70/0.04)^0.1263 = 2.57
# From WHO allometric scaling (50/0.04)^0.3 = 8.5


# Combined Mouse -> Human AF
AF_interTK <- AUC_dose_h.GM/AUC_dose_m
cat("AF_interTK",quantile(AF_interTK,prob=c(0.5,0.05,0.95)),"\n")

### Inter-species TD
# Default TK/TD has GM=1, GSD=1.95
# Assume TK and TD are equal and independent
Inter_sigmaTD <- log(1.95)/sqrt(2)
AF_interTD <- rlnorm(10000,meanlog = 0, sdlog=Inter_sigmaTD)
cat("AF_interTD",quantile(AF_interTD,prob=c(0.5,0.05,0.95)),"\n")

### Intra-species TK
TK_var.GSD <- rep(AUC_dose_h.df$GSD_AUC_dose,20)[sample.int(10000)]
cat("TK_var.GSD",quantile(TK_var.GSD,prob=c(0.5,0.05,0.95)),"\n")

### Intra-species TD
TD_var.df <- fread(file.path(datafolder,"DON-EC10_pop_samples.csv"))
TD_var.GSD <- TD_var.df$EC10.GSD[sample.int(4000,size = 10000,replace=TRUE)]
cat("TD_var.GSD",quantile(TK_var.GSD,prob=c(0.5,0.05,0.95)),"\n")

### Intra-species TK+TD
Intra_var.GSD <- exp(sqrt(log(TK_var.GSD)^2+log(TD_var.GSD)^2))
AF_intraTKTD <- Intra_var.GSD^qnorm(0.99)
cat("AF_intraTKTD",quantile(AF_intraTKTD,prob=c(0.5,0.05,0.95)),"\n")

### HDMI
HDMI <- bmd / (AF_interTK*AF_interTD*AF_intraTKTD)
cat("HDMI",quantile(HDMI,prob=c(0.5,0.05,0.95)),"\n")

HDMI.df <- data.frame(bmd,AUC_dose_h.GM,AUC_dose_m,AF_interTK,AF_interTD,
                      TK_var.GSD,TD_var.GSD,Intra_var.GSD,AF_intraTKTD,HDMI)
fwrite(HDMI.df,"HDMI.samples.csv")
