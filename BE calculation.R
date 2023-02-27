library(data.table)

load("multicheck_model.calib.Rdata")
parms.samp <- multicheck$parms.samp
posteriors.samp <- parms.samp[,2:11]

datafolder <- "InputData"
set.seed(3.14159)

### POD from BBMD
bmd_m.df <- fread(file.path(datafolder,"don-male-bmds.csv"))
bmd_f.df <- fread(file.path(datafolder,"don-female-bmds.csv"))
bmd_m <- bmd_m.df$model_average
bmd_f <- bmd_f.df$model_average
bmd <- 1000*c(bmd_m,bmd_f)[sample.int(30000,10000)]
cat("bmd",quantile(bmd,prob=c(0.5,0.05,0.95)),"\n")

#Cavg_a
#Cavg/dose = (AUC/dose)/24 hr
#Cavg_a = bmd * (Cavg/dose)
#Cavg_dose_a.GM <- 0.8/24
#Cavg_dose_a.GSD <- 3^(1/1.96)
Cavg_dose_a <- rlnorm(10000,meanlog = log(0.033), sdlog=log(1.75))
cat("Cavg_dose_a",quantile(Cavg_dose_a,prob=c(0.5,0.05,0.95)),"\n")
AF_interTKblood <- 1/Cavg_dose_a
cat("AF_interTKblood", quantile(AF_interTKblood,prob=c(0.5,0.05,0.95)),"\n")

### Inter-species TK
## AUC per unit dose calculation
# AUC = Dose * Fgutabs /(Vd * ktot)
#convert to human TK
#U/Cavg
M_lnU_Cavg <- posteriors.samp$M_lnktot.1. + log(1.24*24) # Vd fixed; k unit: 1/hr->per day
SD_lnU_Cavg <- posteriors.samp$SD_lnktot.1.
GM_U_Cavg <- exp(M_lnU_Cavg)
U_Cavg_h.GM <- rep(GM_U_Cavg,20)[sample.int(10000)]
cat("U_Cavg_h.GM",quantile(U_Cavg_h.GM,prob=c(0.5,0.05,0.95)),"\n")
AF_interTKurine <- 1/(U_Cavg_h.GM*Cavg_dose_a)
cat("AF_interTKurine",quantile(AF_interTKurine,prob=c(0.5,0.05,0.95)),"\n")

### Inter-species TD
# Default TK/TD has GM=1, GSD=1.95
# Assume TK and TD are equal and independent
Inter_sigmaTD <- log(1.95)/sqrt(2)
AF_interTD <- rlnorm(10000,meanlog = 0, sdlog=Inter_sigmaTD)
cat("AF_interTD",quantile(AF_interTD,prob=c(0.5,0.05,0.95)),"\n")

#Intra-species TK
GSD_U_Cavg <- exp(SD_lnU_Cavg)
TK_var.Cavg_h <- rep(GSD_U_Cavg,20)[sample.int(10000)]
cat("TK_var.Cavg_h", quantile(TK_var.Cavg_h, prob=c(0.5, 0.05, 0.95)),"\n")

### Intra-species TD
TD_var.df <- fread(file.path(datafolder,"DON-EC10_pop_samples.csv"))
TD_var.GSD <- TD_var.df$EC10.GSD[sample.int(4000,size = 10000,replace=TRUE)]
AF_intraTD <- TD_var.GSD^qnorm(0.99)
cat("TD_var.GSD",quantile(TD_var.GSD,prob=c(0.5,0.05,0.95)),"\n")
cat("AF_intraTD",quantile(AF_intraTD,prob=c(0.5,0.05,0.95)),"\n")

### Intra-species TK+TD
Intra.GSD <- exp(sqrt(log(TK_var.Cavg_h)^2+log(TD_var.GSD)^2))
AF_intraTKTD.be <- Intra.GSD^qnorm(0.99)
cat("AF_intraTKTD.be",quantile(AF_intraTKTD.be,prob=c(0.5,0.05,0.95)),"\n")

#### Urine BE
BE.urine <- (bmd*Cavg_dose_a*U_Cavg_h.GM) / (AF_interTD*AF_intraTKTD.be)
cat("BE.urine",quantile(BE.urine,prob=c(0.5,0.05,0.95)),"\n")

BE.urine.df <- data.frame(bmd, Cavg_dose_a, GM_U_Cavg, U_Cavg_h.GM, AF_interTKurine,
                    AF_interTD,GSD_U_Cavg, TK_var.Cavg_h,
                    TD_var.GSD, Intra.GSD, AF_intraTKTD.be, BE.urine)
fwrite(BE.urine.df,"Urine BE.samples.csv")

BE.urine.quan <- sapply(BE.urine.df, quantile, probs=c(0.5,0.05,0.95))
write.csv(BE.urine.quan, file="BE.urine.quantiles.csv")

#### Blood BE
BE.blood <- (bmd*Cavg_dose_a) / (AF_interTD*AF_intraTD)
cat("BE.blood",quantile(BE.blood,prob=c(0.5,0.05,0.95)),"\n")

BE.blood.df <- data.frame(bmd, Cavg_dose_a, AF_interTKblood, AF_interTD,
                          TD_var.GSD, AF_intraTD, BE.blood)
fwrite(BE.blood.df,"Blood BE.samples.csv")

BE.blood.quan <- sapply(BE.blood.df, quantile, probs=c(0.5,0.05,0.95))
write.csv(BE.blood.quan, file="BE.blood.quantiles.csv")
