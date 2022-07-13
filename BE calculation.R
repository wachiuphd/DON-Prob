library(data.table)

load("multicheck_model.calib.Rdata")
parms.samp <- multicheck$parms.samp
posteriors.samp <- parms.samp[,2:11]

datafolder <- "InputData"
### POD from BBMD
bmd_m.df <- fread(file.path(datafolder,"iverson-f-et-al-1995-male-bmds.csv"))
bmd_f.df <- fread(file.path(datafolder,"iverson-f-et-al-1995-female-bmds.csv"))
bmd_m <- bmd_m.df$model_average
bmd_f <- bmd_f.df$model_average
bmd <- 1000*c(bmd_m,bmd_f)[sample.int(10000)]
cat("bmd",quantile(bmd,prob=c(0.5,0.05,0.95)),"\n")

### Inter-species TK
## AUC per unit dose calculation
# AUC = Dose * Fgutabs /(Vd * ktot)
#convert to human TK
M_lnAUC <- posteriors.samp$M_lnktot.1. + log(1.24) # Vd fixed
SD_lnAUC <- posteriors.samp$SD_lnktot.1.
GM_AUC <- exp(M_lnAUC)
AUC_h.GM <- rep(GM_AUC,20)[sample.int(10000)]
cat("AUC_h.GM",quantile(AUC_h.GM,prob=c(0.5,0.05,0.95)),"\n")

### Inter-species TD
# Default TK/TD has GM=1, GSD=1.95
# Assume TK and TD are equal and independent
Inter_sigmaTD <- log(1.95)/sqrt(2)
AF_interTD <- rlnorm(10000,meanlog = 0, sdlog=Inter_sigmaTD)
cat("AF_interTD",quantile(AF_interTD,prob=c(0.5,0.05,0.95)),"\n")

#Intra-species TK
TK_var.u.GSD <- exp(SD_lnAUC)
cat("TK_var.u.GSD", quantile(TK_var.u.GSD, prob=c(0.5, 0.05, 0.95)),"\n")

### Intra-species TD
TD_var.df <- fread(file.path(datafolder,"DON-EC10_pop_samples.csv"))
TD_var.GSD <- TD_var.df$EC10.GSD[sample.int(4000,size = 10000,replace=TRUE)]
cat("TD_var.GSD",quantile(TD_var.GSD,prob=c(0.5,0.05,0.95)),"\n")

### Intra-species TK+TD
Intra.GSD <- exp(sqrt(log(TK_var.u.GSD)^2+log(TD_var.GSD)^2))
AF_intraTKTD.be <- Intra.GSD^qnorm(0.99)
cat("AF_intraTKTD.be",quantile(AF_intraTKTD.be,prob=c(0.5,0.05,0.95)),"\n")

### BE
#Vd=1.24
BE <- bmd / (AUC_h.GM*AF_interTD*AF_intraTKTD.be)
cat("BE",quantile(BE,prob=c(0.5,0.05,0.95)),"\n")

BE.df <- data.frame(bmd,GM_AUC, AUC_h.GM,AF_interTD,TK_var.u.GSD,
                    TD_var.GSD, AF_intraTKTD.be,BE)
fwrite(BE.df,"BE.samples.csv")
