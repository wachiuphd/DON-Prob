datafolder <- "BBMD plot"
female_mice <- fread(file.path(datafolder,"iverson-f-et-al-1995-female-parameters.csv"))
male_mice <- fread(file.path(datafolder,"iverson-f-et-al-1995-male-parameters.csv"))

#incidence of effect
#AED_i <- d_H * AF_interTK * Intra_var.GSD^qnorm(0.99)

#set function
frachange.1 <- function(a,b,dose){
  resp <- (a*exp(b*dose))/(a*exp(b*0))-1
  return(resp)
}
frachange.2 <- function(a,b,g,dose){
  resp <- (a*exp(b*(dose^g)))/(a*exp(b*(0^g)))-1
  return(resp)
}
frachange.3 <- function(a,b,c,dose){
  resp <- (a*(c-(c-1)*exp(-b*dose)))/(a*(c-(c-1)*exp(-b*0)))-1
  return(resp)
}
frachange.4 <- function(a,b,c,g,dose){
  resp <- (a*(c-(c-1)*exp(-(b*dose)^g)))/(a*(c-(c-1)*exp(-(b*0)^g)))-1
  return(resp)
}
frachange.5 <- function(a,b,c,g,dose){
  resp <- (a+(b*dose^g)/(c^g+dose^g))/(a+(b*0^g)/(c^g+0^g))-1
  return(resp)
}
frachange.6 <- function(a,b,g,dose){
  resp <- (a+b*(dose^g))/(a+b*(0^g))-1
  return(resp)
}
frachange.7 <- function(a,b,c,dose){
  resp <- (a+(b*dose)/(c+dose))/(a+(b*0)/(c+0))-1
  return(resp)
}
frachange.8 <- function(a,b,dose){
  resp <- (a+b*dose)/(a+b*0)-1
  return(resp)
}

#human exposure dose
d_H <- seq(0,10)

AED <- d_H[2] * AF_interTK * AF_interTD * (Intra_var.GSD^qnorm(0.99))

#model weight
w_k_female <- c(0.271, 0.077, 0.216, 0.052, 0.058, 0.031, 0.174, 0.12)
w_k_male <- c(0.224, 0.12, 0.09, 0.064, 0.065, 0.11, 0.058, 0.27)
w_k5000_female <- 5000*w_k_female
w_k5000_male <- 5000*w_k_male

df_female <- matrix(NA, nrow=10000, ncol=5000)

samp.parms.1 <- female_mice[sample(c(1:5000), size=w_k5000_female[1]),]
df_female[c(1:10000), c(1:1355)] <- frachange.1(a=samp.parms.1$a, b=samp.parms.1$b, dose=AED)

samp.parms.2 <- female_mice[sample(c(5001:10000), size=w_k5000_female[2]),]
df_female[c(1:11), c(1356:1740)] <- frachange.2(a=samp.parms.2$a, b=samp.parms.2$b, g=samp.parms.2$g, dose=AED)
samp.parms.3 <- female_mice[sample(c(10001:15000), size=w_k5000_female[3]),]
df_female[c(1:11), c(1741:2820)] <- frachange.3(a=samp.parms.3$a, b=samp.parms.3$b, c=samp.parms.3$c, dose=AED)
samp.parms.4 <- female_mice[sample(c(15001:20000), size=w_k5000_female[4]),]
df_female[c(1:11), c(2821:3080)] <- frachange.4(a=samp.parms.4$a, b=samp.parms.4$b, c=samp.parms.4$c, g=samp.parms.4$g, dose=AED)
samp.parms.5 <- female_mice[sample(c(20001:25000), size=w_k5000_female[5]),]
df_female[c(1:11), c(3081:3370)] <- frachange.5(a=samp.parms.5$a, b=samp.parms.5$b, c=samp.parms.5$c, g=samp.parms.5$g, dose=AED)
samp.parms.6 <- female_mice[sample(c(25001:30000), size=w_k5000_female[6]),]
df_female[c(1:11), c(3371:3525)] <- frachange.6(a=samp.parms.6$a, b=samp.parms.6$b, g=samp.parms.6$g, dose=AED)
samp.parms.7 <- female_mice[sample(c(30001:35000), size=w_k5000_female[7]),]
df_female[c(1:11), c(3526:4394)] <- frachange.7(a=samp.parms.7$a, b=samp.parms.7$b, c=samp.parms.7$c, dose=AED)
samp.parms.8 <- female_mice[sample(c(35001:40000), size=w_k5000_female[8]),]
df_female[c(1:11), c(4396:4995)] <- frachange.8(a=samp.parms.8$a, b=samp.parms.8$b, dose=AED)

df_female_quan <- apply(df_female, 1, quantile, prob = c(0.05,0.5,0.95), na.rm = TRUE)
df_female_quan <- t(df_female_quan)
