library(data.table)
library(ggpubr)

datafolder <- "BBMD plot"
female_mice <- fread(file.path(datafolder,"iverson-f-et-al-1995-female-parameters.csv"))
male_mice <- fread(file.path(datafolder,"iverson-f-et-al-1995-male-parameters.csv"))

#human exposure dose levels
d_H <- c(0,10^seq(-4,1,0.1)) # revised

# incidence of BMR=5%
HDMI.df <- fread("HDMI.samples.csv")
# This is dose (back in original units) where median person has effect = BMR
HDM50 <- (HDMI.df$bmd/1000) / (HDMI.df$AF_interTK*HDMI.df$AF_interTD)
incBMR <- numeric() # incidence of BMR
for (j in 1:length(d_H)) {
  # This is the fraction of population affected at dose d_H[j]
  incBMR <- cbind(incBMR,plnorm(d_H[j],meanlog = log(HDM50), sdlog = log(Intra_var.GSD)))
}
incBMR <- as.data.frame(incBMR)
names(incBMR) <- paste("dose",1:length(d_H),sep=".")
incBMRconf <- as.data.frame(t(apply(incBMR,2,quantile,prob=c(0.05,0.5,0.95))))
incBMRconf$d_H <- d_H


incBMR5 <- ggplot(incBMRconf)+
  geom_ribbon(aes(x=d_H,ymin =`5%`,ymax=`95%`),alpha=0.2)+
  geom_line(aes(x=d_H,y=`50%`))+
  geom_line(aes(x=d_H,y=`5%`),linetype="dotted",alpha=0.5)+
  geom_line(aes(x=d_H,y=`95%`),linetype="dotted",alpha=0.5)+
  xlab("Dose (mg/kg/day)")+ylab("Incidence")+
  theme_bw()+
  coord_cartesian(xlim=c(1e-3,1),ylim=c(0,1))+
  scale_x_log10()+
  theme(text = element_text(size = 16))
print(incBMR5)

incBMR5.log <- ggplot(incBMRconf)+
  geom_ribbon(aes(x=d_H,ymin =`5%`,ymax=`95%`),alpha=0.2)+
  geom_line(aes(x=d_H,y=`50%`))+
  geom_line(aes(x=d_H,y=`5%`),linetype="dotted",alpha=0.5)+
  geom_line(aes(x=d_H,y=`95%`),linetype="dotted",alpha=0.5)+
  xlab("Dose (mg/kg/day)")+ylab("Incidence")+
  theme_bw()+
  coord_cartesian(xlim=c(1e-3,1),ylim=c(1e-4,1))+
  scale_x_log10()+scale_y_log10()+
  theme(text = element_text(size = 16))
print(incBMR5.log)

# dose-response (fractional change) at different population percentiles

#set function 
frachange.1 <- function(a,b,dose, dose.max = 1.5){
  resp <- (a*exp(b*dose/dose.max))/(a*exp(b*0))-1
  return(resp)
}
frachange.2 <- function(a,b,g,dose, dose.max = 1.5){
  resp <- (a*exp(b*((dose/dose.max)^g)))/(a*exp(b*(0^g)))-1
  return(resp)
}
frachange.3 <- function(a,b,c,dose, dose.max = 1.5){
  resp <- (a*(c-(c-1)*exp(-b*dose/dose.max)))/(a*(c-(c-1)*exp(-b*0)))-1
  return(resp)
}
frachange.4 <- function(a,b,c,g,dose, dose.max = 1.5){
  resp <- (a*(c-(c-1)*exp(-(b*dose/dose.max)^g)))/(a*(c-(c-1)*exp(-(b*0)^g)))-1
  return(resp)
}
frachange.5 <- function(a,b,c,g,dose, dose.max = 1.5){
  resp <- (a+(b*(dose/dose.max)^g)/(c^g+(dose/dose.max)^g))/(a+(b*0^g)/(c^g+0^g))-1
  return(resp)
}
frachange.6 <- function(a,b,g,dose, dose.max = 1.5){
  resp <- (a+b*((dose/dose.max)^g))/(a+b*(0^g))-1
  return(resp)
}
frachange.7 <- function(a,b,c,dose, dose.max = 1.5){
  resp <- (a+(b*dose/dose.max)/(c+dose/dose.max))/(a+(b*0)/(c+0))-1
  return(resp)
}
frachange.8 <- function(a,b,dose, dose.max = 1.5){
  resp <- (a+b*dose/dose.max)/(a+b*0)-1
  return(resp)
}

#Individual percentiles

yconf.df <- data.frame()
for (inc in c(0.01,0.5,0.99)) { # individual percentiles
  z_i <- qnorm(1-inc) # percentile
  y<-numeric() # blank 
  for (j in 1:length(d_H)) {
    AED <- d_H[j] * AF_interTK * AF_interTD * (Intra_var.GSD^z_i)
    
    #model weight
    w_k_female <- c(0.271, 0.077, 0.216, 0.052, 0.058, 0.031, 0.174, 0.12)
    w_k_female <- w_k_female/sum(w_k_female) # Make sure adds up to 1
    w_k_female_cum <- cumsum(w_k_female) # Cumulative
    w_k_male <- c(0.224, 0.12, 0.09, 0.064, 0.065, 0.11, 0.058, 0.27)
    w_k_male <- w_k_male/sum(w_k_male) # Make sure adds up to 1
    w_k_male_cum <- cumsum(w_k_male) # Cumulative
    
    w_k5000_female_cum <- c(0,round(5000*w_k_female_cum)) # start/end index numbers
    w_k5000_male_cum <- 5000+c(0,round(5000*w_k_male_cum)) # start/end index numbers
    
    w_k5000_female <- diff(w_k5000_female_cum) # Difference is size
    w_k5000_male <- diff(w_k5000_male_cum) # Difference is size
    
    samp.parms.1 <- female_mice[sample(c(1:5000), size=w_k5000_female[1]),]
    samp.parms.2 <- female_mice[sample(c(5001:10000), size=w_k5000_female[2]),]
    samp.parms.3 <- female_mice[sample(c(10001:15000), size=w_k5000_female[3]),]
    samp.parms.4 <- female_mice[sample(c(15001:20000), size=w_k5000_female[4]),]
    samp.parms.5 <- female_mice[sample(c(20001:25000), size=w_k5000_female[5]),]
    samp.parms.6 <- female_mice[sample(c(25001:30000), size=w_k5000_female[6]),]
    samp.parms.7 <- female_mice[sample(c(30001:35000), size=w_k5000_female[7]),]
    samp.parms.8 <- female_mice[sample(c(35001:40000), size=w_k5000_female[8]),]

    # Use subsamples of AED
    y_female <- c(
      frachange.1(a=samp.parms.1$a, b=samp.parms.1$b, dose=AED[(1+w_k5000_female_cum[1]):w_k5000_female_cum[2]], dose.max=1.5)
      ,frachange.2(a=samp.parms.2$a, b=samp.parms.2$b, g=samp.parms.2$g, dose=AED[(1+w_k5000_female_cum[2]):w_k5000_female_cum[3]], dose.max=1.5)
      ,frachange.3(a=samp.parms.3$a, b=samp.parms.3$b, c=samp.parms.3$c, dose=AED[(1+w_k5000_female_cum[3]):w_k5000_female_cum[4]], dose.max=1.5)
      ,frachange.4(a=samp.parms.4$a, b=samp.parms.4$b, c=samp.parms.4$c, g=samp.parms.4$g, dose=AED[(1+w_k5000_female_cum[4]):w_k5000_female_cum[5]], dose.max=1.5)
      ,frachange.5(a=samp.parms.5$a, b=samp.parms.5$b, c=samp.parms.5$c, g=samp.parms.5$g, dose=AED[(1+w_k5000_female_cum[5]):w_k5000_female_cum[6]], dose.max=1.5)
      ,frachange.6(a=samp.parms.6$a, b=samp.parms.6$b, g=samp.parms.6$g, dose=AED[(1+w_k5000_female_cum[6]):w_k5000_female_cum[7]], dose.max=1.5)
      ,frachange.7(a=samp.parms.7$a, b=samp.parms.7$b, c=samp.parms.7$c, dose=AED[(1+w_k5000_female_cum[7]):w_k5000_female_cum[8]], dose.max=1.5)
      ,frachange.8(a=samp.parms.8$a, b=samp.parms.8$b, dose=AED[(1+w_k5000_female_cum[8]):w_k5000_female_cum[9]], dose.max=1.5)
    )
    
    samp.parms.1 <- male_mice[sample(c(1:5000), size=w_k5000_male[1]),]
    samp.parms.2 <- male_mice[sample(c(5001:10000), size=w_k5000_male[2]),]
    samp.parms.3 <- male_mice[sample(c(10001:15000), size=w_k5000_male[3]),]
    samp.parms.4 <- male_mice[sample(c(15001:20000), size=w_k5000_male[4]),]
    samp.parms.5 <- male_mice[sample(c(20001:25000), size=w_k5000_male[5]),]
    samp.parms.6 <- male_mice[sample(c(25001:30000), size=w_k5000_male[6]),]
    samp.parms.7 <- male_mice[sample(c(30001:35000), size=w_k5000_male[7]),]
    samp.parms.8 <- male_mice[sample(c(35001:40000), size=w_k5000_male[8]),]
    
    y_male <- c(
      frachange.1(a=samp.parms.1$a, b=samp.parms.1$b, dose=AED[(1+w_k5000_male_cum[1]):w_k5000_male_cum[2]], dose.max=1.1)
      ,frachange.2(a=samp.parms.2$a, b=samp.parms.2$b, g=samp.parms.2$g, dose=AED[(1+w_k5000_male_cum[2]):w_k5000_male_cum[3]], dose.max=1.1)
      ,frachange.3(a=samp.parms.3$a, b=samp.parms.3$b, c=samp.parms.3$c, dose=AED[(1+w_k5000_male_cum[3]):w_k5000_male_cum[4]], dose.max=1.1)
      ,frachange.4(a=samp.parms.4$a, b=samp.parms.4$b, c=samp.parms.4$c, g=samp.parms.4$g, dose=AED[(1+w_k5000_male_cum[4]):w_k5000_male_cum[5]], dose.max=1.1)
      ,frachange.5(a=samp.parms.5$a, b=samp.parms.5$b, c=samp.parms.5$c, g=samp.parms.5$g, dose=AED[(1+w_k5000_male_cum[5]):w_k5000_male_cum[6]], dose.max=1.1)
      ,frachange.6(a=samp.parms.6$a, b=samp.parms.6$b, g=samp.parms.6$g, dose=AED[(1+w_k5000_male_cum[6]):w_k5000_male_cum[7]], dose.max=1.1)
      ,frachange.7(a=samp.parms.7$a, b=samp.parms.7$b, c=samp.parms.7$c, dose=AED[(1+w_k5000_male_cum[7]):w_k5000_male_cum[8]], dose.max=1.1)
      ,frachange.8(a=samp.parms.8$a, b=samp.parms.8$b, dose=AED[(1+w_k5000_male_cum[8]):w_k5000_male_cum[9]], dose.max=1.1)
    )
    y <- cbind(y,c(y_female,y_male))
  }
  y <- as.data.frame(y)
  names(y) <- paste("dose",1:length(d_H),sep=".")
  yconf <- as.data.frame(t(apply(y,2,quantile,prob=c(0.05,0.5,0.95))))
  yconf$d_H <- d_H
  yconf$inc <- inc
  yconf.df <- rbind(yconf.df,yconf)
}
yconf.df$poppct <- factor(paste0(100*yconf.df$inc,"%ile"))

percent01 <- ggplot(subset(yconf.df, poppct %in% "1%ile"))+
  geom_ribbon(aes(x=d_H,ymin =`5%`,ymax=`95%`),alpha=0.2)+
  geom_line(aes(x=d_H,y=`50%`))+
  geom_line(aes(x=d_H,y=`5%`),linetype="dotted",alpha=0.5)+
  geom_line(aes(x=d_H,y=`95%`),linetype="dotted",alpha=0.5)+
  theme_bw()+
  ggtitle("1%ile of Population") + 
  xlab("Dose (mg/kg/day)")+ylab("Fractional Change")+
  coord_cartesian(xlim=c(1e-3,1),ylim=c(-0.25,0))+
  scale_x_log10()+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 16))
print(percent01)

percent50 <- ggplot(subset(yconf.df, poppct %in% "50%ile"))+
  geom_ribbon(aes(x=d_H,ymin =`5%`,ymax=`95%`),alpha=0.2)+
  geom_line(aes(x=d_H,y=`50%`))+
  geom_line(aes(x=d_H,y=`5%`),linetype="dotted",alpha=0.5)+
  geom_line(aes(x=d_H,y=`95%`),linetype="dotted",alpha=0.5)+
  theme_bw()+
  ggtitle("50%ile of Population") + 
  xlab("Dose (mg/kg/day)")+ylab("Fractional Change")+
  coord_cartesian(xlim=c(1e-3,1),ylim=c(-0.25,0))+
  scale_x_log10()+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 16))
print(percent50)

percent99 <- ggplot(subset(yconf.df, poppct %in% "99%ile"))+
  geom_ribbon(aes(x=d_H,ymin =`5%`,ymax=`95%`),alpha=0.2)+
  geom_line(aes(x=d_H,y=`50%`))+
  geom_line(aes(x=d_H,y=`5%`),linetype="dotted",alpha=0.5)+
  geom_line(aes(x=d_H,y=`95%`),linetype="dotted",alpha=0.5)+
  theme_bw()+
  ggtitle("99%ile of Population") + 
  xlab("Dose (mg/kg/day)")+ylab("Fractional Change")+
  coord_cartesian(xlim=c(1e-3,1),ylim=c(-0.25,0))+
  scale_x_log10()+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 16))
print(percent99)

#ggplot(yconf.df)+
#  geom_ribbon(aes(x=d_H,ymin =`5%`,ymax=`95%`,fill=poppct),alpha=0.2)+
#  geom_line(aes(x=d_H,y=`50%`,group=poppct))+
#  geom_line(aes(x=d_H,y=`5%`,color=poppct),linetype="dotted",alpha=0.5)+
#  geom_line(aes(x=d_H,y=`95%`,color=poppct),linetype="dotted",alpha=0.5)+
#  xlab("Dose (mg/kg/day)")+ylab("Fractional Change")+
#  coord_cartesian(xlim=c(1e-3,1),ylim=c(-0.25,0))+scale_x_log10()


### Population mean

yconfmean.df <- data.frame()
ymean<-numeric() # blank 
npop <- 100
zrand <- rnorm(npop) # npop random individuals
for (j in 1:length(d_H)) {
  y<-numeric() # blank 
  for (k in 1:npop) {
    z_i <- zrand[k]
    AED <- d_H[j] * AF_interTK * AF_interTD * (Intra_var.GSD^z_i)
    
    #model weight
    w_k_female <- c(0.271, 0.077, 0.216, 0.052, 0.058, 0.031, 0.174, 0.12)
    w_k_female <- w_k_female/sum(w_k_female) # Make sure adds up to 1
    w_k_female_cum <- cumsum(w_k_female) # Cumulative
    w_k_male <- c(0.224, 0.12, 0.09, 0.064, 0.065, 0.11, 0.058, 0.27)
    w_k_male <- w_k_male/sum(w_k_male) # Make sure adds up to 1
    w_k_male_cum <- cumsum(w_k_male) # Cumulative
    
    w_k5000_female_cum <- c(0,round(5000*w_k_female_cum)) # start/end index numbers
    w_k5000_male_cum <- 5000+c(0,round(5000*w_k_male_cum)) # start/end index numbers
    
    w_k5000_female <- diff(w_k5000_female_cum) # Difference is size
    w_k5000_male <- diff(w_k5000_male_cum) # Difference is size
    
    samp.parms.1 <- female_mice[sample(c(1:5000), size=w_k5000_female[1]),]
    samp.parms.2 <- female_mice[sample(c(5001:10000), size=w_k5000_female[2]),]
    samp.parms.3 <- female_mice[sample(c(10001:15000), size=w_k5000_female[3]),]
    samp.parms.4 <- female_mice[sample(c(15001:20000), size=w_k5000_female[4]),]
    samp.parms.5 <- female_mice[sample(c(20001:25000), size=w_k5000_female[5]),]
    samp.parms.6 <- female_mice[sample(c(25001:30000), size=w_k5000_female[6]),]
    samp.parms.7 <- female_mice[sample(c(30001:35000), size=w_k5000_female[7]),]
    samp.parms.8 <- female_mice[sample(c(35001:40000), size=w_k5000_female[8]),]

    # Use subsamples of AED
    y_female <- c(
      frachange.1(a=samp.parms.1$a, b=samp.parms.1$b, dose=AED[(1+w_k5000_female_cum[1]):w_k5000_female_cum[2]])
      ,frachange.2(a=samp.parms.2$a, b=samp.parms.2$b, g=samp.parms.2$g, dose=AED[(1+w_k5000_female_cum[2]):w_k5000_female_cum[3]])
      ,frachange.3(a=samp.parms.3$a, b=samp.parms.3$b, c=samp.parms.3$c, dose=AED[(1+w_k5000_female_cum[3]):w_k5000_female_cum[4]])
      ,frachange.4(a=samp.parms.4$a, b=samp.parms.4$b, c=samp.parms.4$c, g=samp.parms.4$g, dose=AED[(1+w_k5000_female_cum[4]):w_k5000_female_cum[5]])
      ,frachange.5(a=samp.parms.5$a, b=samp.parms.5$b, c=samp.parms.5$c, g=samp.parms.5$g, dose=AED[(1+w_k5000_female_cum[5]):w_k5000_female_cum[6]])
      ,frachange.6(a=samp.parms.6$a, b=samp.parms.6$b, g=samp.parms.6$g, dose=AED[(1+w_k5000_female_cum[6]):w_k5000_female_cum[7]])
      ,frachange.7(a=samp.parms.7$a, b=samp.parms.7$b, c=samp.parms.7$c, dose=AED[(1+w_k5000_female_cum[7]):w_k5000_female_cum[8]])
      ,frachange.8(a=samp.parms.8$a, b=samp.parms.8$b, dose=AED[(1+w_k5000_female_cum[8]):w_k5000_female_cum[9]])
    )
    
    samp.parms.1 <- male_mice[sample(c(1:5000), size=w_k5000_male[1]),]
    samp.parms.2 <- male_mice[sample(c(5001:10000), size=w_k5000_male[2]),]
    samp.parms.3 <- male_mice[sample(c(10001:15000), size=w_k5000_male[3]),]
    samp.parms.4 <- male_mice[sample(c(15001:20000), size=w_k5000_male[4]),]
    samp.parms.5 <- male_mice[sample(c(20001:25000), size=w_k5000_male[5]),]
    samp.parms.6 <- male_mice[sample(c(25001:30000), size=w_k5000_male[6]),]
    samp.parms.7 <- male_mice[sample(c(30001:35000), size=w_k5000_male[7]),]
    samp.parms.8 <- male_mice[sample(c(35001:40000), size=w_k5000_male[8]),]
    
    y_male <- c(
      frachange.1(a=samp.parms.1$a, b=samp.parms.1$b, dose=AED[(1+w_k5000_male_cum[1]):w_k5000_male_cum[2]])
      ,frachange.2(a=samp.parms.2$a, b=samp.parms.2$b, g=samp.parms.2$g, dose=AED[(1+w_k5000_male_cum[2]):w_k5000_male_cum[3]])
      ,frachange.3(a=samp.parms.3$a, b=samp.parms.3$b, c=samp.parms.3$c, dose=AED[(1+w_k5000_male_cum[3]):w_k5000_male_cum[4]])
      ,frachange.4(a=samp.parms.4$a, b=samp.parms.4$b, c=samp.parms.4$c, g=samp.parms.4$g, dose=AED[(1+w_k5000_male_cum[4]):w_k5000_male_cum[5]])
      ,frachange.5(a=samp.parms.5$a, b=samp.parms.5$b, c=samp.parms.5$c, g=samp.parms.5$g, dose=AED[(1+w_k5000_male_cum[5]):w_k5000_male_cum[6]])
      ,frachange.6(a=samp.parms.6$a, b=samp.parms.6$b, g=samp.parms.6$g, dose=AED[(1+w_k5000_male_cum[6]):w_k5000_male_cum[7]])
      ,frachange.7(a=samp.parms.7$a, b=samp.parms.7$b, c=samp.parms.7$c, dose=AED[(1+w_k5000_male_cum[7]):w_k5000_male_cum[8]])
      ,frachange.8(a=samp.parms.8$a, b=samp.parms.8$b, dose=AED[(1+w_k5000_male_cum[8]):w_k5000_male_cum[9]])
    )
    y <- cbind(y,c(y_female,y_male))
  }
  ymean <- cbind(ymean,apply(y,1,mean))
}
ymean<-as.data.frame(ymean)
names(ymean) <- paste("dose",1:length(d_H),sep=".")
yconfmean <- as.data.frame(t(apply(ymean,2,quantile,prob=c(0.05,0.5,0.95))))
yconfmean$d_H <- d_H

popmean <- ggplot(yconfmean)+
  geom_ribbon(aes(x=d_H,ymin =`5%`,ymax=`95%`),alpha=0.2)+
  geom_line(aes(x=d_H,y=`50%`))+
  geom_line(aes(x=d_H,y=`5%`),linetype="dotted",alpha=0.5)+
  geom_line(aes(x=d_H,y=`95%`),linetype="dotted",alpha=0.5)+
  theme_bw()+
  ggtitle("Population Mean") + 
  xlab("Dose (mg/kg/day)")+ylab("Fractional Change")+
  coord_cartesian(xlim=c(1e-3,1),ylim=c(-0.25,0))+
  scale_x_log10()+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 16))
print(popmean)

doseresp.BMR5 <- ggarrange(incBMR5, incBMR5.log, labels = c("A","B"),
                          nrow = 2, ncol = 1,
                          font.label = list(size = 18))
doseresp.pop <- ggarrange(popmean, percent01, percent50, percent99,
                          labels = c("C","D","E","F"),
                          nrow = 4, ncol = 1,
                          font.label = list(size = 18))
doseresp <- ggarrange(doseresp.BMR5, doseresp.pop,
                      nrow=1, ncol=2)
ggsave("doseresp.pdf",plot=doseresp, width=20, height=20)
