library(data.table)
library(ggplot2)
library(ggpubr)

datafolder <- "InputData"
don.urine <- fread(file.path(datafolder,"don exposure.csv"))

#iMOE
don.urine$exp <- don.urine$`PDI DON`*0.64
urine.q <- quantile(don.urine$exp, prob=c(0.5,0.95))
#Median: 0.134 ug/kg-d = 0.000134 mg/kg-d
#P95: 1.014 ug/kg-d = 0.001014 mg/kg-d
#P95/P50=7.567
#GM: 0.000134; GSD: 7.567^(1/1.96)=2.808216
ue <- rlnorm(10000, meanlog = log(0.000134), sdlog=log(2.808216))

moe <- BE.urine/ue
p <- ecdf(log(moe))

p(1) #fraction

pdf("Margin of exposure.pdf", width=10, height=8)
plot(p, main="Margin of Exposure")
abline(v=p(1), col="red")
mtext("fraction = 0.0116", side = 3)

dev.off()

#urine UEi
datafolder <- "BBMD plot"
female_mice <- fread(file.path(datafolder,"don-female-parameters.csv"))
male_mice <- fread(file.path(datafolder,"don-male-parameters.csv"))

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

set.seed(1234)

### Population random sample 
yran<-numeric() # blank 
npop <- length(ue)
zrand <- rnorm(npop) # npop random individuals
for (j in 1:length(ue)) {
  y<-numeric() # blank 
  #for (k in 1:npop) {
    z_i <- zrand[j]
    AED <- (ue[j] * AF_interTD * (Intra.GSD^z_i))/(Cavg_dose_a*U_Cavg_h.GM)
    
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
    
    samp.parms.1 <- female_mice[sample(c(1:15000), size=w_k5000_female[1]),]
    samp.parms.2 <- female_mice[sample(c(15001:30000), size=w_k5000_female[2]),]
    samp.parms.3 <- female_mice[sample(c(30001:45000), size=w_k5000_female[3]),]
    samp.parms.4 <- female_mice[sample(c(45001:60000), size=w_k5000_female[4]),]
    samp.parms.5 <- female_mice[sample(c(60001:75000), size=w_k5000_female[5]),]
    samp.parms.6 <- female_mice[sample(c(75001:90000), size=w_k5000_female[6]),]
    samp.parms.7 <- female_mice[sample(c(90001:105000), size=w_k5000_female[7]),]
    samp.parms.8 <- female_mice[sample(c(105001:120000), size=w_k5000_female[8]),]
    
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
    
    samp.parms.1 <- male_mice[sample(c(1:15000), size=w_k5000_male[1]),]
    samp.parms.2 <- male_mice[sample(c(15001:30000), size=w_k5000_male[2]),]
    samp.parms.3 <- male_mice[sample(c(30001:45000), size=w_k5000_male[3]),]
    samp.parms.4 <- male_mice[sample(c(45001:60000), size=w_k5000_male[4]),]
    samp.parms.5 <- male_mice[sample(c(60001:75000), size=w_k5000_male[5]),]
    samp.parms.6 <- male_mice[sample(c(75001:90000), size=w_k5000_male[6]),]
    samp.parms.7 <- male_mice[sample(c(90001:105000), size=w_k5000_male[7]),]
    samp.parms.8 <- male_mice[sample(c(105001:120000), size=w_k5000_male[8]),]
    
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
    y <- c(y_female,y_male)
  yran <- cbind(yran,y)
}
yran<-as.data.frame(yran)
names(yran) <- paste("ue",1:length(ue),sep=".")
write.csv(yran, file="UE fraction of response.csv", row.names=FALSE)

yranpop <- as.data.frame(apply(yran, 1, sample, size=1))
colnames(yranpop) <- "random"

ymean <- as.data.frame(apply(yran,1,mean))
colnames(ymean) <- "mean"

y1perc <- as.data.frame(apply(yran,1,quantile,prob=0.01))
colnames(y1perc) <- "1%pop"

ymed <- as.data.frame(apply(yran,1,quantile,prob=0.5))
colnames(ymed) <- "50%pop"

y99perc <- as.data.frame(apply(yran,1,quantile,prob=0.99))
colnames(y99perc) <- "99%pop"

yran.df <- cbind(yranpop, ymean, y1perc, ymed, y99perc)
write.csv(yran.df, file="Perc population UE fraction of response.csv", row.names=FALSE)

random <- ggplot(yran.df, aes(x=random))+
  geom_density()+
  theme_bw()+
  ggtitle("random individual")+
  xlab("fraction of response (% BW change)")+
  theme(axis.title.y = element_blank())

mean <- ggplot(yran.df, aes(x=mean))+
  geom_density()+
  theme_bw()+
  ggtitle("population mean")+
  xlab("fraction of response (% BW change)")+
  theme(axis.title.y = element_blank())

perc01 <- ggplot(yran.df, aes(x=`1%pop`))+
  geom_density()+
  theme_bw()+
  ggtitle("1%ile of Population")+
  xlab("fraction of response (% BW change)")+
  theme(axis.title.y = element_blank())

perc50 <- ggplot(yran.df, aes(x=`50%pop`))+
  geom_density()+
  theme_bw()+
  ggtitle("50%ile of Population")+
  xlab("fraction of response (% BW change)")+
  theme(axis.title.y = element_blank())

perc99 <- ggplot(yran.df, aes(x=`99%pop`))+
  geom_density()+
  theme_bw()+
  ggtitle("99%ile of Population")+
  xlab("fraction of response (% BW change)")+
  theme(axis.title.y = element_blank())

frac <- ggarrange(random, mean, perc01, perc50, perc99, nrow=5, ncol=1)
ggsave(frac, file = "Density plot_Fraction of response.pdf", width=8, height=10)

quan <- apply(yran.df, 2, quantile, prob=c(0.5, 0.05, 0.95))
write.csv(quan, file="Quantiles of Perc population UE fraction of response.csv")

##boxplot
library(reshape2)

yran.df <- read.csv("Perc population UE fraction of response.csv")
melt.yran <- melt(yran.df[,1:5])

expboxplot <- ggplot(melt.yran, aes(x=value, y=variable))+
  geom_boxplot(aes(fill=variable),
               outlier.shape=NA,
               show.legend = FALSE)+
  scale_y_discrete(limits=rev)+
  theme_bw()+
  xlab("fraction of response (% BW change)")+
  theme(axis.title.y = element_blank())+
  xlim(-0.0075,0)

ggsave(expboxplot, file="Boxplot UE Fraction of response.pdf", width=8, height=10)
