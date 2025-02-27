library(readr)
library(ggplot2)
library(ggpubr)

dose_mice <- read.csv("iverson dose mice.csv")

##female
female_mice <- read.csv("don-female-parameters.csv")
dose_mice_female <- dose_mice[5:8,]
dose <- data.frame(dose=seq(from=0, to=2, by=0.04))
set.seed(3.14159)

#Exp 2
#y=a*exp(b*dose)
Exp2_female <- as.data.frame(matrix(NA, 51, 15000))
for (i in 1:15000){
  a <- female_mice$a[i]
  b <- female_mice$b[i]
  Exp2_female[,i] <- a*exp(b*dose$dose/1.5)
}
#Plot
Exp2_female_quan <- apply(Exp2_female, 1, quantile, prob = c(0.05,0.5,0.95), na.rm = TRUE)
tExp2_female_quan <- t(Exp2_female_quan)
tExp2_female_quan <- cbind(dose, tExp2_female_quan)
Exp2_femaleplot <- ggplot(tExp2_female_quan, aes(x=dose))+
  geom_ribbon(data=tExp2_female_quan, aes(x=dose, y=tExp2_female_quan$`50%`, 
                                          ymin=tExp2_female_quan$`5%`, 
                                          ymax=tExp2_female_quan$`95%`), alpha=0.2) +
  geom_line(data=tExp2_female_quan, aes(x=dose, y=tExp2_female_quan$`50%`))+
  geom_point(data=dose_mice_female, aes(x=dose, y=mean))+
  geom_errorbar(data=dose_mice_female, aes(x=dose, ymin=mean-SD, ymax=mean+SD)) +
  ggtitle("Exponential 2 (26.9%)") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=16),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12)) +
  xlab("Dose (mg/kg-day)") + ylab("Body weight (g)")+
  scale_x_continuous(breaks=seq(0,2,by=0.25))+
  scale_y_continuous(breaks=seq(20,50,by=5))

#Exp 3
#y=a*exp(b*(dose^g))
Exp3_female <- as.data.frame(matrix(NA, 51, 15000))
for (i in 1:15000){
  a <- female_mice$a[i+15000]
  b <- female_mice$b[i+15000]
  g <- female_mice$g[i+15000]
  Exp3_female[,i] <- a*exp(b*((dose$dose/1.5)^g))
}
#Plot
Exp3_female_quan <- apply(Exp3_female, 1, quantile, prob = c(0.05,0.5,0.95), na.rm = TRUE)
tExp3_female_quan <- t(Exp3_female_quan)
tExp3_female_quan <- cbind(dose, tExp3_female_quan)
Exp3_femaleplot <- ggplot(tExp3_female_quan, aes(x=dose))+
  geom_ribbon(data=tExp3_female_quan, aes(x=dose, y=tExp3_female_quan$`50%`, 
                                          ymin=tExp3_female_quan$`5%`, 
                                          ymax=tExp3_female_quan$`95%`), alpha=0.2) +
  geom_line(data=tExp3_female_quan, aes(x=dose, y=tExp3_female_quan$`50%`))+
  geom_point(data=dose_mice_female, aes(x=dose, y=mean))+
  geom_errorbar(data=dose_mice_female, aes(x=dose, ymin=mean-SD, ymax=mean+SD)) +
  ggtitle("Exponential 3 (7.7%)") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=16),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12)) +
  xlab("Dose (mg/kg-day)") + ylab("Body weight (g)")+
  scale_x_continuous(breaks=seq(0,2,by=0.25))+
  scale_y_continuous(breaks=seq(20,50,by=5))

#Exp 4
#y=a*[c-(c-1)*exp(-b*dose)]
Exp4_female <- as.data.frame(matrix(NA, 51, 15000))
for (i in 1:15000){
  a <- female_mice$a[i+30000]
  b <- female_mice$b[i+30000]
  c <- female_mice$c[i+30000]
  Exp4_female[,i] <- a*(c-(c-1)*exp(-b*(dose$dose/1.5)))
}
#Plot
Exp4_female_quan <- apply(Exp4_female, 1, quantile, prob = c(0.05,0.5,0.95), na.rm = TRUE)
tExp4_female_quan <- t(Exp4_female_quan)
tExp4_female_quan <- cbind(dose, tExp4_female_quan)
Exp4_femaleplot <- ggplot(tExp4_female_quan, aes(x=dose))+
  geom_ribbon(data=tExp4_female_quan, aes(x=dose, y=tExp4_female_quan$`50%`, 
                                          ymin=tExp4_female_quan$`5%`, 
                                          ymax=tExp4_female_quan$`95%`), alpha=0.2) +
  geom_line(data=tExp4_female_quan, aes(x=dose, y=tExp4_female_quan$`50%`))+
  geom_point(data=dose_mice_female, aes(x=dose, y=mean))+
  geom_errorbar(data=dose_mice_female, aes(x=dose, ymin=mean-SD, ymax=mean+SD))+
  ggtitle("Exponential 4 (21.2%)") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=16),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12)) +
  xlab("Dose (mg/kg-day)") + ylab("Body weight (g)")+
  scale_x_continuous(breaks=seq(0,2,by=0.25))+
  scale_y_continuous(breaks=seq(20,50,by=5))

#Exp 5
#y=a*[c-(c-1)*exp(-(b*dose)^g)]
Exp5_female <- as.data.frame(matrix(NA, 51, 15000))
for (i in 1:15000){
  a <- female_mice$a[i+45000]
  b <- female_mice$b[i+45000]
  c <- female_mice$c[i+45000]
  g <- female_mice$g[i+45000]
  Exp5_female[,i] <- a*(c-(c-1)*exp(-(b*(dose$dose/1.5))^g))
}
#Plot
Exp5_female_quan <- apply(Exp5_female, 1, quantile, prob = c(0.05,0.5,0.95), na.rm = TRUE)
tExp5_female_quan <- t(Exp5_female_quan)
tExp5_female_quan <- cbind(dose, tExp5_female_quan)
Exp5_femaleplot <- ggplot(tExp5_female_quan, aes(x=dose))+
  geom_ribbon(data=tExp5_female_quan, aes(x=dose, y=tExp5_female_quan$`50%`, 
                                          ymin=tExp5_female_quan$`5%`, 
                                          ymax=tExp5_female_quan$`95%`), alpha=0.2) +
  geom_line(data=tExp5_female_quan, aes(x=dose, y=tExp5_female_quan$`50%`))+
  geom_point(data=dose_mice_female, aes(x=dose, y=mean))+
  geom_errorbar(data=dose_mice_female, aes(x=dose, ymin=mean-SD, ymax=mean+SD))+
  ggtitle("Exponential 5 (5.2%)") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=16),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12)) +
  xlab("Dose (mg/kg-day)") + ylab("Body weight (g)")+
  scale_x_continuous(breaks=seq(0,2,by=0.25))+
  scale_y_continuous(breaks=seq(20,50,by=5))

#Hill
#y=a+(b*dose^g)/(c^g+dose^g)
Hill_female <- as.data.frame(matrix(NA, 51, 15000))
for (i in 1:15000){
  a <- female_mice$a[i+60000]
  b <- female_mice$b[i+60000]
  c <- female_mice$c[i+60000]
  g <- female_mice$g[i+60000]
  Hill_female[,i] <- a+(b*(dose$dose/1.5)^g)/(c^g+(dose$dose/1.5)^g)
}
#Plot
Hill_female_quan <- apply(Hill_female, 1, quantile, prob = c(0.05,0.5,0.95), na.rm = TRUE)
tHill_female_quan <- t(Hill_female_quan)
tHill_female_quan <- cbind(dose, tHill_female_quan)
Hill_femaleplot <- ggplot(tHill_female_quan, aes(x=dose))+
  geom_ribbon(data=tHill_female_quan, aes(x=dose, y=tHill_female_quan$`50%`, 
                                          ymin=tHill_female_quan$`5%`, 
                                          ymax=tHill_female_quan$`95%`), alpha=0.2) +
  geom_line(data=tHill_female_quan, aes(x=dose, y=tHill_female_quan$`50%`))+
  geom_point(data=dose_mice_female, aes(x=dose, y=mean))+
  geom_errorbar(data=dose_mice_female, aes(x=dose, ymin=mean-SD, ymax=mean+SD))+
  ggtitle("Hill (5.9%)") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=16),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12)) +
  xlab("Dose (mg/kg-day)") + ylab("Body weight (g)")+
  scale_x_continuous(breaks=seq(0,2,by=0.25))+
  scale_y_continuous(breaks=seq(20,50,by=5))

#Power
#y=a+b*(dose^g)
Power_female <- as.data.frame(matrix(NA, 51, 15000))
for (i in 1:15000){
  a <- female_mice$a[i+75000]
  b <- female_mice$b[i+75000]
  g <- female_mice$g[i+75000]
  Power_female[,i] <- a+b*((dose$dose/1.5)^g)
}
#Plot
Power_female_quan <- apply(Power_female, 1, quantile, prob = c(0.05,0.5,0.95), na.rm = TRUE)
tPower_female_quan <- t(Power_female_quan)
tPower_female_quan <- cbind(dose, tPower_female_quan)
Power_femaleplot <- ggplot(tPower_female_quan, aes(x=dose))+
  geom_ribbon(data=tPower_female_quan, aes(x=dose, y=tPower_female_quan$`50%`, 
                                          ymin=tPower_female_quan$`5%`, 
                                          ymax=tPower_female_quan$`95%`), alpha=0.2) +
  geom_line(data=tPower_female_quan, aes(x=dose, y=tPower_female_quan$`50%`))+
  geom_point(data=dose_mice_female, aes(x=dose, y=mean))+
  geom_errorbar(data=dose_mice_female, aes(x=dose, ymin=mean-SD, ymax=mean+SD))+
  ggtitle("Power (3.2%)") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=16),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12)) +
  xlab("Dose (mg/kg-day)") + ylab("Body weight (g)")+
  scale_x_continuous(breaks=seq(0,2,by=0.25))+
  scale_y_continuous(breaks=seq(20,50,by=5))

#MichaelisMenten
#y=a+(b*dose)/(c+dose)
Mich_female <- as.data.frame(matrix(NA, 51, 15000))
for (i in 1:15000){
  a <- female_mice$a[i+90000]
  b <- female_mice$b[i+90000]
  c <- female_mice$c[i+90000]
  Mich_female[,i] <- a+(b*(dose$dose/1.5))/(c+dose$dose/1.5)
}
#Plot
Mich_female_quan <- apply(Mich_female, 1, quantile, prob = c(0.05,0.5,0.95), na.rm = TRUE)
tMich_female_quan <- t(Mich_female_quan)
tMich_female_quan <- cbind(dose, tMich_female_quan)
Mich_femaleplot <- ggplot(tMich_female_quan, aes(x=dose))+
  geom_ribbon(data=tMich_female_quan, aes(x=dose, y=tMich_female_quan$`50%`, 
                                          ymin=tMich_female_quan$`5%`, 
                                          ymax=tMich_female_quan$`95%`), alpha=0.2) +
  geom_line(data=tMich_female_quan, aes(x=dose, y=tMich_female_quan$`50%`))+
  geom_point(data=dose_mice_female, aes(x=dose, y=mean))+
  geom_errorbar(data=dose_mice_female, aes(x=dose, ymin=mean-SD, ymax=mean+SD))+
  ggtitle("Michaelis-Menten (17.9%)") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=16),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12)) +
  xlab("Dose (mg/kg-day)") + ylab("Body weight (g)")+
  scale_x_continuous(breaks=seq(0,2,by=0.25))+
  scale_y_continuous(breaks=seq(20,50,by=5))

#Linear
#y=a+b*dose
Linear_female <- as.data.frame(matrix(NA, 51, 15000))
for (i in 1:15000){
  a <- female_mice$a[i+105000]
  b <- female_mice$b[i+105000]
  Linear_female[,i] <- a+b*dose$dose/1.5
}
#Plot
Linear_female_quan <- apply(Linear_female, 1, quantile, prob = c(0.05,0.5,0.95), na.rm = TRUE)
tLinear_female_quan <- t(Linear_female_quan)
tLinear_female_quan <- cbind(dose, tLinear_female_quan)
Linear_femaleplot <- ggplot(tLinear_female_quan, aes(x=dose))+
  geom_ribbon(data=tLinear_female_quan, aes(x=dose, y=tLinear_female_quan$`50%`, 
                                          ymin=tLinear_female_quan$`5%`, 
                                          ymax=tLinear_female_quan$`95%`), alpha=0.2) +
  geom_line(data=tLinear_female_quan, aes(x=dose, y=tLinear_female_quan$`50%`))+
  geom_point(data=dose_mice_female, aes(x=dose, y=mean))+
  geom_errorbar(data=dose_mice_female, aes(x=dose, ymin=mean-SD, ymax=mean+SD))+
  ggtitle("Linear (12.0%)") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=16),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12)) +
  xlab("Dose (mg/kg-day)") + ylab("Body weight (g)")+
  scale_x_continuous(breaks=seq(0,2,by=0.25))+
  scale_y_continuous(breaks=seq(20,50,by=5))

#Combine plots
female_combine <- ggarrange(Exp2_femaleplot, Exp3_femaleplot, Exp4_femaleplot, Exp5_femaleplot,
         Hill_femaleplot, Power_femaleplot, Mich_femaleplot, Linear_femaleplot,
         ncol=2, nrow=4)
female_combine <- annotate_figure(female_combine, top = text_grob("Female mice", size = 22))

##male
male_mice <- read.csv("don-male-parameters.csv")
dose_mice_male <- dose_mice[1:4,]
dose <- data.frame(dose=seq(from=0, to=2, by=0.04))
set.seed(3.14159)

#Exp 2
#y=a*exp(b*dose)
Exp2_male <- as.data.frame(matrix(NA, 51, 15000))
for (i in 1:15000){
  a <- male_mice$a[i]
  b <- male_mice$b[i]
  Exp2_male[,i] <- a*exp(b*dose$dose/1.1)
}
#Plot
Exp2_male_quan <- apply(Exp2_male, 1, quantile, prob = c(0.05,0.5,0.95), na.rm = TRUE)
tExp2_male_quan <- t(Exp2_male_quan)
tExp2_male_quan <- cbind(dose, tExp2_male_quan)
Exp2_maleplot <- ggplot(tExp2_male_quan, aes(x=dose))+
  geom_ribbon(data=tExp2_male_quan, aes(x=dose, y=tExp2_male_quan$`50%`, 
                                          ymin=tExp2_male_quan$`5%`, 
                                          ymax=tExp2_male_quan$`95%`), alpha=0.2) +
  geom_line(data=tExp2_male_quan, aes(x=dose, y=tExp2_male_quan$`50%`))+
  geom_point(data=dose_mice_male, aes(x=dose, y=mean))+
  geom_errorbar(data=dose_mice_male, aes(x=dose, ymin=mean-SD, ymax=mean+SD)) +
  ggtitle("Exponential 2 (22.4%)") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=16),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12)) +
  xlab("Dose (mg/kg-day)") + ylab("Body weight (g)")+
  scale_x_continuous(breaks=seq(0,2,by=0.25))+
  scale_y_continuous(limits=c(25,47),breaks=c(25,30,35,40,45),labels=c(25,30,35,40,45))

#Exp 3
#y=a*exp(b*(dose^g))
Exp3_male <- as.data.frame(matrix(NA, 51, 15000))
for (i in 1:15000){
  a <- male_mice$a[i+15000]
  b <- male_mice$b[i+15000]
  g <- male_mice$g[i+15000]
  Exp3_male[,i] <- a*exp(b*((dose$dose/1.1)^g))
}
#Plot
Exp3_male_quan <- apply(Exp3_male, 1, quantile, prob = c(0.05,0.5,0.95), na.rm = TRUE)
tExp3_male_quan <- t(Exp3_male_quan)
tExp3_male_quan <- cbind(dose, tExp3_male_quan)
Exp3_maleplot <- ggplot(tExp3_male_quan, aes(x=dose))+
  geom_ribbon(data=tExp3_male_quan, aes(x=dose, y=tExp3_male_quan$`50%`, 
                                          ymin=tExp3_male_quan$`5%`, 
                                          ymax=tExp3_male_quan$`95%`), alpha=0.2) +
  geom_line(data=tExp3_male_quan, aes(x=dose, y=tExp3_male_quan$`50%`))+
  geom_point(data=dose_mice_male, aes(x=dose, y=mean))+
  geom_errorbar(data=dose_mice_male, aes(x=dose, ymin=mean-SD, ymax=mean+SD)) +
  ggtitle("Exponential 3 (12.1%)") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=16),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12)) +
  xlab("Dose (mg/kg-day)") + ylab("Body weight (g)")+
  scale_x_continuous(breaks=seq(0,2,by=0.25))+
  scale_y_continuous(limits=c(23,47),breaks=c(25,30,35,40,45),labels=c(25,30,35,40,45))

#Exp 4
#y=a*[c-(c-1)*exp(-b*dose)]
Exp4_male <- as.data.frame(matrix(NA, 51, 15000))
for (i in 1:15000){
  a <- male_mice$a[i+30000]
  b <- male_mice$b[i+30000]
  c <- male_mice$c[i+30000]
  Exp4_male[,i] <- a*(c-(c-1)*exp(-b*(dose$dose/1.1)))
}
#Plot
Exp4_male_quan <- apply(Exp4_male, 1, quantile, prob = c(0.05,0.5,0.95), na.rm = TRUE)
tExp4_male_quan <- t(Exp4_male_quan)
tExp4_male_quan <- cbind(dose, tExp4_male_quan)
Exp4_maleplot <- ggplot(tExp4_male_quan, aes(x=dose))+
  geom_ribbon(data=tExp4_male_quan, aes(x=dose, y=tExp4_male_quan$`50%`, 
                                          ymin=tExp4_male_quan$`5%`, 
                                          ymax=tExp4_male_quan$`95%`), alpha=0.2) +
  geom_line(data=tExp4_male_quan, aes(x=dose, y=tExp4_male_quan$`50%`))+
  geom_point(data=dose_mice_male, aes(x=dose, y=mean))+
  geom_errorbar(data=dose_mice_male, aes(x=dose, ymin=mean-SD, ymax=mean+SD))+
  ggtitle("Exponential 4 (8.9%)") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=16),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12)) +
  xlab("Dose (mg/kg-day)") + ylab("Body weight (g)")+
  scale_x_continuous(breaks=seq(0,2,by=0.25))+
  scale_y_continuous(limits=c(25,47),breaks=c(25,30,35,40,45),labels=c(25,30,35,40,45))

#Exp 5
#y=a*[c-(c-1)*exp(-(b*dose)^g)]
Exp5_male <- as.data.frame(matrix(NA, 51, 15000))
for (i in 1:15000){
  a <- male_mice$a[i+45000]
  b <- male_mice$b[i+45000]
  c <- male_mice$c[i+45000]
  g <- male_mice$g[i+45000]
  Exp5_male[,i] <- a*(c-(c-1)*exp(-(b*(dose$dose/1.1))^g))
}
#Plot
Exp5_male_quan <- apply(Exp5_male, 1, quantile, prob = c(0.05,0.5,0.95), na.rm = TRUE)
tExp5_male_quan <- t(Exp5_male_quan)
tExp5_male_quan <- cbind(dose, tExp5_male_quan)
Exp5_maleplot <- ggplot(tExp5_male_quan, aes(x=dose))+
  geom_ribbon(data=tExp5_male_quan, aes(x=dose, y=tExp5_male_quan$`50%`, 
                                          ymin=tExp5_male_quan$`5%`, 
                                          ymax=tExp5_male_quan$`95%`), alpha=0.2) +
  geom_line(data=tExp5_male_quan, aes(x=dose, y=tExp5_male_quan$`50%`))+
  geom_point(data=dose_mice_male, aes(x=dose, y=mean))+
  geom_errorbar(data=dose_mice_male, aes(x=dose, ymin=mean-SD, ymax=mean+SD))+
  ggtitle("Exponential 5 (6.3%)") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=16),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12)) +
  xlab("Dose (mg/kg-day)") + ylab("Body weight (g)")+
  scale_x_continuous(breaks=seq(0,2,by=0.25))+
  scale_y_continuous(limits=c(25,47),breaks=c(25,30,35,40,45),labels=c(25,30,35,40,45))

#Hill
#y=a+(b*dose^g)/(c^g+dose^g)
Hill_male <- as.data.frame(matrix(NA, 51, 15000))
for (i in 1:15000){
  a <- male_mice$a[i+60000]
  b <- male_mice$b[i+60000]
  c <- male_mice$c[i+60000]
  g <- male_mice$g[i+60000]
  Hill_male[,i] <- a+(b*(dose$dose/1.1)^g)/(c^g+(dose$dose/1.1)^g)
}
#Plot
Hill_male_quan <- apply(Hill_male, 1, quantile, prob = c(0.05,0.5,0.95), na.rm = TRUE)
tHill_male_quan <- t(Hill_male_quan)
tHill_male_quan <- cbind(dose, tHill_male_quan)
Hill_maleplot <- ggplot(tHill_male_quan, aes(x=dose))+
  geom_ribbon(data=tHill_male_quan, aes(x=dose, y=tHill_male_quan$`50%`, 
                                          ymin=tHill_male_quan$`5%`, 
                                          ymax=tHill_male_quan$`95%`), alpha=0.2) +
  geom_line(data=tHill_male_quan, aes(x=dose, y=tHill_male_quan$`50%`))+
  geom_point(data=dose_mice_male, aes(x=dose, y=mean))+
  geom_errorbar(data=dose_mice_male, aes(x=dose, ymin=mean-SD, ymax=mean+SD))+
  ggtitle("Hill (6.6%)") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=16),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12)) +
  xlab("Dose (mg/kg-day)") + ylab("Body weight (g)")+
  scale_x_continuous(breaks=seq(0,2,by=0.25))+
  scale_y_continuous(limits=c(25,47),breaks=c(25,30,35,40,45),labels=c(25,30,35,40,45))

#Power
#y=a+b*(dose^g)
Power_male <- as.data.frame(matrix(NA, 51, 15000))
for (i in 1:15000){
  a <- male_mice$a[i+75000]
  b <- male_mice$b[i+75000]
  g <- male_mice$g[i+75000]
  Power_male[,i] <- a+b*((dose$dose/1.1)^g)
}
#Plot
Power_male_quan <- apply(Power_male, 1, quantile, prob = c(0.05,0.5,0.95), na.rm = TRUE)
tPower_male_quan <- t(Power_male_quan)
tPower_male_quan <- cbind(dose, tPower_male_quan)
Power_maleplot <- ggplot(tPower_male_quan, aes(x=dose))+
  geom_ribbon(data=tPower_male_quan, aes(x=dose, y=tPower_male_quan$`50%`, 
                                           ymin=tPower_male_quan$`5%`, 
                                           ymax=tPower_male_quan$`95%`), alpha=0.2) +
  geom_line(data=tPower_male_quan, aes(x=dose, y=tPower_male_quan$`50%`))+
  geom_point(data=dose_mice_male, aes(x=dose, y=mean))+
  geom_errorbar(data=dose_mice_male, aes(x=dose, ymin=mean-SD, ymax=mean+SD))+
  ggtitle("Power (11.0%)") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=16),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12)) +
  xlab("Dose (mg/kg-day)") + ylab("Body weight (g)")+
  scale_x_continuous(breaks=seq(0,2,by=0.25))+
  scale_y_continuous(limits=c(20,47),breaks=c(25,30,35,40,45),labels=c(25,30,35,40,45))

#MichaelisMenten
#y=a+(b*dose)/(c+dose)
Mich_male <- as.data.frame(matrix(NA, 51, 15000))
for (i in 1:15000){
  a <- male_mice$a[i+90000]
  b <- male_mice$b[i+90000]
  c <- male_mice$c[i+90000]
  Mich_male[,i] <- a+(b*(dose$dose/1.1))/(c+dose$dose/1.1)
}
#Plot
Mich_male_quan <- apply(Mich_male, 1, quantile, prob = c(0.05,0.5,0.95), na.rm = TRUE)
tMich_male_quan <- t(Mich_male_quan)
tMich_male_quan <- cbind(dose, tMich_male_quan)
Mich_maleplot <- ggplot(tMich_male_quan, aes(x=dose))+
  geom_ribbon(data=tMich_male_quan, aes(x=dose, y=tMich_male_quan$`50%`, 
                                          ymin=tMich_male_quan$`5%`, 
                                          ymax=tMich_male_quan$`95%`), alpha=0.2) +
  geom_line(data=tMich_male_quan, aes(x=dose, y=tMich_male_quan$`50%`))+
  geom_point(data=dose_mice_male, aes(x=dose, y=mean))+
  geom_errorbar(data=dose_mice_male, aes(x=dose, ymin=mean-SD, ymax=mean+SD))+
  ggtitle("Michaelis-Menten (5.8%)") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=16),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12)) +
  xlab("Dose (mg/kg-day)") + ylab("Body weight (g)")+
  scale_x_continuous(breaks=seq(0,2,by=0.25))+
  scale_y_continuous(limits=c(25,47),breaks=c(25,30,35,40,45),labels=c(25,30,35,40,45))

#Linear
#y=a+b*dose
Linear_male <- as.data.frame(matrix(NA, 51, 15000))
for (i in 1:15000){
  a <- male_mice$a[i+105000]
  b <- male_mice$b[i+105000]
  Linear_male[,i] <- a+b*dose$dose/1.1
}
#Plot
Linear_male_quan <- apply(Linear_male, 1, quantile, prob = c(0.05,0.5,0.95), na.rm = TRUE)
tLinear_male_quan <- t(Linear_male_quan)
tLinear_male_quan <- cbind(dose, tLinear_male_quan)
Linear_maleplot <- ggplot(tLinear_male_quan, aes(x=dose))+
  geom_ribbon(data=tLinear_male_quan, aes(x=dose, y=tLinear_male_quan$`50%`, 
                                            ymin=tLinear_male_quan$`5%`, 
                                            ymax=tLinear_male_quan$`95%`), alpha=0.2) +
  geom_line(data=tLinear_male_quan, aes(x=dose, y=tLinear_male_quan$`50%`))+
  geom_point(data=dose_mice_male, aes(x=dose, y=mean))+
  geom_errorbar(data=dose_mice_male, aes(x=dose, ymin=mean-SD, ymax=mean+SD))+
  ggtitle("Linear (26.9%)") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=16),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12)) +
  xlab("Dose (mg/kg-day)") + ylab("Body weight (g)")+
  scale_x_continuous(breaks=seq(0,2,by=0.25))+
  scale_y_continuous(limits=c(25,47),breaks=c(25,30,35,40,45),labels=c(25,30,35,40,45))

#Combine plots
male_combine <- ggarrange(Exp2_maleplot, Exp3_maleplot, Exp4_maleplot, Exp5_maleplot,
                            Hill_maleplot, Power_maleplot, Mich_maleplot, Linear_maleplot,
                            ncol=2, nrow=4)
male_combine <- annotate_figure(male_combine, top = text_grob("Male mice", size = 22))

#combine male and female models
doseresponse <- ggarrange(male_combine, female_combine, ncol = 2, nrow = 1)

#posterior model average
bmds_f <- read.csv("don-female-bmds.csv")
bmds_m <- read.csv("don-male-bmds.csv")

quantile(bmds_f$model_average, prob=c(0.5, 0.05, 0.95))
quantile(bmds_m$model_average, prob=c(0.5, 0.05, 0.95))

bmds.average_f <- ggplot(bmds_f, aes(x=model_average)) + geom_density() +
  ggtitle("Model Average") + 
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=16),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12)) +
  annotate("text", x=0.22, y=1, label= "0.22 [0.14-0.32]", size=8) +
  annotate("pointrange", x = 0.22, y = 0, xmin = quantile(bmds_f$model_average, prob=c(0.05)), xmax = quantile(bmds_f$model_average, prob=c(0.95)), size =1) +
  scale_x_continuous(limit=c(0,0.5), breaks=seq(0,0.5,by=0.1))+
  scale_y_continuous(limits=c(0,12.5), breaks=seq(0,12,by=2))+
  xlab("Dose (mg/kg-day)")+ylab("Density")

bmds.average_m <- ggplot(bmds_m, aes(x=model_average)) + geom_density() +
  ggtitle("Model Average") + 
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=16),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12)) +
  annotate("text", x=0.3, y=1, label= "0.27 [0.22-0.42]", size=8) +
  annotate("pointrange", x = 0.27, y = 0, xmin = quantile(bmds_m$model_average, prob=c(0.05)), xmax = quantile(bmds_m$model_average, prob=c(0.95)), size =1) +
  scale_x_continuous(limit=c(0,0.5), breaks=seq(0,0.5,by=0.1))+
  scale_y_continuous(limits=c(0,12.5), breaks=seq(0,12,by=2))+
  xlab("Dose (mg/kg-day)")+ylab("Density")

modelaverage <- ggarrange(bmds.average_m, bmds.average_f, ncol = 2, nrow = 1, align="hv")

#combine into one graph
BBMD <- ggarrange(doseresponse, modelaverage, ncol=1, nrow=2, heights = c(12, 6), align="hv")
ggsave("BBMD.pdf",plot=BBMD, width=20, height=18, scale=0.8)
