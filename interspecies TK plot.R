library(ggplot2)
dat <- data.frame(bw=c(0.135,27),auc=c(0.4881,1.972))

res <- lm(log10(auc) ~ log10(bw), data=dat)
#res <- lm(auc ~ bw, data=dat)

df <- data.frame(x=c(0.135, 27, 0.0316, 0.0246), y=c(0.4881, 1.972, 1.15794, 0.62), type=c("free","free","total","elisa"))
  

#P95/P50: 0.62/0.312 = 1.987179


plot <- ggplot() + 
  geom_point(data=df, aes(x=x, y=y, shape=type),size=4) +
  geom_abline(slope = coef(res)[["log10(bw)"]], intercept = coef(res)[["(Intercept)"]])+
  geom_abline(slope = coef(res)[["log10(bw)"]], intercept = coef(res)[["(Intercept)"]]+log10(1.98), linetype="dashed")+
  geom_abline(slope = coef(res)[["log10(bw)"]], intercept = coef(res)[["(Intercept)"]]-log10(1.98), linetype="dashed")+
  geom_vline(xintercept=0.0363, color="red")+ ##Iverson et al. 1995 study: BW=(0.0373+0.0353)/2=0.0363
  theme_bw()+
  scale_x_log10(limits=c(0.01,30))+  scale_y_log10(limits=c(0.1,10))+
  xlab("Body weight (kg)")+ylab("AUC/dose (mg-h/L per mg/kg-day)")+
  scale_shape_manual(name = "DON measurements",
                     labels = c("Total DON", "Free DON", "ELISA method"),
                     values = c("total"=18, "free"=16, "elisa"=15))
print(plot)

ggsave(plot, file="interspecies TK plot.pdf", width=10, height=6, scale=0.7)
