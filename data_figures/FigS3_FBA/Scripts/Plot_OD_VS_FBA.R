rm(list=ls())
library(ggplot2)
library(ggpmisc)
library(data.table)
library(viridis)


a0= fread('../Data/Final_OD_data.csv')
#Lactate precipitated so we removed it from the data
a0 = a0[a0$variable != "L-Lactate"]
b0 = fread('../Data/FBA_Growth_Rates.csv')

b0$variable = b0$Carbon_Source
b0[is.na(b0$Growth_Rate_E_coli_Anaerobic)]$Growth_Rate_E_coli_Anaerobic =0
for(j in unique(a0$Transfer)){
  a = a0[Transfer==j]
  b = merge(a,b0,by='variable')
  formula <- y~x
  p1 <- ggplot(b,aes(x=Pred_Growth_Rate,y=value,col=Backbone)) + 
    labs(x = 'FBA Growth Rate (Ecoli) /Cmoles', y ='OD') +
    facet_wrap(~paste('Inoculum', Inoculum,sep =' '),ncol=4 ) + 
    geom_point() +theme_classic()+
    geom_smooth(method='lm',formula=formula,se =F) +
    stat_poly_eq(aes(label = paste(..rr.label..)), 
                 formula = formula, parse = TRUE, size = 3) +
    scale_y_continuous(breaks=c(0,0.6)) + 
    scale_x_continuous(breaks=c(0.01,0.04))
  ggsave(paste('../Plots/FBAGrowthVS_OOD',j,'.png',sep=''),p1)
}

final = fread('../Data/Final_OD_data.csv')[Transfer==7]
ctrl = final[final$Backbon=='M9BHI' & final$variable=='ctrl',]
b0 = fread('../Data/FBA_Growth_Rates.csv')
final =final[final$variable %in% b0$Carbon_Source & final$Backbone=='M9BHI',]

final$Predicted =b0$Pred_Growth_Rate[match(final$variable,b0$Carbon_Source)]
final[is.na(final$Predicted),]$Predicted =0
#Remove lactate because of
final= final[final$variable!= 'L-Lactate',]
final$Observed = final$value
final = final[,list(mean=mean(value)),by=list(Backbone,variable,Inoculum,Predicted)]
ctrl = ctrl[,list(mean=mean(value)),by=list(Inoculum)]

final$Ctrl =ctrl$mean[match(final$Inoculum,ctrl$Inoculum)]
final$Corrected =  final$mean-final$Ctrl
final[final$Corrected<0]$Corrected = 0
p1 <-ggplot(final,mapping = aes(x=Predicted,y=Corrected,col=as.factor(Inoculum))) +
  geom_point(size=2,shape=1,stroke=1) +theme_classic() +
  labs(x='FBA Predicted Yield',y = expression(Delta*OD),col='Donor')+     
  geom_smooth(method='lm',formula=formula,se =F,linetype=2) +
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "left",
               formula = formula, parse = TRUE, size = 4) +scale_color_brewer(palette = "Dark2") +scale_y_continuous(breaks=c(0,0.5)) +
  scale_x_continuous(breaks=c(0,0.005),limits=c(0,0.005))+
  scale_y_continuous(limits=c(0,0.5))
  theme(axis.title=element_text(size=12),axis.text=element_text(size=10),axis.line = element_line(size=1))
ggsave('../Plots/FigureS3.png',p1,height=4,width=5)