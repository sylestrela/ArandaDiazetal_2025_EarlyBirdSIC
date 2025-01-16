"Script takes report file from Skanit software and returns two one Contains the mean 620 reading for each well to be printed please specify timepoint and it will save file  For this to
work you will need to install the packages data.table and all input files must be in the same working directory as the script"
rm(list=ls())
library(data.table)

extract_data <- function(t_filename){
  raw_data = readLines(t_filename)
  rs = raw_data[(grep("Results",raw_data)[1]+2):(length(raw_data))]
  writeLines(rs,'Temp.txt')
  rs = fread('Temp.txt')
  
  file.remove('Temp.txt')
  rs = rs[,1:(ncol(rs)-1)]
  rs2 = cbind(rs[1:4,2:13],rs[5:8,2:13])
  colnames(rs2) = c('Lactose','D-Glucose','D-Galactose','Sucrose','D-Fructose','L-Rhamnose',
                    'Glycerol','D-Ribose','Pyruvate','L-Lactate','Oxaloacetate', 'Citrate',
                    'Maltose','Melibiose','D-Mannose','D-Mannitol','D-Sorbitol','Galactitol',
                    'L-Arabinose','L-Fucose','Acetate','Succinate','Fumarate','ctrl')
  rs2$Replicate = seq(1,4)
  rs2$Transfer = as.numeric(substring(strsplit(t_filename,split='_')[[1]][4][1],2))
  rs2$Inoculum = as.numeric(substring(strsplit(strsplit(t_filename,split='_')[[1]][8][1],'[.]')[[1]][1],2))
  rs2$Backbone = strsplit(t_filename,split='_')[[1]][7]
  rs2 = melt(rs2,id.vars=c('Replicate','Transfer','Inoculum','Backbone'))
  return(rs2)
}

files = list.files(path ='../Data/Raw_OD_Data/',pattern='.txt',full.names=TRUE)
dat = data.table()
for(i in files){
  dat = rbind(dat,extract_data(i))
}

fwrite(dat,'../Data/Final_OD_data.csv')