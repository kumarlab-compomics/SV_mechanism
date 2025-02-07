#!/usr/bin/Rscript

library(DNAshapeR)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

pred_SV = getShape(args[1])
pred_PRE = getShape(args[2])
pred_POST = getShape(args[3])

needtosplit = toString(args[1])
#stringerID = strsplit(strsplit(last(strsplit(needtosplit, split = "/")[[1]]), split = 'SV_')[[1]][2], split = '.fa')[[1]]
stringerID = toString(args[4])

#stringerCSV = paste0(strsplit(needtosplit, split = 'SV')[[1]][1], 'allDNA_ShapeR.csv')
stringerCSV = paste0(strsplit(needtosplit, split = '.fa'), '.csv')

print(stringerID)
print(stringerCSV)

df = data.frame(matrix(ncol = 0, nrow = 1))

df = df %>%
  mutate(
    ID = stringerID,

    SV_MGW_mean = mean(pred_SV$MGW, na.rm = TRUE),
     SV_MGW_sd = sd(pred_SV$MGW, na.rm = TRUE),
    SV_HelT_mean = mean(pred_SV$HelT, na.rm = TRUE),
     SV_HelT_sd = sd(pred_SV$HelT, na.rm = TRUE),  
      SV_ProT_mean = mean(pred_SV$ProT, na.rm = TRUE),
     SV_ProT_sd = sd(pred_SV$ProT, na.rm = TRUE),      
      SV_Roll_mean = mean(pred_SV$Roll, na.rm = TRUE),
     SV_Roll_sd = sd(pred_SV$Roll, na.rm = TRUE),    
          SV_EP_mean = mean(pred_SV$EP, na.rm = TRUE),
     SV_EP_sd = sd(pred_SV$EP, na.rm = TRUE),

    PRE_MGW_mean = mean(pred_PRE$MGW, na.rm = TRUE),
     PRE_MGW_sd = sd(pred_PRE$MGW, na.rm = TRUE),
    PRE_HelT_mean = mean(pred_PRE$HelT, na.rm = TRUE),
     PRE_HelT_sd = sd(pred_PRE$HelT, na.rm = TRUE),
      PRE_ProT_mean = mean(pred_PRE$ProT, na.rm = TRUE),
     PRE_ProT_sd = sd(pred_PRE$ProT, na.rm = TRUE),
      PRE_Roll_mean = mean(pred_PRE$Roll, na.rm = TRUE),
     PRE_Roll_sd = sd(pred_PRE$Roll, na.rm = TRUE),
      PRE_EP_mean = mean(pred_PRE$EP, na.rm = TRUE),
     PRE_EP_sd = sd(pred_PRE$EP, na.rm = TRUE),

    POST_MGW_mean = mean(pred_POST$MGW, na.rm = TRUE),
     POST_MGW_sd = sd(pred_POST$MGW, na.rm = TRUE),
    POST_HelT_mean = mean(pred_POST$HelT, na.rm = TRUE),
     POST_HelT_sd = sd(pred_POST$HelT, na.rm = TRUE),
      POST_ProT_mean = mean(pred_POST$ProT, na.rm = TRUE),
     POST_ProT_sd = sd(pred_POST$ProT, na.rm = TRUE),
      POST_Roll_mean = mean(pred_POST$Roll, na.rm = TRUE),
     POST_Roll_sd = sd(pred_POST$Roll, na.rm = TRUE),
      POST_EP_mean = mean(pred_POST$EP, na.rm = TRUE),
     POST_EP_sd = sd(pred_POST$EP, na.rm = TRUE) 
         )

#write(df, file = stringerCSV, append = TRUE, sep = "\t")

write.csv(df, file = stringerCSV, row.names = FALSE, )

