#This file reads in ALL the %-signal-change values, per-participant, per-parcel, per-contrast,
# Those %-signal-change calculations are produced by the awesome toolbox analyses, and represent a single overall calculation
#derived for the whole parcel region (not individual voxels, as mk sometimes forgets)

rm(list = ls())
library(bootstrap)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

##SET YOUR DIRECTORY!!!

########
#READ IN DATA
########
#Here, we read in all those files, calculate a whole passle of mean and standard error bars, and then make graphs

# Add in the contrast and ROI names so it's not just numbers!!!!! (This ordering comes from the 
# standard ordering produced by the 2nd level analyses; we'll arrange differently in the plots)

RHLangROI.Names = c('RPost Temp', 'RAnt Temp', 'RAngG', 'RIFG',      'RMFG',     'RIFG orb');
LangROI.Names = c('LPost Temp', 'LAnt Temp', 'LAngG', 'LIFG',      'LMFG',     'LIFG orb');

MDROI.Names = c('LIFG op',  'RIFG op', 'LMFG',    'RMFG',    'LMFG orb',
                      'RMFG orb', 'LPrecG', 'RPrecG',  'LInsula', 'RInsula',
                      'LSMA',    'RSMA',   'LPar Inf', 'RPar Inf', 'LPar Sup',
                      'RPar Sup', 'LACC',   'RACC');

ToMROI.Names = c('DM PFC', 'LTPJ',  'MM PFC', 'PC',
                      'RTPJ',  'VM PFC', 'RSTS');

normal.contrasts = c('agt', 'pat')


# myResults = read.csv('RHLangfROIsrespNonlitJokes.csv')%>%
#   mutate(ROIName = RHLangROI.Names[ROI]) %>%
#   mutate(contrastName = normal.contrasts[Contrast])%>%
#   mutate(Group = 'RHLang')
# allSigChange = myResults

myResults = read.csv('toolbox_outputs/langloc_resp_agpat_20170309.csv') %>%
  mutate(ROIName = LangROI.Names[ROI]) %>%
  mutate(contrastName = normal.contrasts[Contrast])%>%
  mutate(Group = 'LHLang')
#allSigChange = rbind(allSigChange, myResults)
allSigChange = myResults

myResults = read.csv('toolbox_outputs/MDloc_resp_agpat_20170317.csv') %>%
  mutate(ROIName = MDROI.Names[ROI]) %>%
  mutate(contrastName = normal.contrasts[Contrast]) %>%
  mutate(Group = 'MDAll')
allSigChange = rbind(allSigChange, myResults)

#Little extra thing here, rename MD to split by L and R hemisphere!
allSigChange[(allSigChange$Group == 'MDAll') & (allSigChange$ROI %%2 == 1),]$Group = 'MDLeft'
allSigChange[(allSigChange$Group == 'MDAll') & (allSigChange$ROI %%2 == 0),]$Group = 'MDRight'

myResults = read.csv('toolbox_outputs/ToMloc_resp_agpat_20170317.csv')%>%
  mutate(ROIName = ToMROI.Names[ROI]) %>%
  mutate(contrastName = normal.contrasts[Contrast]) %>%
  mutate(Group = 'ToM')
allSigChange = rbind(allSigChange, myResults)

#########
# TRANSFORMATIONS
#########

#First, in addition to the by-region signal changes, we are going to give each person an average signal change value for each localizer 
avgSigChange = aggregate(allSigChange$sigChange, by=list(allSigChange$Group,allSigChange$SubjectNumber,allSigChange$contrastName), mean)
names(avgSigChange) = c('Group','SubjectNumber', 'contrastName','sigChange')
avgSigChange$ROIName = 'LocalizerAverage'
avgSigChange$ROI = 0

allSigChange <- allSigChange %>%
  dplyr::select(one_of(c('Group','ROIName', 'ROI','SubjectNumber', 'contrastName','sigChange')))

allSigChange <- rbind(allSigChange, avgSigChange)

#Drop the contrasts we're not interested in...
toGraph = allSigChange %>%
  filter(contrastName %in% c('agt','pat'))

#Next, get the table that we'll be making the graphs from: for each region (including the average region), take all 
#the individual signal changes and calculate a mean, a standard error (incase we want it) and bootstrapped CIs
sterr <- function(mylist){
  my_se = sd(mylist)/sqrt(length(mylist)) 
  
  return(my_se)
}

mystats = aggregate(toGraph$sigChange, by=list(toGraph$Group, toGraph$ROIName, toGraph$ROI,toGraph$contrastName), mean)
names(mystats) = c('Group','ROIName', 'ROI','contrastName', 'themean')
myster = aggregate(toGraph$sigChange, by=list(toGraph$Group, toGraph$ROIName, toGraph$ROI,toGraph$contrastName), sterr)
names(myster) = c('Group','ROIName', 'ROI','contrastName', 'sterr')

mystats = merge(mystats,myster)
mystats$se_up = mystats$themean + mystats$sterr
mystats$se_down = mystats$themean - mystats$sterr

#Print out a simple summary for mega-graphs
avgz <- filter(mystats, ROIName == 'LocalizerAverage')
write.csv(avgz, 'jokes_localizer_avg.csv')

#bootstrapped 95% confidence intervals! calculate them from allSigChange
#then merge into mystats
bootup <- function(mylist){
  foo <- bootstrap(mylist, 1000, mean)
  return(quantile(foo$thetastar, 0.975)[1])
}
bootdown <- function(mylist){
  foo <- bootstrap(mylist, 1000, mean)
  return(quantile(foo$thetastar, 0.025)[1])
}
mybootup = aggregate(toGraph$sigChange, by=list(toGraph$Group, toGraph$ROIName, toGraph$ROI, toGraph$contrastName), bootup)
names(mybootup) = c('Group','ROIName', 'ROI','contrastName', 'bootup')
mybootdown = aggregate(toGraph$sigChange, by=list(toGraph$Group, toGraph$ROIName, toGraph$ROI, toGraph$contrastName), bootdown)
names(mybootdown) = c('Group','ROIName', 'ROI','contrastName', 'bootdown')

mystats = merge(mystats,mybootup)
mystats = merge(mystats,mybootdown)

#########
# Effect size reports
#########
#report  a simple measure of effect size: the
#mean signal change in each system. Here they are:
eff <- mystats %>%
  filter(ROIName == 'LocalizerAverage') %>%
  filter(contrastName == 'agt' | contrastName == 'pat') %>%
  dplyr::select(Group, contrastName,themean) %>%
  spread(contrastName, themean) %>%
  mutate(sigChange = agt-pat)
  
            
#########
# Graphs!
#########

#Now we can use the information stored in mystats to make pretty graphs! This could be done in excel too by printing mystats
#Change to figs output folder
#setwd("~/Dropbox/_Projects/Jokes - fMRI/Jokes-Analysis Repository/Analyses_paper/reproducible analyses/figs")


#Select the rows we want for each graph, and order them how we want! For now, localizerAverage will just come first in all sets
mystats$contNo <- 1
mystats[mystats$contrastName == 'agt',]$contNo <- 1
mystats[mystats$contrastName == 'pat',]$contNo <- 2
mystats = arrange(mystats, ROI)
mystats = arrange(mystats, contNo)

#Add a new col grouping to separate out the localizer average
mystats$ROIGroup <- ""
mystats[mystats$ROIName == "LocalizerAverage",]$ROIGroup <- "across fROIs"
mystats = arrange(mystats, desc(ROIGroup))

#Changes for prettiness
mystats[mystats$ROIName=="LocalizerAverage",]$ROIName <- "average across fROIs"
mystats$ROIName <- str_wrap(mystats$ROIName, width = 4)

mystats$contrastLabel <- mystats$contrastName
mystats[mystats$contrastName == "agt",]$contrastLabel <- "Agent highlight\n  "
mystats[mystats$contrastName == "pat",]$contrastLabel <- "Patient highlight\n   "


#Subsets & Ordering (elaborate code, probably can condense these; ggplot is finicky at orders)
# RHLang = filter(mystats, Group == 'RHLang')
# RHLang <- RHLang[order(RHLang$ROI),]
# RHLang$PresOrder = c(13,14, 9,10, 7,8, 11,12, 3,4,5,6,1,2) #Reorder for standard presentation!
# RHLang <- RHLang[order(RHLang$PresOrder),]
# RHLang = arrange(RHLang, desc(ROIGroup))


LHLang = filter(mystats, Group == 'LHLang')
LHLang <- LHLang[order(LHLang$ROI),]
LHLang$PresOrder = c(13,14, 9,10, 7,8, 11,12, 3,4,5,6,1,2)
LHLang <- LHLang[order(LHLang$PresOrder),]
LHLang = arrange(LHLang, desc(ROIGroup))


MDLeft = filter(mystats, Group == 'MDLeft')
MDLeft <- MDLeft[order(MDLeft$ROI),]
MDLeft = arrange(MDLeft, desc(ROIGroup))

MDRight = filter(mystats, Group == 'MDRight')
MDRight <- MDRight[order(MDRight$ROI),]
MDRight = arrange(MDRight, desc(ROIGroup))

ToM = filter(mystats, Group == 'ToM')
ToM <- ToM[order(ToM$ROI),]
ToM$PresOrder = c(1,2,3,4,9,10,5,6,7,8,11,12,13,14,15,16)
ToM <- ToM[order(ToM$PresOrder),]
ToM = arrange(ToM, desc(ROIGroup))

#Graphing function!

makeBar = function(plotData,ylow=-0.5,yhigh=2.5, mycolors = c("gray35", "gray60")) {

  #freeze factor orders
  plotData$ROIName <- factor(plotData$ROIName, levels = unique(plotData$ROIName))
  plotData$ROIGroup <- factor(plotData$ROIGroup, levels = unique(plotData$ROIGroup))
  plotData$contrastLabel <- factor(plotData$contrastLabel, levels = unique(plotData$contrastLabel))
  myfi = paste(plotData$Group[1], '.jpg', sep="")#filename
  print(myfi)

ggplot(data=plotData, aes(x=ROIName, y=themean, fill=contrastLabel)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=bootdown, ymax=bootup), colour="black", width=.1, position=position_dodge(.9)) +
  coord_cartesian(ylim=c(ylow,yhigh)) +
  scale_y_continuous(breaks = seq(-0.5, 2.5, 0.5))+
  xlab('') +
  ylab(str_wrap('% signal change over fixation', width=18)) +
  scale_fill_manual(name="", values=mycolors) +
  theme_bw() +
  theme(legend.key = element_blank()) +
  #theme(text = element_text(size = 40)) +
  theme(text = element_text(size = 25)) +
  facet_grid(~ROIGroup, scale='free_x', space='free_x') +
  theme(strip.background = element_blank()) +
  theme(strip.text = element_blank()) 
  # Optional, remove for RHLang and ToMCustom since we want the legend there...
  #+ theme(legend.position="none")
 

  ggsave(filename=myfi, width=length(unique(plotData$ROIName))*2.2, height=6.1)
  
}

makeBar(LHLang)
makeBar(MDLeft)
makeBar(MDRight)
makeBar(ToM, -0.5, 1)
