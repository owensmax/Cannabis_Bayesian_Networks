#BUILD FINAL DATASET

#load libraries for analyses
library("FactoMineR") #used for principal components anlayses
library("psych") #used for general stats
library("bnlearn") #used for Bayesian network analyses
library("lmerTest") #used for residualizing
source("/home/max/Documents/Bayes/correlation_matrix.R") #used for format Table 1

#reset plotting in R
par(mar=c(1,1,1,1))

#read in data
ROIs <-read.csv("/home/max/Documents/Bayes/UPDATED_SPCs_WITH_BSL_AND_FUP_THICK.csv")
BIS_fup <-read.csv("/home/max/Documents/Bayes/FU2_799_BIS.csv")
names(BIS_fup)[19:21]=c('BIS_SUM_Attention_FU2','BIS_SUM_Motor_FU2','BIS_SUM_Planning_FU2')
candat <-read.csv("/home/max/Documents/Bayes/cannabis_data.csv")
nic <-read.csv("/home/max/Documents/Bayes/lifetime_nic.csv")
audit <-read.csv("/home/max/Documents/Bayes/AUDIT_C.csv")
prs<-read.csv("/home/max/Documents/Bayes/IMAGEN_Cann-PRS_MDS.csv") #WRONG S1-S8, THOUGH C1-C4 IS OK
prs3<-read.csv("/home/max/Documents/Bayes/PRS_CANNnew.csv")
prs4<-read.csv("/home/max/Documents/Bayes/PRS-CUDnew.csv")

#rename columns and merge datasets
colnames(prs)[1]<-'Subject'
colnames(prs3)[1]<-'Subject'
colnames(prs4)[1]<-'Subject'
candat<-cbind(candat,nic[c('Life_Nicotine')],audit[c('AUDIT_C')])
candat_noprs <- candat
candat<-merge(candat,prs[c('Subject','C1','C2','C3','C4')],by="Subject")
candat<-merge(candat,prs3,by="Subject")
candat<-merge(candat,prs4,by="Subject")
candat_bsl <-candat[ which(candat$Time=='BSL'),]
candat_fup <-candat[ which(candat$Time=='FUP'),]
candat2<-merge(candat_bsl,candat_fup,by='Subject')

#rename columns to improve clarity
colnames(candat2) <- gsub('.x','_BSL',colnames(candat2))
colnames(candat2) <- gsub(".y", "_FU2",colnames(candat2))
colnames(candat2) <- gsub('S_BSL_BSL','Sex',colnames(candat2))
colnames(candat2) <- gsub('S_BSL_FU2','Sex_FU2',colnames(candat2))
colnames(candat2) <- gsub('Mod_to_Hea_FU2_at_FU2_FU2','Mod_to_Heavy_at_FU2',colnames(candat2))
colnames(candat2) <- gsub('Mod_to_Hea_FU2_at_FU2_BSL','Mod_to_Heavy_at_BSL',colnames(candat2))
colnames(candat2) <- gsub('Sex','Sex',colnames(candat2))
colnames(candat2) <- gsub('Hand_BSL','Handedness',colnames(candat2))
colnames(candat2) <- gsub('Age_BSL','Baseline_Age',colnames(candat2))
colnames(candat2) <- gsub('Site_BSL','Site',colnames(candat2))

#standardize polygenic risk scores
candat2$C1_FU2<-scale(candat2$C1_FU2)
candat2$C2_FU2<-scale(candat2$C2_FU2)
candat2$C3_FU2<-scale(candat2$C3_FU2)
candat2$C4_FU2<-scale(candat2$C4_FU2)
candat2$S1CANNnew_FU2<-scale(candat2$S1CANNnew_FU2)
candat2$S2CANNnew_FU2<-scale(candat2$S2CANNnew_FU2)
candat2$S3CANNnew_FU2<-scale(candat2$S3CANNnew_FU2)
candat2$S4CANNnew_FU2<-scale(candat2$S4CANNnew_FU2)
candat2$S5CANNnew_FU2<-scale(candat2$S5CANNnew_FU2)
candat2$S6CANNnew_FU2<-scale(candat2$S6CANNnew_FU2)
candat2$S7CANNnew_FU2<-scale(candat2$S7CANNnew_FU2)
candat2$S8CANNnew_FU2<-scale(candat2$S8CANNnew_FU2)
candat2$S1CUDnew_FU2<-scale(candat2$S1CUDnew_FU2)
candat2$S2CUDnew_FU2<-scale(candat2$S2CUDnew_FU2)
candat2$S3CUDnew_FU2<-scale(candat2$S3CUDnew_FU2)
candat2$S4CUDnew_FU2<-scale(candat2$S4CUDnew_FU2)
candat2$S5CUDnew_FU2<-scale(candat2$S5CUDnew_FU2)
candat2$S6CUDnew_FU2<-scale(candat2$S6CUDnew_FU2)
candat2$S7CUDnew_FU2<-scale(candat2$S7CUDnew_FU2)
candat2$S8CUDnew_FU2<-scale(candat2$S8CUDnew_FU2)

#make reduced dataset
imagen_data=Reduce(function(x,y) merge(x = x, y = y, by = "Subject"), 
       list(ROIs, candat2[c('Subject','Site','Sex','Handedness','Baseline_Age','Life_BSL','Mod_to_Heavy_at_BSL','Site_FU2','Sex_FU2','Hand_FU2','Age_FU2','Life_FU2','Mod_to_Heavy_at_FU2','Life_Nicotine_BSL','Life_Nicotine_FU2','AUDIT_C_FU2','AUDIT_C_BSL','S5CANNnew_FU2','S5CUDnew_FU2','C1_FU2','C2_FU2','C3_FU2','C4_FU2')]))
imagen_data<-imagen_data[complete.cases(imagen_data), ]

#create change scores
imagen_data$rdmpfc_change<-imagen_data$right_dmPFC_FUP-imagen_data$right_dmPFC_BSL
imagen_data$ldmpfc_change<-imagen_data$left_dmPFC_FUP-imagen_data$left_dmPFC_BSL
imagen_data$Age_change<-imagen_data$Age_FU2-imagen_data$Baseline_Age
imagen_data$Life_change<-imagen_data$Life_FU2-imagen_data$Life_BSL
imagen_data$AUDIT_change<-imagen_data$AUDIT_C_FU2-imagen_data$AUDIT_C_BSL
imagen_data$Life_Nicotine_change<-imagen_data$Life_Nicotine_FU2-imagen_data$Life_Nicotine_BSL

#MAKE PCA SCORES OF BRAIN VARS
brain.pca.bl <- PCA(imagen_data[,c("right_dmPFC_BSL","left_dmPFC_BSL")], scale.unit=TRUE, graph=T)
imagen_data$brain.pca.bl<-brain.pca.bl$ind$coord[,1]
brain.pca.fu2 <- PCA(imagen_data[,c("right_dmPFC_FUP","left_dmPFC_FUP")], scale.unit=TRUE, graph=T)
imagen_data$brain.pca.fu2<-brain.pca.fu2$ind$coord[,1]
brain.pca.change <- PCA(imagen_data[,c("rdmpfc_change","ldmpfc_change")], scale.unit=TRUE, graph=T)
imagen_data$brain.pca.change<-brain.pca.change$ind$coord[,1]
brain.pca.spc <- PCA(imagen_data[,c("spc_left_dmPFC","spc_right_dmPFC")], scale.unit=TRUE, graph=T)
imagen_data$brain.pca.spc<-brain.pca.spc$ind$coord[,1]

#TURN VARIABLES OF INTEREST INTO RESIDUALS
#brain interecept (BL)
model1<-lmer(brain.pca.bl~ (1|Site) + (1|Sex), data=imagen_data)
Baseline_dPFC_Thickness <- as.matrix(residuals(model1))

#brain change resid
model1<-lmer(brain.pca.change~ (1|Site) + (1|Sex), data=imagen_data)
Change_in_dPFC_Thickness <- as.matrix(residuals(model1))

#canabis growth
model1<-lmer(Life_change~ (1|Site) + (1|Sex), data=imagen_data)
Cannabis_Use_by_Age_19<- as.matrix(residuals(model1))

#alcohol intercept
model1<-lmer(AUDIT_C_BSL~ (1|Site) + (1|Sex), data=imagen_data)
Baseline_Alcohol_Use<- as.matrix(residuals(model1))

#alcohol growth
model1<-lmer(AUDIT_change~ (1|Site) + (1|Sex), data=imagen_data)
Change_in_Alcohol_Use<- as.matrix(residuals(model1))

#nicotine intercept
model1<-lmer(Life_Nicotine_BSL~ (1|Site) + (1|Sex), data=imagen_data)
Baseline_Tobacco_Use<- as.matrix(residuals(model1))

#nicotine growth
model1<-lmer(Life_Nicotine_change~ (1|Site) + (1|Sex), data=imagen_data)
Change_in_Tobacco_Use<- as.matrix(residuals(model1))

#age
model1<-lmer(Baseline_Age~ (1|Site) + (1|Sex), data=imagen_data)
Baseline_Age<- as.matrix(residuals(model1))

#handedness
model1<-lmer(Handedness~(1|Site) + (1|Sex), data=imagen_data)
Handedness<- as.matrix(residuals(model1))

#Cannabis Use Lifetime PRS
model1<-lmer(S5CANNnew_FU2~ C1_FU2+C2_FU2+C3_FU2+C4_FU2+(1|Site), data=imagen_data)
S5CANNnew_FU2<- as.matrix(residuals(model1))

#create corr table with significance values (Table 1)
write.csv(correlation_matrix(imagen_data[c('Handedness','Baseline_Age','Life_Nicotine_BSL','AUDIT_C_BSL','S5CANNnew_FU2','Life_change','AUDIT_change','Life_Nicotine_change','brain.pca.bl','brain.pca.change')], 
           type = "pearson",
           digits = 2, 
           decimal.mark = ".",
           use = "lower", 
           show_significance = TRUE, 
           replace_diagonal = TRUE, 
           replacement = ""), "/home/max/Documents/Bayes/bayes_cor_table_sigs.csv")

#make final BN network dataset
data.regressed <- cbind(Baseline_dPFC_Thickness,Change_in_dPFC_Thickness,Cannabis_Use_by_Age_19,Baseline_Tobacco_Use,Change_in_Tobacco_Use,Baseline_Alcohol_Use,Change_in_Alcohol_Use,Handedness,Baseline_Age,S5CANNnew_FU2)## assembles the residuals

#make final dataset into a dataframe
data.regressed <- as.data.frame(scale(data.regressed, center = T, scale = T))## center+scale the data 
colnames(data.regressed) <- c('Baseline_dPFC_Thickness','Change_in_dPFC_Thickness','Cannabis_Use_by_Age_19','Baseline_Tobacco_Use','Change_in_Tobacco_Use','Baseline_Alcohol_Use','Change_in_Alcohol_Use','Handedness','Baseline_Age','Cannabis_Use_PRS')## assigns variable names, in the same order as above

## generate the blacklist of the directions that you are not interested in
blacklist = data.frame(
  from =c(
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Baseline_dPFC_Thickness',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Baseline_dPFC_Thickness',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Baseline_dPFC_Thickness',
    'Baseline_Alcohol_Use',
    'Baseline_Tobacco_Use',
    'Handedness',
    'Baseline_Age'
    ), 
to = c(
  'Baseline_dPFC_Thickness',
  'Baseline_dPFC_Thickness',
  'Baseline_dPFC_Thickness',
  'Baseline_dPFC_Thickness',
  'Baseline_Alcohol_Use',
  'Baseline_Alcohol_Use',
  'Baseline_Alcohol_Use',
  'Baseline_Alcohol_Use',
  'Baseline_Tobacco_Use',
  'Baseline_Tobacco_Use',
  'Baseline_Tobacco_Use',
  'Baseline_Tobacco_Use',
  'Handedness',
  'Handedness',
  'Handedness',
  'Handedness',
  'Handedness',
  'Baseline_Age',
  'Baseline_Age',
  'Baseline_Age',
  'Baseline_Age',
  'Baseline_Age',
  'Cannabis_Use_PRS',
  'Cannabis_Use_PRS',
  'Cannabis_Use_PRS',
  'Cannabis_Use_PRS',
  'Cannabis_Use_PRS',
  'Cannabis_Use_PRS',
  'Cannabis_Use_PRS',
  'Cannabis_Use_PRS',
  'Cannabis_Use_PRS'
)) 

## generate the blacklist that includes cortical thickness to cannabis use on blacklist
blacklist_Brain2Cannabis = data.frame(
  from =c(
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Baseline_dPFC_Thickness',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Baseline_dPFC_Thickness',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Baseline_dPFC_Thickness',
    'Baseline_Alcohol_Use',
    'Baseline_Tobacco_Use',
    'Handedness',
    'Baseline_Age',
    'Change_in_dPFC_Thickness'
  ), 
  to = c(
    'Baseline_dPFC_Thickness',
    'Baseline_dPFC_Thickness',
    'Baseline_dPFC_Thickness',
    'Baseline_dPFC_Thickness',
    'Baseline_Alcohol_Use',
    'Baseline_Alcohol_Use',
    'Baseline_Alcohol_Use',
    'Baseline_Alcohol_Use',
    'Baseline_Tobacco_Use',
    'Baseline_Tobacco_Use',
    'Baseline_Tobacco_Use',
    'Baseline_Tobacco_Use',
    'Handedness',
    'Handedness',
    'Handedness',
    'Handedness',
    'Handedness',
    'Baseline_Age',
    'Baseline_Age',
    'Baseline_Age',
    'Baseline_Age',
    'Baseline_Age',
    'Cannabis_Use_PRS',
    'Cannabis_Use_PRS',
    'Cannabis_Use_PRS',
    'Cannabis_Use_PRS',
    'Cannabis_Use_PRS',
    'Cannabis_Use_PRS',
    'Cannabis_Use_PRS',
    'Cannabis_Use_PRS',
    'Cannabis_Use_PRS',
    'Cannabis_Use_by_Age_19'
  )) 

##  generate the blacklist that includes cannabis use to cortical thickness on blacklist
blacklist_Cannabis2Brain = data.frame(
  from =c(
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Baseline_dPFC_Thickness',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Baseline_dPFC_Thickness',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Baseline_dPFC_Thickness',
    'Baseline_Alcohol_Use',
    'Baseline_Tobacco_Use',
    'Handedness',
    'Baseline_Age',
    'Cannabis_Use_by_Age_19'
  ), 
  to = c(
    'Baseline_dPFC_Thickness',
    'Baseline_dPFC_Thickness',
    'Baseline_dPFC_Thickness',
    'Baseline_dPFC_Thickness',
    'Baseline_Alcohol_Use',
    'Baseline_Alcohol_Use',
    'Baseline_Alcohol_Use',
    'Baseline_Alcohol_Use',
    'Baseline_Tobacco_Use',
    'Baseline_Tobacco_Use',
    'Baseline_Tobacco_Use',
    'Baseline_Tobacco_Use',
    'Handedness',
    'Handedness',
    'Handedness',
    'Handedness',
    'Handedness',
    'Baseline_Age',
    'Baseline_Age',
    'Baseline_Age',
    'Baseline_Age',
    'Baseline_Age',
    'Cannabis_Use_PRS',
    'Cannabis_Use_PRS',
    'Cannabis_Use_PRS',
    'Cannabis_Use_PRS',
    'Cannabis_Use_PRS',
    'Cannabis_Use_PRS',
    'Cannabis_Use_PRS',
    'Cannabis_Use_PRS',
    'Cannabis_Use_PRS',
    'Change_in_dPFC_Thickness'
  )) 

### the following is using different algorithms provided by the bnlearn package to get the BCN, 

# bootstrap hill climbing

#Use BIC for scoring
boot = boot.strength(data = data.regressed, R = 10000, algorithm = "hc", 
                     algorithm.args = list(start=NULL,score ="bic-g", blacklist = blacklist))

boot[(boot$strength > 0.5) & (boot$direction >= 0.5), ]

avg.boot_hc = averaged.network(boot, threshold = 0.80)
plot(avg.boot_hc,main="Hill Climbing")

##aic instead of bic
boot = boot.strength(data = data.regressed, R = 10000, algorithm = "hc", 
                     algorithm.args = list(start=NULL,score ="aic-g", blacklist = blacklist))

boot[(boot$strength > 0.8) & (boot$direction >= 0.5), ]

#set up random starting places for the hill climb
random_graphs_starts=function(num){
  start=random.graph(nodes = names(data.regressed), num = num)
  return(start)
}

bn_multi_start_str=function(start){
  boot=boot.strength(data = data.regressed, R = 100, algorithm = "hc", 
                     algorithm.args = list(start=start,score ="bic-g",  
                                           blacklist = blacklist))
  return(boot[24,3])
}


bn_multi_start_dir=function(start){
  boot=boot.strength(data = data.regressed, R = 100, algorithm = "hc", 
                     algorithm.args = list(start=start,score ="bic-g",  
                                           blacklist = blacklist))
  return(boot[24,4])
}

start=random_graphs_starts(100)

direction=sapply(start,bn_multi_start_dir)
strength=sapply(start,bn_multi_start_str)

mean(direction)
mean(strength)

#run hill climb with random perturbations introducted
boot = boot.strength(data = data.regressed, R = 10000, algorithm = "hc", 
                     algorithm.args = list(start=NULL,score ="bic-g", blacklist = blacklist,perturb=10))

boot[(boot$strength > 0.8) & (boot$direction >= 0.5), ]

################test with one direction blacklisted and see which does better#######
#block brain to cannabis directed arc
boot = boot.strength(data = data.regressed, R = 10000, algorithm = "hc", 
                     algorithm.args = list(start=NULL,score ="bic-g", blacklist = blacklist_Cannabis2Brain))
boot[(boot$strength > 0.9) & (boot$direction >= 0.5), ]

bn.hc.1 <- hc(data.regressed, blacklist = blacklist_Cannabis2Brain)

score(bn.hc.1,data.regressed, type='bic-g') 

#block cannabis to brain directed arc
boot = boot.strength(data = data.regressed, R = 10000, algorithm = "hc", 
                     algorithm.args = list(start=NULL,score ="bic-g", blacklist = blacklist_Brain2Cannabis))
boot[(boot$strength > 0.9) & (boot$direction >= 0.5), ]

bn.hc.2 <- hc(data.regressed, blacklist = blacklist_Brain2Cannabis)
plot(bn.hc)
bn.hc

score(bn.hc.2,data.regressed, type='bic-g')

###############tabu (a score based method based on hill climb)###########

#test some new scoring methods
boot = boot.strength(data = data.regressed, R = 10000, algorithm = "tabu", 
                     algorithm.args = list(start=NULL,score ="bic-g", blacklist = blacklist))

boot[(boot$strength > 0.9) & (boot$direction >= 0.5), ]

avg.boot_hc = averaged.network(boot, threshold = 0.90)
plot(avg.boot_hc,main="Tabu")

########### Grow shrink###########

#primary grow shrink
boot = boot.strength(data = data.regressed, R = 10000, algorithm = "gs", 
                     algorithm.args = list(blacklist = blacklist, test= 'mi-g'))

boot[(boot$strength > 0.8) & (boot$direction >= 0.5), ]

avg.boot_gs = averaged.network(boot, threshold = 0.90)
plot(avg.boot_gs,main="Grow Shrink")


###GS with different independence tests
boot = boot.strength(data = data.regressed, R = 10000, algorithm = "gs", debug = FALSE,
                     algorithm.args = list(blacklist = blacklist, test= 'mc-zf', alpha = 0.05, B=10))

boot[(boot$strength > 0.8) & (boot$direction >= 0.7), ]

boot = boot.strength(data = data.regressed, R = 10000, algorithm = "gs", debug = FALSE,
                     algorithm.args = list(blacklist = blacklist, test= 'cor', alpha = 0.05, B=10))

boot[(boot$strength > 0.8) & (boot$direction >= 0.7), ]

avg.boot_gs = averaged.network(boot, threshold = 0.90)
plot(avg.boot_gs,main="Grow Shrink")

# bootstrap iamb
# ver 1: 100000 and 0.90

boot = boot.strength(data = data.regressed, R = 10000, algorithm = "iamb", 
                     algorithm.args = list(blacklist = blacklist))

boot[(boot$strength > 0.8) & (boot$direction >= 0.5), ]

avg.boot_iamb = averaged.network(boot, threshold = 0.90)
plot(avg.boot_iamb,main="IAMB")

################ min max hill climb (hybrid algorithm)#######
# ver 1: 100000 and 0.90

bn.mmhc <- mmhc(data.regressed, blacklist = blacklist)
plot(bn.mmhc)

boot = boot.strength(data = data.regressed, R = 10000, algorithm = "mmhc", 
                     algorithm.args = list(blacklist = blacklist))

boot[(boot$strength > 0.7) & (boot$direction >= 0.5), ]

###############rsmax2 (hybrid)###########

#test some new scoring methods
boot = boot.strength(data = data.regressed, R = 10000, algorithm = "rsmax2", 
                     algorithm.args = list(blacklist = blacklist))

boot[(boot$strength > 0.7) & (boot$direction >= 0.5), ]
avg.boot_hc = averaged.network(boot, threshold = 0.90)
plot(avg.boot_hc,main="Max Min Hill Climbing")


##########################make version without prs for comparison###########################################
candat_bsl_noprs <-candat_noprs[ which(candat_noprs$Time=='BSL'),]
candat_fup_noprs <-candat_noprs[ which(candat_noprs$Time=='FUP'),]
candat2_noprs<-merge(candat_bsl_noprs,candat_fup_noprs,by='Subject')

colnames(candat2_noprs) <- gsub('.x','_BSL',colnames(candat2_noprs))
colnames(candat2_noprs) <- gsub(".y", "_FU2",colnames(candat2_noprs))
colnames(candat2_noprs) <- gsub('S_BSL_BSL','Sex',colnames(candat2_noprs))
colnames(candat2_noprs) <- gsub('S_BSL_FU2','Sex_FU2',colnames(candat2_noprs))
colnames(candat2_noprs) <- gsub('Mod_to_Hea_FU2_at_FU2_FU2','Mod_to_Heavy_at_FU2',colnames(candat2_noprs))
colnames(candat2_noprs) <- gsub('Mod_to_Hea_FU2_at_FU2_BSL','Mod_to_Heavy_at_BSL',colnames(candat2_noprs))
colnames(candat2_noprs) <- gsub('Sex','Sex',colnames(candat2_noprs))
colnames(candat2_noprs) <- gsub('Hand_BSL','Handedness',colnames(candat2_noprs))
colnames(candat2_noprs) <- gsub('Age_BSL','Baseline_Age',colnames(candat2_noprs))
colnames(candat2_noprs) <- gsub('Site_BSL','Site',colnames(candat2_noprs))

imagen_data_noprs=Reduce(function(x,y) merge(x = x, y = y, by = "Subject"), 
                         list(ROIs, candat2_noprs[c('Subject','Site','Sex','Handedness','Baseline_Age','Life_BSL','Mod_to_Heavy_at_BSL','Site_FU2','Sex_FU2','Hand_FU2','Age_FU2','Life_FU2','Mod_to_Heavy_at_FU2','Life_Nicotine_BSL','Life_Nicotine_FU2','AUDIT_C_FU2','AUDIT_C_BSL')]))
imagen_data_noprs<-imagen_data_noprs[complete.cases(imagen_data_noprs), ]

imagen_data_noprs$rdmpfc_change<-imagen_data_noprs$right_dmPFC_FUP-imagen_data_noprs$right_dmPFC_BSL
imagen_data_noprs$ldmpfc_change<-imagen_data_noprs$left_dmPFC_FUP-imagen_data_noprs$left_dmPFC_BSL
imagen_data_noprs$Age_change<-imagen_data_noprs$Age_FU2-imagen_data_noprs$Baseline_Age
imagen_data_noprs$Life_change<-imagen_data_noprs$Life_FU2-imagen_data_noprs$Life_BSL
imagen_data_noprs$AUDIT_change<-imagen_data_noprs$AUDIT_C_FU2-imagen_data_noprs$AUDIT_C_BSL
imagen_data_noprs$Life_Nicotine_change<-imagen_data_noprs$Life_Nicotine_FU2-imagen_data_noprs$Life_Nicotine_BSL

imagen_data_noprs$rdmpfc_change<-imagen_data_noprs$right_dmPFC_FUP-imagen_data_noprs$right_dmPFC_BSL
imagen_data_noprs$ldmpfc_change<-imagen_data_noprs$left_dmPFC_FUP-imagen_data_noprs$left_dmPFC_BSL
imagen_data_noprs$Age_change<-imagen_data_noprs$Age_FU2-imagen_data_noprs$Baseline_Age
imagen_data_noprs$Life_change<-imagen_data_noprs$Life_FU2-imagen_data_noprs$Life_BSL
imagen_data_noprs$AUDIT_change<-imagen_data_noprs$AUDIT_C_FU2-imagen_data_noprs$AUDIT_C_BSL
imagen_data_noprs$Life_Nicotine_change<-imagen_data_noprs$Life_Nicotine_FU2-imagen_data_noprs$Life_Nicotine_BSL

#MAKE PCA SCORES OF BRAIN VARS
brain.pca.bl <- PCA(imagen_data_noprs[,c("right_dmPFC_BSL","left_dmPFC_BSL")], scale.unit=TRUE, graph=T)
imagen_data_noprs$brain.pca.bl<-brain.pca.bl$ind$coord[,1]
brain.pca.fu2 <- PCA(imagen_data_noprs[,c("right_dmPFC_FUP","left_dmPFC_FUP")], scale.unit=TRUE, graph=T)
imagen_data_noprs$brain.pca.fu2<-brain.pca.fu2$ind$coord[,1]
brain.pca.change <- PCA(imagen_data_noprs[,c("rdmpfc_change","ldmpfc_change")], scale.unit=TRUE, graph=T)
imagen_data_noprs$brain.pca.change<-brain.pca.change$ind$coord[,1]
brain.pca.spc <- PCA(imagen_data_noprs[,c("spc_left_dmPFC","spc_right_dmPFC")], scale.unit=TRUE, graph=T)
imagen_data_noprs$brain.pca.spc<-brain.pca.spc$ind$coord[,1]

#TURN VARIABLES OF INTEREST INTO RESIDUALS
#brain interecept (BL)
model1<-lmer(brain.pca.bl~ (1|Site) + (1|Sex), data=imagen_data_noprs)
Baseline_dPFC_Thickness <- as.matrix(residuals(model1))

#brain change resid
model1<-lmer(brain.pca.change~ (1|Site) + (1|Sex), data=imagen_data_noprs)
Change_in_dPFC_Thickness <- as.matrix(residuals(model1))

#canabis growth
model1<-lmer(Life_change~ (1|Site) + (1|Sex), data=imagen_data_noprs)
Cannabis_Use_by_Age_19<- as.matrix(residuals(model1))

#alcohol intercept
model1<-lmer(AUDIT_C_BSL~ (1|Site) + (1|Sex), data=imagen_data_noprs)
Baseline_Alcohol_Use<- as.matrix(residuals(model1))

#alcohol growth
model1<-lmer(AUDIT_change~ (1|Site) + (1|Sex), data=imagen_data_noprs)
Change_in_Alcohol_Use<- as.matrix(residuals(model1))

#nicotine intercept
model1<-lmer(Life_Nicotine_BSL~ (1|Site) + (1|Sex), data=imagen_data_noprs)
Baseline_Tobacco_Use<- as.matrix(residuals(model1))

#nicotine growth
model1<-lmer(Life_Nicotine_change~ (1|Site) + (1|Sex), data=imagen_data_noprs)
Change_in_Tobacco_Use<- as.matrix(residuals(model1))

#age
model1<-lmer(Baseline_Age~ (1|Site) + (1|Sex), data=imagen_data_noprs)
Baseline_Age<- as.matrix(residuals(model1))

#handedness
model1<-lmer(Handedness~(1|Site) + (1|Sex), data=imagen_data_noprs)
Handedness<- as.matrix(residuals(model1))

#make final BN network dataset
data.regressed_noprs <- cbind(Baseline_dPFC_Thickness,Change_in_dPFC_Thickness,Cannabis_Use_by_Age_19,Baseline_Tobacco_Use,Change_in_Tobacco_Use,Baseline_Alcohol_Use,Change_in_Alcohol_Use,Handedness,Baseline_Age)## assembles the residuals

data.regressed_noprs <- as.data.frame(scale(data.regressed_noprs, center = T, scale = T))## center+scale the data 
colnames(data.regressed_noprs) <- c('Baseline_dPFC_Thickness','Change_in_dPFC_Thickness','Cannabis_Use_by_Age_19','Baseline_Tobacco_Use','Change_in_Tobacco_Use','Baseline_Alcohol_Use','Change_in_Alcohol_Use','Handedness','Baseline_Age')

## generate the blacklist of the directions that you are not interested in
blacklist = data.frame(
  from =c(
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Baseline_dPFC_Thickness',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Baseline_dPFC_Thickness'
  ), 
  to = c(
    'Baseline_dPFC_Thickness',
    'Baseline_dPFC_Thickness',
    'Baseline_dPFC_Thickness',
    'Baseline_dPFC_Thickness',
    'Baseline_Alcohol_Use',
    'Baseline_Alcohol_Use',
    'Baseline_Alcohol_Use',
    'Baseline_Alcohol_Use',
    'Baseline_Tobacco_Use',
    'Baseline_Tobacco_Use',
    'Baseline_Tobacco_Use',
    'Baseline_Tobacco_Use',
    'Handedness',
    'Handedness',
    'Handedness',
    'Handedness',
    'Handedness',
    'Baseline_Age',
    'Baseline_Age',
    'Baseline_Age',
    'Baseline_Age',
    'Baseline_Age'
  )) 

## generate the blacklist of the directions that you are not interested in
blacklist_Brain2Cannabis = data.frame(
  from =c(
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Baseline_dPFC_Thickness',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Baseline_dPFC_Thickness',
    'Change_in_dPFC_Thickness'
  ), 
  to = c(
    'Baseline_dPFC_Thickness',
    'Baseline_dPFC_Thickness',
    'Baseline_dPFC_Thickness',
    'Baseline_dPFC_Thickness',
    'Baseline_Alcohol_Use',
    'Baseline_Alcohol_Use',
    'Baseline_Alcohol_Use',
    'Baseline_Alcohol_Use',
    'Baseline_Tobacco_Use',
    'Baseline_Tobacco_Use',
    'Baseline_Tobacco_Use',
    'Baseline_Tobacco_Use',
    'Handedness',
    'Handedness',
    'Handedness',
    'Handedness',
    'Handedness',
    'Baseline_Age',
    'Baseline_Age',
    'Baseline_Age',
    'Baseline_Age',
    'Baseline_Age',
    'Cannabis_Use_by_Age_19'
  )) 

## generate the blacklist of the directions that you are not interested in
blacklist_Cannabis2Brain = data.frame(
  from =c(
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Baseline_dPFC_Thickness',
    'Change_in_Alcohol_Use',
    'Change_in_dPFC_Thickness',
    'Cannabis_Use_by_Age_19',
    'Change_in_Tobacco_Use',
    'Baseline_dPFC_Thickness',
    'Cannabis_Use_by_Age_19'
  ), 
  to = c(
    'Baseline_dPFC_Thickness',
    'Baseline_dPFC_Thickness',
    'Baseline_dPFC_Thickness',
    'Baseline_dPFC_Thickness',
    'Baseline_Alcohol_Use',
    'Baseline_Alcohol_Use',
    'Baseline_Alcohol_Use',
    'Baseline_Alcohol_Use',
    'Baseline_Tobacco_Use',
    'Baseline_Tobacco_Use',
    'Baseline_Tobacco_Use',
    'Baseline_Tobacco_Use',
    'Handedness',
    'Handedness',
    'Handedness',
    'Handedness',
    'Handedness',
    'Baseline_Age',
    'Baseline_Age',
    'Baseline_Age',
    'Baseline_Age',
    'Baseline_Age',
    'Change_in_dPFC_Thickness'
  )) 

### the following is using different algorithms provided by the bnlearn package to get the BCN, 

# bootstrap hill climbing

#Use BIC for scoring
boot = boot.strength(data = data.regressed_noprs, R = 10000, algorithm = "hc", 
                     algorithm.args = list(start=NULL,score ="bic-g", blacklist = blacklist))

boot[(boot$strength > 0.5) & (boot$direction >= 0.5), ]

avg.boot_hc = averaged.network(boot, threshold = 0.80)
plot(avg.boot_hc,main="Hill Climbing")

##aic instead of bic
boot = boot.strength(data = data.regressed_noprs, R = 10000, algorithm = "hc", 
                     algorithm.args = list(start=NULL,score ="aic-g", blacklist = blacklist))

boot[(boot$strength > 0.8) & (boot$direction >= 0.5), ]

#set up random starting places for the hill climb
random_graphs_starts=function(num){
  start=random.graph(nodes = names(data.regressed_noprs), num = num)
  return(start)
}

bn_multi_start_str=function(start){
  boot=boot.strength(data = data.regressed_noprs, R = 100, algorithm = "hc", 
                     algorithm.args = list(start=start,score ="bic-g",  
                                           blacklist = blacklist))
  return(boot[24,3])
}


bn_multi_start_dir=function(start){
  boot=boot.strength(data = data.regressed_noprs, R = 100, algorithm = "hc", 
                     algorithm.args = list(start=start,score ="bic-g",  
                                           blacklist = blacklist))
  return(boot[24,4])
}

start=random_graphs_starts(100)

direction=sapply(start,bn_multi_start_dir)
strength=sapply(start,bn_multi_start_str)

mean(direction)
mean(strength)

#run hill climb with random perturbations introducted
boot = boot.strength(data = data.regressed_noprs, R = 10000, algorithm = "hc", 
                     algorithm.args = list(start=NULL,score ="bic-g", blacklist = blacklist,perturb=10))

boot[(boot$strength > 0.8) & (boot$direction >= 0.5), ]

################test with one direction blacklisted and see which does better#######
#block brain to cannabis directed arc
boot = boot.strength(data = data.regressed_noprs, R = 10000, algorithm = "hc", 
                     algorithm.args = list(start=NULL,score ="bic-g", blacklist = blacklist_Cannabis2Brain))
boot[(boot$strength > 0.9) & (boot$direction >= 0.5), ]

bn.hc.1 <- hc(data.regressed_noprs, blacklist = blacklist_Cannabis2Brain)

score(bn.hc.1,data.regressed_noprs, type='bic-g') 

#block cannabis to brain directed arc
boot = boot.strength(data = data.regressed_noprs, R = 10000, algorithm = "hc", 
                     algorithm.args = list(start=NULL,score ="bic-g", blacklist = blacklist_Brain2Cannabis))
boot[(boot$strength > 0.9) & (boot$direction >= 0.5), ]

bn.hc.2 <- hc(data.regressed_noprs, blacklist = blacklist_Brain2Cannabis)
plot(bn.hc)
bn.hc

score(bn.hc.2,data.regressed_noprs, type='bic-g')

###############tabu (a score based method based on hill climb)###########

#test some new scoring methods
boot = boot.strength(data = data.regressed_noprs, R = 10000, algorithm = "tabu", 
                     algorithm.args = list(start=NULL,score ="bic-g", blacklist = blacklist))

boot[(boot$strength > 0.9) & (boot$direction >= 0.5), ]

avg.boot_hc = averaged.network(boot, threshold = 0.90)
plot(avg.boot_hc,main="Tabu")

########### Grow shrink###########

#primary grow shrink
boot = boot.strength(data = data.regressed_noprs, R = 10000, algorithm = "gs", 
                     algorithm.args = list(blacklist = blacklist, test= 'mi-g'))

boot[(boot$strength > 0.8) & (boot$direction >= 0.5), ]

warnings()

avg.boot_gs = averaged.network(boot, threshold = 0.90)
plot(avg.boot_gs,main="Grow Shrink")


###GS with different independence tests
boot = boot.strength(data = data.regressed_noprs, R = 10000, algorithm = "gs", debug = FALSE,
                     algorithm.args = list(blacklist = blacklist, test= 'mc-zf', alpha = 0.05, B=10))

boot[(boot$strength > 0.8) & (boot$direction >= 0.7), ]

boot = boot.strength(data = data.regressed_noprs, R = 10000, algorithm = "gs", debug = FALSE,
                     algorithm.args = list(blacklist = blacklist, test= 'cor', alpha = 0.05, B=10))

boot[(boot$strength > 0.8) & (boot$direction >= 0.7), ]

avg.boot_gs = averaged.network(boot, threshold = 0.90)
plot(avg.boot_gs,main="Grow Shrink")

# bootstrap iamb

boot = boot.strength(data = data.regressed_noprs, R = 10000, algorithm = "iamb", 
                     algorithm.args = list(blacklist = blacklist))

boot[(boot$strength > 0.8) & (boot$direction >= 0.5), ]

avg.boot_iamb = averaged.network(boot, threshold = 0.90)
plot(avg.boot_iamb,main="IAMB")

################ min max hill climb (hybrid algorithm)#######

bn.mmhc <- mmhc(data.regressed_noprs, blacklist = blacklist)
plot(bn.mmhc)

boot = boot.strength(data = data.regressed_noprs, R = 10000, algorithm = "mmhc", 
                     algorithm.args = list(blacklist = blacklist))
warnings()

boot[(boot$strength > 0.7) & (boot$direction >= 0.5), ]

###############rsmax2 (hybrid)###########

#test some new scoring methods
boot = boot.strength(data = data.regressed_noprs, R = 10000, algorithm = "rsmax2", 
                     algorithm.args = list(blacklist = blacklist))

boot[(boot$strength > 0.7) & (boot$direction >= 0.5), ]
avg.boot_hc = averaged.network(boot, threshold = 0.90)
plot(avg.boot_hc,main="Max Min Hill Climbing")