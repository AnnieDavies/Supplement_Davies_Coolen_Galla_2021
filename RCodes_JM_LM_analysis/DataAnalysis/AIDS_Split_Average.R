# Code to perform data analysis on AIDS data set in Davies, Coolen and Galla (2021)
# for joint models and landmarking models.
# Written by A. Davies (2021),
# aided by code by D. Rizopoulos: https://github.com/drizopoulos/jm_and_lm
# Data set is split in half 20 times into a training and test data set.
# The data is saved at each iteration to be read in to the C++ codes for the 
# retarded kernel models.
# A joint model is fitted to the training data set.
# At each landmark (base) time a landmark model is fitted to the training data.
# First we perform the fixed base time analysis with base time t = 6 months:
# The prediction time u is varied in steps of 0.2 months from 6 to 18 months
# Prediction error calculated using the fitted model and test data 
# for each combination of t and u is stored at each iteration.
# Then we perform the fixed prediction window analysis for three windows:
# w1=6 months, w2=9 months, w3=12 months
# At w1 we use base times t=0,2,6,12 months
# At w2 & w3 we use base times t=0,2,6 months
# Again, prediction error for each scenario is stored for each iteration.
# An average of PE over the 20 iterations is performed by C++ code Average_Split.cpp
# For me, this code took about 10 hours to run.

library("JMbayes")
library("splines")
library("xtable")
library("lattice")

# internal JMbayes function to create landmark data sets
dataLM <- JMbayes:::dataLM 
# version of dataLM that creates a data set keeping event times >= landmark time for use
# in PE.AD.coxph only (i.e. to create the data set over which we sum in the PE eq):
source("dataLM.AD.R")
# edited versions of prederrJM (from JMbayes) so that it exactly matches the PE eq
# for coxph objects:
source("PE.AD.coxph.R")
# for JMbayes objects:
source("PE.AD.JM.R")

#load original data
#AIDS data with all observations of covariates
data(aids, package="JMbayes")
#AIDS data with only the baseline observations of covariates
data(aids.id, package="JMbayes")

#encode fixed covariates for saving data 
#drug (ddI=1 ddC=0)
aids$drug_code <- as.numeric(aids$drug == "ddI")
aids.id$drug_code <- as.numeric(aids.id$drug == "ddI")

#gender (male=1 female=0)
aids$gender_code <- as.numeric(aids$gender == "male")
aids.id$gender_code <- as.numeric(aids.id$gender == "male")

#PrevOI (AIDS=1 noAIDS=0)
aids$prevOI_code <- as.numeric(aids$prevOI == "AIDS")
aids.id$prevOI_code <- as.numeric(aids.id$prevOI == "AIDS")

#Stratum (AZT failure=1 intolerance=0)
aids$AZT_code <- as.numeric(aids$AZT == "failure")
aids.id$AZT_code <- as.numeric(aids.id$AZT == "failure")

#base folder in which to save test and training data in each loop
fn.base <- "~....\\AIDS_DATA\\LOOP"

r <- 20 #number of repeats
n <- 467 #number of subjects in orig data
V <- 2 #number of ways to split data

# base time for fixed base time analysis
timeLM.aids <- 6.0

# prediction windows for fixed window analysis
w1.aids <- 6.0
w2.aids <- 9.0
w3.aids <- 12.0


set.seed(312)
for(k in 1:r){
  
  print(paste("*****************loop", k))
  print(Sys.time())
  
  #split data into training and test data---------------------------------------
  splits <- split(seq_len(n), sample(rep(seq_len(V), length.out = n)))
  vec.test.id <- unlist(splits[1], recursive = TRUE, use.names = TRUE)
  
  #create test and train data---------------------------------------------------
  
  #train data: select ids not in vec.test.id
  #train_data from aids (all observations)
  train_data <- aids[!aids$patient %in% vec.test.id, ]
  train_data$patient <- match(train_data$patient, unique(train_data$patient))
  #train_data.id from aids.id (just baseline observations)
  train_data.id <- aids.id[!aids.id$patient %in% vec.test.id, ]
  train_data.id$patient <- match(train_data.id$patient, unique(train_data.id$patient))
  
  #test data: select ids in vec.test.id
  #test_data from aids (all observations)
  test_data <- aids[aids$patient %in% vec.test.id, ]
  test_data$patient <- match(test_data$patient, unique(test_data$patient))
  #test_data.id from aids.id (just baseline observations)
  test_data.id <- aids.id[aids.id$patient %in% vec.test.id, ]
  test_data.id$patient <- match(test_data.id$patient, unique(test_data.id$patient))
  
  #extract values from test and train data--------------------------------------
  #TRAIN:
  #number of observations per individual
  n.train <- table(train_data$patient)
  #fixed variables - use train_data.id
  events.train <-c(train_data.id[['Time']])
  censor.train <-c(train_data.id[['death']])
  Zf.drug.train <-c(train_data.id[['drug_code']])
  Zf.gender.train <-c(train_data.id[['gender_code']])
  Zf.prevOI.train <-c(train_data.id[['prevOI_code']])
  Zf.AZT.train <-c(train_data.id[['AZT_code']])
  #longitudinal variables - use train data
  t.obs.train <- c(train_data[['obstime']])
  Zl.CD4.train <-c(train_data[['CD4']])
  
  #TEST:
  #number of observations per individual 
  n.test <- table(test_data$patient)
  #fixed variables - use test_data.id
  events.test <-c(test_data.id[['Time']])
  censor.test <-c(test_data.id[['death']])
  Zf.drug.test <-c(test_data.id[['drug_code']])
  Zf.gender.test <-c(test_data.id[['gender_code']])
  Zf.prevOI.test <-c(test_data.id[['prevOI_code']])
  Zf.AZT.test <-c(test_data.id[['AZT_code']])
  #longitudinal variables - use test data
  t.obs.test <- c(test_data[['obstime']])
  Zl.CD4.test <-c(test_data[['CD4']]) 
  
  
  #create filenames and write to file-------------------------------------------
  #TEST ID (from split)
  fn.test.id <- paste(fn.base, k, "\\test_ids.txt", sep='')
  write(vec.test.id, fn.test.id, ncolumns = 1)
  
  #TRAIN
  fn.n.train <- paste(fn.base, k, "\\TRAIN\\n_obs.txt", sep='')
  fn.events.train <- paste(fn.base, k, "\\TRAIN\\Events.txt", sep='')
  fn.censor.train <- paste(fn.base, k, "\\TRAIN\\Censor.txt", sep='')
  fn.Zf.drug.train <- paste(fn.base, k, "\\TRAIN\\Zfix_drug.txt", sep='')
  fn.Zf.gender.train <- paste(fn.base, k, "\\TRAIN\\Zfix_gender.txt", sep='')
  fn.Zf.prevOI.train <- paste(fn.base, k, "\\TRAIN\\Zfix_prevOI.txt", sep='')
  fn.Zf.AZT.train <- paste(fn.base, k, "\\TRAIN\\Zfix_AZT.txt", sep='')
  fn.t.obs.train <- paste(fn.base, k, "\\TRAIN\\t_obs.txt", sep='')
  fn.Zl.CD4.train <- paste(fn.base, k, "\\TRAIN\\Zlong_CD4.txt", sep='')
  
  write(formatC(n.train, digits=15, format="fg", flag = "-"), fn.n.train, ncolumns=1)
  write(formatC(events.train, digits=15, format="fg", flag = "-"), fn.events.train, ncolumns=1)
  write(formatC(censor.train, digits=15, format="fg", flag = "-"), fn.censor.train, ncolumns=1)
  write(formatC(Zf.drug.train, digits=15, format="fg", flag = "-"), fn.Zf.drug.train, ncolumns=1)
  write(formatC(Zf.gender.train, digits=15, format="fg", flag = "-"), fn.Zf.gender.train, ncolumns=1)
  write(formatC(Zf.prevOI.train, digits=15, format="fg", flag = "-"), fn.Zf.prevOI.train, ncolumns=1)
  write(formatC(Zf.AZT.train, digits=15, format="fg", flag = "-"), fn.Zf.AZT.train, ncolumns=1)
  write(formatC(t.obs.train, digits=15, format="fg", flag = "-"), fn.t.obs.train, ncolumns=1)
  write(formatC(Zl.CD4.train, digits=15, format="fg", flag = "-"), fn.Zl.CD4.train, ncolumns=1)
  
  #TEST
  fn.n.test <- paste(fn.base, k, "\\TEST\\n_obs.txt", sep='')
  fn.events.test <- paste(fn.base, k, "\\TEST\\Events.txt", sep='')
  fn.censor.test <- paste(fn.base, k, "\\TEST\\Censor.txt", sep='')
  fn.Zf.drug.test <- paste(fn.base, k, "\\TEST\\Zfix_drug.txt", sep='')
  fn.Zf.gender.test <- paste(fn.base, k, "\\TEST\\Zfix_gender.txt", sep='')
  fn.Zf.prevOI.test <- paste(fn.base, k, "\\TEST\\Zfix_prevOI.txt", sep='')
  fn.Zf.AZT.test <- paste(fn.base, k, "\\TEST\\Zfix_AZT.txt", sep='')
  fn.t.obs.test <- paste(fn.base, k, "\\TEST\\t_obs.txt", sep='')
  fn.Zl.CD4.test <- paste(fn.base, k, "\\TEST\\Zlong_CD4.txt", sep='')
  
  write(formatC(n.test, digits=15, format="fg", flag = "-"), fn.n.test, ncolumns=1)
  write(formatC(events.test, digits=15, format="fg", flag = "-"), fn.events.test, ncolumns=1)
  write(formatC(censor.test, digits=15, format="fg", flag = "-"), fn.censor.test, ncolumns=1)
  write(formatC(Zf.drug.test, digits=15, format="fg", flag = "-"), fn.Zf.drug.test, ncolumns=1)
  write(formatC(Zf.gender.test, digits=15, format="fg", flag = "-"), fn.Zf.gender.test, ncolumns=1)
  write(formatC(Zf.prevOI.test, digits=15, format="fg", flag = "-"), fn.Zf.prevOI.test, ncolumns=1)
  write(formatC(Zf.AZT.test, digits=15, format="fg", flag = "-"), fn.Zf.AZT.test, ncolumns=1)
  write(formatC(t.obs.test, digits=15, format="fg", flag = "-"), fn.t.obs.test, ncolumns=1)
  write(formatC(Zl.CD4.test, digits=15, format="fg", flag = "-"), fn.Zl.CD4.test, ncolumns=1)
  
  print("training and test data created and saved")
  print(Sys.time())
  
  #-----------------------------------------------------------------------------
  
  ################
  # Joint Model
  ################
  
  print("JM model begin")
  
  #fit Joint model using train_data
  
  #longitudinal model
  long.train <- lme(CD4 ~ obstime + obstime:drug,
                    random = ~ obstime | patient, data = train_data)
  print("longitudinal done")
  print(Sys.time())
  #survival model
  Surv.train <- coxph(Surv(Time, death) ~ drug + prevOI + AZT + gender, 
                      data = train_data.id, x = TRUE)
  print("survival done")
  #joint model
  Joint.train <- jointModelBayes(long.train, Surv.train, timeVar = "obstime")
  print("joint done")
  print(Sys.time())
  
  ###########################
  # Fix t
  ###########################
  
  print("begin fix t")
  
  ################
  # Landmarking
  ################
  
  # fit landmark model to training data at the fixed base time (t=6)
  
  print("LM model begin")
  
  #Create a landmark data set at the landmark time (timeLM.aids) using training data
  train.LM <- dataLM(train_data, timeLM.aids, respVar = "CD4", timeVar = "obstime", 
                     evTimeVar = "Time", idVar = "patient", summary = "value")
  #Fit a standard Cox model to the landmark data set 
  Cox.train.LM <- coxph(Surv(Time, death) ~ drug + prevOI + AZT + gender + CD4, data = train.LM)
  
  print("done")
  print(Sys.time())
  
  ################
  # Pred Error
  ################
  
  # for fixed t=6 months, iterate over prediction time u = 6->18 months in steps of 0.2
  # calculate PE for LM and JM for each combination of (t,u) using newdata = test data
  
  print("Fix t loop begin")
  
  vec.pe.LM <-c()
  vec.pe.JM <-c()
  for(i in 0:60){
    timeHZ <- timeLM.aids+(i*0.2) #prediction or 'horizon' time u
    
    PE.LM<-PE.AD.coxph(Cox.train.LM, newdata = test_data, Tstart = timeLM.aids, Thoriz = timeHZ,
                  idVar = "patient", timeVar = "obstime", respVar = "CD4", evTimeVar = "Time",
                  lossFun = "square", summary = "value")
    
    PE.JM<-PE.AD.JM(Joint.train, newdata = test_data, Tstart = timeLM.aids, Thoriz = timeHZ, 
                  idVar = "patient", timeVar = "obstime", respVar = "CD4", evTimeVar = "Time",
                  lossFun = "square", summary = "value", simulate = TRUE)
    
    vec.pe.LM <-c(vec.pe.LM, PE.LM["prederr"])
    vec.pe.JM <-c(vec.pe.JM, PE.JM["prederr"])
  }
  #save vector of prediction errors to files
  fn.fixt.LM <- paste(fn.base, k, "\\LM_fixt.txt", sep='')
  fn.fixt.JM <- paste(fn.base, k, "\\JM_fixt.txt", sep='')
  write(formatC(unlist(vec.pe.LM), digits=15, format="fg", flag = "-"), fn.fixt.LM, ncolumns=1)
  write(formatC(unlist(vec.pe.JM), digits=15, format="fg", flag = "-"), fn.fixt.JM, ncolumns=1)
  
  print("fix t loop done")
  print(Sys.time())
 
  ###########################
  # Pred Window
  ###########################
  
  ################
  # Window 1
  ################
  
  # For fixed window w1=6 months work out PE for base time t=0,2,6,12 months
  # For each base time t we must create a new LM data set and fit a corresponding
  # LM model (using the training data)
  
  print("window 1 begin")
  
  vec.pe.LM1 <-c()
  vec.pe.JM1 <-c()
  
  #t = 0, 2, 6, 12
  for(i in 0:3){
    if(i<2){
      t.aids <- i*2
    } else if(i==2){
      t.aids <- 6
    }
    else if(i==3){
      t.aids <- 12
    }
    print(t.aids)
    
    #create LM data set at LM time = t.aids using training data
    aids.LM1 <- dataLM(train_data, t.aids, respVar = "CD4", timeVar = "obstime", 
                       evTimeVar = "Time", idVar = "patient", summary = "value")
    #Fit a standard Cox model to the landmark data set
    Cox.aids.LM1 <- coxph(Surv(Time, death) ~ drug + prevOI + AZT + gender + CD4, 
                          data = aids.LM1)
    
    #PE using test_data
    #LM
    PE.LM1<-PE.AD.coxph(Cox.aids.LM1, newdata = test_data, Tstart = t.aids, Thoriz = t.aids+w1.aids,
                     idVar = "patient", timeVar = "obstime", respVar = "CD4", evTimeVar = "Time",
                     lossFun = "square", summary = "value")
    #JM
    PE.JM1<-PE.AD.JM(Joint.train, newdata = test_data, Tstart = t.aids, Thoriz = t.aids+w1.aids, 
                     idVar = "patient", timeVar = "obstime", respVar = "CD4", evTimeVar = "Time",
                     lossFun = "square", summary = "value", simulate = TRUE)
    
    vec.pe.LM1 <-c(vec.pe.LM1, PE.LM1["prederr"])
    vec.pe.JM1 <-c(vec.pe.JM1, PE.JM1["prederr"])
  }
  #save vector of prediction errors to files
  fn.pw1.LM <- paste(fn.base, k, "\\LM_w1.txt", sep='')
  fn.pw1.JM <- paste(fn.base, k, "\\JM_w1.txt", sep='')
  write(formatC(unlist(vec.pe.LM1), digits=15, format="fg", flag = "-"), fn.pw1.LM, ncolumns=1)
  write(formatC(unlist(vec.pe.JM1), digits=15, format="fg", flag = "-"), fn.pw1.JM, ncolumns=1)
  print("done")
  print(Sys.time())
  
  ################
  # Window 2
  ################
  
  # For fixed window w2=9 months work out PE for base time t=0,2,6 months
  # For each base time t we must create a new LM data set and fit a corresponding
  # LM model (using the training data)
  
  print("window 2 begin")
  
  vec.pe.LM2 <-c()
  vec.pe.JM2 <-c()
  #t = 0, 2, 6
  for(i in 0:2){
    if(i<2){
      t.aids <- i*2
    } else if(i==2){
      t.aids <- 6
    }
    print(t.aids)
    
    #create LM data set at LM time = t.aids using training data
    aids.LM2 <- dataLM(train_data, t.aids, respVar = "CD4", timeVar = "obstime", 
                       evTimeVar = "Time", idVar = "patient", summary = "value")
    #Fit a standard Cox model to the landmark data set
    Cox.aids.LM2 <- coxph(Surv(Time, death) ~ drug + prevOI + AZT + gender + CD4, 
                          data = aids.LM2)
    
    #PE using test_data
    #LM
    PE.LM2<-PE.AD.coxph(Cox.aids.LM2, newdata = test_data, Tstart = t.aids, Thoriz = t.aids+w2.aids,
                      idVar = "patient", timeVar = "obstime", respVar = "CD4", evTimeVar = "Time",
                      lossFun = "square", summary = "value")
    #JM
    PE.JM2<-PE.AD.JM(Joint.train, newdata = test_data, Tstart = t.aids, Thoriz = t.aids+w2.aids, 
                      idVar = "patient", timeVar = "obstime", respVar = "CD4", evTimeVar = "Time",
                      lossFun = "square", summary = "value", simulate = TRUE)
    
    vec.pe.LM2 <-c(vec.pe.LM2, PE.LM2["prederr"])
    vec.pe.JM2 <-c(vec.pe.JM2, PE.JM2["prederr"])
  }
  #save vector of prediction errors to files
  fn.pw2.LM <- paste(fn.base, k, "\\LM_w2.txt", sep='')
  fn.pw2.JM <- paste(fn.base, k, "\\JM_w2.txt", sep='')
  write(formatC(unlist(vec.pe.LM2), digits=15, format="fg", flag = "-"), fn.pw2.LM, ncolumns=1)
  write(formatC(unlist(vec.pe.JM2), digits=15, format="fg", flag = "-"), fn.pw2.JM, ncolumns=1)
  print("done")
  print(Sys.time())
  
  ################
  # Window 3
  ################
  
  # For fixed window w3=12 months work out PE for base time t=0,2,6 months
  # For each base time t we must create a new LM data set and fit a corresponding
  # LM model (using the training data)
  
  print("window 3 begin")
  
  vec.pe.LM3 <-c()
  vec.pe.JM3 <-c()
  #t = 0, 2, 6
  for(i in 0:2){
    if(i<2){
      t.aids <- i*2
    } else if(i==2){
      t.aids <- 6
    }
    print(t.aids)
    
    #create LM data set at LM time = t.aids using training data
    aids.LM3 <- dataLM(train_data, t.aids, respVar = "CD4", timeVar = "obstime", 
                       evTimeVar = "Time", idVar = "patient", summary = "value")
    #Fit a standard Cox model to the landmark data set
    Cox.aids.LM3 <- coxph(Surv(Time, death) ~ drug + prevOI + AZT + gender + CD4, 
                          data = aids.LM3)
    
    #PE using test_data
    #LM
    PE.LM3<-PE.AD.coxph(Cox.aids.LM3, newdata = test_data, Tstart = t.aids, Thoriz = t.aids+w3.aids,
                      idVar = "patient", timeVar = "obstime", respVar = "CD4", evTimeVar = "Time",
                      lossFun = "square", summary = "value")
    #JM
    PE.JM3<-PE.AD.JM(Joint.train, newdata = test_data, Tstart = t.aids, Thoriz = t.aids+w3.aids, 
                      idVar = "patient", timeVar = "obstime", respVar = "CD4", evTimeVar = "Time",
                      lossFun = "square", summary = "value", simulate = TRUE)
    
    vec.pe.LM3 <-c(vec.pe.LM3, PE.LM3["prederr"])
    vec.pe.JM3 <-c(vec.pe.JM3, PE.JM3["prederr"])
  }
  #save vector of prediction errors to files
  fn.pw3.LM <- paste(fn.base, k, "\\LM_w3.txt", sep='')
  fn.pw3.JM <- paste(fn.base, k, "\\JM_w3.txt", sep='')
  write(formatC(unlist(vec.pe.LM3), digits=15, format="fg", flag = "-"), fn.pw3.LM, ncolumns=1)
  write(formatC(unlist(vec.pe.JM3), digits=15, format="fg", flag = "-"), fn.pw3.JM, ncolumns=1)
  print("done")
  print(Sys.time())
  
}
