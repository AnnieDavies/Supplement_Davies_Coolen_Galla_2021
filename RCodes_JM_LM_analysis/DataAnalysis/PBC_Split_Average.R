# Code to perform data analysis on PBC data set in Davies, Coolen and Galla (2021).
# Written by A. Davies (2021),
# aided by code by D. Rizopoulos: https://github.com/drizopoulos/jm_and_lm
# Data set is split in half 20 times into a training and test data set.
# The data is saved at each iteration to be read in to the C++ codes for the 
# retarded kernel models.
# Analysis is performed for data that treats the transplant event as a 
# censoring event (labelled 'cen') and for data that treats the two events (death
# and transplant) as a composite event (labelled 'comp').
# For each version (cen and comp) we fit two joint models (one with a linear
# longitudinal model and one using cubic splines).
# Joint models are fitted to the training data set.
# At each landmark (base) time landmark models are fitted to the training data.
# First we perform the fixed base time analysis with base time t = 3 years:
# The prediction time u is varied in steps of 0.2 years from 3 to 8 years
# Prediction error for each combination of t and u is stored at each iteration.
# Then we perform the fixed prediction window analysis for three windows:
# w1=1 year, w2=2 years, w3=3 years
# At w1 we use base times t=0->9 years in steps of 0.2 years
# At w2 we use base times t=0->8 years in steps of 0.2 years
# At w3 we use base times t=0->7 years in steps of 0.2 years
# Again, prediction error for each scenario is stored for each iteration.
# An average of PE over the 20 iterations is performed by C++ code Average_Split.cpp
# For me, this code took about 3 days to run.


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
# for mvJMbayes objects:
source("PE.AD.mvJM.R")

#load original data
#PBC data with all observations of covariates
data(pbc2, package="JMbayes")
#PBC data with only the baseline observations of covariates
data(pbc2.id, package="JMbayes")

#In pbc2 and pbc2.id, status2 indicates the event "death" treating the event
#'transplant' as a censoring event i.e.
pbc2.id$status2 <- as.numeric(pbc2.id$status == "dead")
pbc2$status2 <- as.numeric(pbc2$status == "dead")

#create the indicator status3 for the composite event (death or transplant)
pbc2.id$status3 <- as.numeric(pbc2.id$status != "alive")
pbc2$status3 <- as.numeric(pbc2$status != "alive")

#base folder in which to save test and training data in each loop
fn.base <- "~...\\PBC_DATA\\LOOP"

r <- 20 #number of repeats
n <- 312 #number of subjects in orig data
V <- 2 #number of ways to split data

# base time for fixed base time analysis
timeLM <- 3.0

# prediction windows for fixed window analysis
w1.pbc <- 1.0
w2.pbc <- 2.0
w3.pbc <- 3.0

set.seed(231)
for(k in 1:r){
  
  print(paste("*****************loop", k))
  print(Sys.time())
  
  #split data into training and test data---------------------------------------
  splits <- split(seq_len(n), sample(rep(seq_len(V), length.out = n)))
  vec.test.id <- unlist(splits[1], recursive = TRUE, use.names = TRUE)
  
  #create test and train data---------------------------------------------------
  
  #train data not in vec.test.id
  #train_data from pbc2 (all observations)
  train_data <- pbc2[!pbc2$id %in% vec.test.id, ]
  train_data$id <- match(train_data$id, unique(train_data$id))
  #train_data.id from pbc2.id (just baseline observations)
  train_data.id <- pbc2.id[!pbc2.id$id %in% vec.test.id, ]
  train_data.id$id <- match(train_data.id$id, unique(train_data.id$id))
  
  #test data is in vec.test.id
  #test_data from pbc2 (all observations)
  test_data <- pbc2[pbc2$id %in% vec.test.id, ]
  test_data$id <- match(test_data$id, unique(test_data$id))
  #test_data.id from pbc2.id (just baseline observations)
  test_data.id <- pbc2.id[pbc2.id$id %in% vec.test.id, ]
  test_data.id$id <- match(test_data.id$id, unique(test_data.id$id))
  
  #extract values from test and train data--------------------------------------
  #TRAIN:
  #number of observations per individual
  n.train <- table(train_data$id)  
  #fixed variables - use train_data.id
  events.train <-c(train_data.id[['years']])
  censor.train <-c(train_data.id[['status2']])
  censor.comp.train <-c(train_data.id[['status3']])
  Zf.age.train <-c(train_data.id[['age']])
  #longitudinal variables - use train data
  t.obs.train <- c(train_data[['year']])
  Zl.bili.train <-c(train_data[['serBilir']])
  Zl.alb.train <-c(train_data[['albumin']])
  Zl.proth.train <-c(train_data[['prothrombin']])
  
  #TEST:
  #number of observations per individual
  n.test <- table(test_data$id)  
  #fixed variables - use test_data.id
  events.test <-c(test_data.id[['years']])
  censor.test <-c(test_data.id[['status2']])
  censor.comp.test <-c(test_data.id[['status3']])
  Zf.age.test <-c(test_data.id[['age']])
  #longitudinal variables - use test data
  t.obs.test <- c(test_data[['year']])
  Zl.bili.test <-c(test_data[['serBilir']])
  Zl.alb.test <-c(test_data[['albumin']])
  Zl.proth.test <-c(test_data[['prothrombin']])
  
  #create filenames and write to file-------------------------------------------
  #TEST ID (from split)
  fn.test.id <- paste(fn.base, k, "\\test_ids.txt", sep='')
  write(vec.test.id, fn.test.id, ncolumns = 1)
  
  #TRAIN
  fn.n.train <- paste(fn.base, k, "\\TRAIN\\n_obs.txt", sep='')
  fn.events.train <- paste(fn.base, k, "\\TRAIN\\Events.txt", sep='')
  fn.censor.train <- paste(fn.base, k, "\\TRAIN\\Censor.txt", sep='')
  fn.censor.comp.train <- paste(fn.base, k, "\\TRAIN\\Censor_comp.txt", sep='')
  fn.Zf.age.train <- paste(fn.base, k, "\\TRAIN\\Zfix_age.txt", sep='')
  fn.t.obs.train <- paste(fn.base, k, "\\TRAIN\\t_obs.txt", sep='')
  fn.Zl.bili.train <- paste(fn.base, k, "\\TRAIN\\Zlong_bili.txt", sep='')
  fn.Zl.alb.train <- paste(fn.base, k, "\\TRAIN\\Zlong_alb.txt", sep='')
  fn.Zl.proth.train <- paste(fn.base, k, "\\TRAIN\\Zlong_proth.txt", sep='')
  
  write(formatC(n.train, digits=15, format="fg", flag = "-"), fn.n.train, ncolumns=1)
  write(formatC(events.train, digits=15, format="fg", flag = "-"), fn.events.train, ncolumns=1)
  write(formatC(censor.train, digits=15, format="fg", flag = "-"), fn.censor.train, ncolumns=1)
  write(formatC(censor.comp.train, digits=15, format="fg", flag = "-"), fn.censor.comp.train, ncolumns=1)
  write(formatC(Zf.age.train, digits=15, format="fg", flag = "-"), fn.Zf.age.train, ncolumns=1)
  write(formatC(t.obs.train, digits=15, format="fg", flag = "-"), fn.t.obs.train, ncolumns=1)
  write(formatC(Zl.bili.train, digits=15, format="fg", flag = "-"), fn.Zl.bili.train, ncolumns=1)
  write(formatC(Zl.alb.train, digits=15, format="fg", flag = "-"), fn.Zl.alb.train, ncolumns=1)
  write(formatC(Zl.proth.train, digits=15, format="fg", flag = "-"), fn.Zl.proth.train, ncolumns=1)
  
  #TEST
  fn.n.test <- paste(fn.base, k, "\\TEST\\n_obs.txt", sep='')
  fn.events.test <- paste(fn.base, k, "\\TEST\\Events.txt", sep='')
  fn.censor.test <- paste(fn.base, k, "\\TEST\\Censor.txt", sep='')
  fn.censor.comp.test <- paste(fn.base, k, "\\TEST\\Censor_comp.txt", sep='')
  fn.Zf.age.test <- paste(fn.base, k, "\\TEST\\Zfix_age.txt", sep='')
  fn.t.obs.test <- paste(fn.base, k, "\\TEST\\t_obs.txt", sep='')
  fn.Zl.bili.test <- paste(fn.base, k, "\\TEST\\Zlong_bili.txt", sep='')
  fn.Zl.alb.test <- paste(fn.base, k, "\\TEST\\Zlong_alb.txt", sep='')
  fn.Zl.proth.test <- paste(fn.base, k, "\\TEST\\Zlong_proth.txt", sep='')
  
  write(formatC(n.test, digits=15, format="fg", flag = "-"), fn.n.test, ncolumns=1)
  write(formatC(events.test, digits=15, format="fg", flag = "-"), fn.events.test, ncolumns=1)
  write(formatC(censor.test, digits=15, format="fg", flag = "-"), fn.censor.test, ncolumns=1)
  write(formatC(censor.comp.test, digits=15, format="fg", flag = "-"), fn.censor.comp.test, ncolumns=1)
  write(formatC(Zf.age.test, digits=15, format="fg", flag = "-"), fn.Zf.age.test, ncolumns=1)
  write(formatC(t.obs.test, digits=15, format="fg", flag = "-"), fn.t.obs.test, ncolumns=1)
  write(formatC(Zl.bili.test, digits=15, format="fg", flag = "-"), fn.Zl.bili.test, ncolumns=1)
  write(formatC(Zl.alb.test, digits=15, format="fg", flag = "-"), fn.Zl.alb.test, ncolumns=1)
  write(formatC(Zl.proth.test, digits=15, format="fg", flag = "-"), fn.Zl.proth.test, ncolumns=1)
  
  print("training and test data created and saved")
  print(Sys.time())
  
  #-----------------------------------------------------------------------------
  
 
  ################
  # Joint Model
  ################
  
  #fit Joint models using train_data
  
  # Longitudinal models --------------------------------------------------------
  # spline longitudinal model
  print("spline longitudinal model begin")
  long.spline <- mvglmer(list(log(serBilir) ~ ns(year,2,B=c(0,14.4)) + (ns(year,2,B=c(0,14.4)) | id),
                              log(albumin) ~ ns(year,2,B=c(0,14.4)) + (ns(year,2,B=c(0,14.4)) | id),
                              log(prothrombin) ~ ns(year,2,B=c(0,14.4)) + (ns(year,2,B=c(0,14.4)) | id)),
                         data = train_data, families = list(gaussian, gaussian, gaussian))
  print("done")
  print(Sys.time())
  
  # linear longitudinal model
  print("linear longitudinal model begin")
  long.linear <- mvglmer(list(log(serBilir) ~ year + (year | id),
                              log(albumin) ~ year + (year | id),
                              log(prothrombin) ~ year + (year | id)),
                         data = train_data, families = list(gaussian, gaussian, gaussian))
  print("done")
  print(Sys.time())
  
  # Survival models ------------------------------------------------------------
  print("survival models")
  
  # censor event survival model
  surv.cen <- coxph(Surv(years, status2) ~ age, data = train_data.id, model = TRUE)
  
  # composite event survival model
  surv.comp <- coxph(Surv(years, status3) ~ age, data = train_data.id, model = TRUE)
  
  print("done")
  print(Sys.time())
  
  # Joint models ---------------------------------------------------------------
  # JM: Spline Censor event
  print("JM: Spline Censor event begin")
  JM.spline.cen <- mvJointModelBayes(long.spline, surv.cen, timeVar = "year")
  print("done")
  print(Sys.time())
  
  # JM: Linear Censor event
  print("JM: Linear Censor event begin")
  JM.linear.cen <- mvJointModelBayes(long.linear, surv.cen, timeVar = "year")
  print("done")
  print(Sys.time())
  
  # JM: Spline Comp event
  print("JM: Spline Comp event begin")
  JM.spline.comp <- mvJointModelBayes(long.spline, surv.comp, timeVar = "year")
  print("done")
  print(Sys.time())
  
  # JM: Linear Comp event
  print("JM: Linear Comp event begin")
  JM.linear.comp <- mvJointModelBayes(long.linear, surv.comp, timeVar = "year")
  print("done")
  print(Sys.time())
  
  ###########################
  # Fix t
  ###########################
  
  print("begin fix t")
  
  ################
  # Landmarking 
  ################
  
  # fit landmark models to training data at the fixed base time (t=3)
  
  # Create a landmark data set at the landmark time (timeLM) using training data
  # We fit three longitudinal covariates (serBilir, albumin and prothrombin)
  # because we use summary = "value" we need only specify one longitudinal covariate
  # (arbitrarily) as 'respVar' (here we use serBilir)
  # NB: this would not hold for other arguments for summary
  pbc.LM <- dataLM(train_data, timeLM, respVar = "serBilir", timeVar = "year", 
                   idVar = "id", evTimeVar = "years", summary = "value")
  
  #we use the same LM data set for both the censor event and composite event models
  
  ################
  # Landmarking CENSOR EVENT
  ################
  
  print("LM censor model begin")
  
  #Fit a standard Cox model to the landmark data set (for the censor event)
  Cox.LM.cen <- coxph(Surv(years, status2) ~ age + log(serBilir) + log(albumin) + 
                        log(prothrombin), data = pbc.LM)
  
  print("done")
  
  ################
  # Landmarking COMP EVENT
  ################
  
  print("LM comp model begin")
  
  #Fit a standard Cox model to the landmark data set (for the composite event)
  Cox.LM.comp <- coxph(Surv(years, status3) ~ age + log(serBilir) + log(albumin) + 
                         log(prothrombin), data = pbc.LM)
  
  print("done")
  print(Sys.time())
  
 
  ################
  # Pred Error
  ################
  
  # for fixed t=3 years, iterate over prediction time u = 3->8 years in steps of 0.2
  # calculate PE for LMs and JMs for each combination of (t,u) using newdata = test data
  
  print("Fix t loop begin")
  
  vec.pe.LM.cen <-c()
  vec.pe.JM.spline.cen <-c()
  vec.pe.JM.linear.cen <-c()
  
  vec.pe.LM.comp <-c()
  vec.pe.JM.spline.comp <-c()
  vec.pe.JM.linear.comp <-c()
  
  for(i in 0:25){
    timeHZ <- timeLM + (i*0.2) #prediction or 'horizon' time u
    
    #censor event
    PE.LM.cen <- PE.AD.coxph(Cox.LM.cen, newdata = test_data, Tstart = timeLM, Thoriz = timeHZ, 
                            idVar = "id", timeVar = "year", respVar = "serBilir",
                            evTimeVar = "years", lossFun = "square", summary = "value")
    PE.JM.spline.cen <- PE.AD.mvJM(JM.spline.cen, newdata = test_data, Tstart = timeLM, 
                                  Thoriz = timeHZ, lossFun = "square", idVar = "id", simulate=TRUE)
    PE.JM.linear.cen <- PE.AD.mvJM(JM.linear.cen, newdata = test_data, Tstart = timeLM, 
                                  Thoriz = timeHZ, lossFun = "square", idVar = "id", simulate=TRUE)
    
    #comp event
    PE.LM.comp <- PE.AD.coxph(Cox.LM.comp, newdata = test_data, Tstart = timeLM, Thoriz = timeHZ, 
                              idVar = "id", timeVar = "year", respVar = "serBilir",
                              evTimeVar = "years", lossFun = "square", summary = "value")
    PE.JM.spline.comp <- PE.AD.mvJM(JM.spline.comp, newdata = test_data, Tstart = timeLM, 
                                    Thoriz = timeHZ, lossFun = "square", idVar = "id", simulate=TRUE)
    PE.JM.linear.comp <- PE.AD.mvJM(JM.linear.comp, newdata = test_data, Tstart = timeLM, 
                                    Thoriz = timeHZ, lossFun = "square", idVar = "id", simulate=TRUE)
    
    #censor event
    vec.pe.LM.cen <-c(vec.pe.LM.cen, PE.LM.cen["prederr"])
    vec.pe.JM.spline.cen <-c(vec.pe.JM.spline.cen, PE.JM.spline.cen["prederr"])
    vec.pe.JM.linear.cen <-c(vec.pe.JM.linear.cen, PE.JM.linear.cen["prederr"])
    #comp event
    vec.pe.LM.comp <-c(vec.pe.LM.comp, PE.LM.comp["prederr"])
    vec.pe.JM.spline.comp <-c(vec.pe.JM.spline.comp, PE.JM.spline.comp["prederr"])
    vec.pe.JM.linear.comp <-c(vec.pe.JM.linear.comp, PE.JM.linear.comp["prederr"])
    
  }
  # save vector of prediction errors to files 
  # censor event
  fn.fixt.LM.cen <- paste(fn.base, k, "\\LM_cen_fixt.txt", sep='')
  fn.fixt.JM.spline.cen <- paste(fn.base, k, "\\JM_spline_cen_fixt.txt", sep='')
  fn.fixt.JM.linear.cen <- paste(fn.base, k, "\\JM_linear_cen_fixt.txt", sep='')
  write(formatC(unlist(vec.pe.LM.cen), digits=15, format="fg", flag = "-"), fn.fixt.LM.cen, ncolumns=1)
  write(formatC(unlist(vec.pe.JM.spline.cen), digits=15, format="fg", flag = "-"), fn.fixt.JM.spline.cen, ncolumns=1)
  write(formatC(unlist(vec.pe.JM.linear.cen), digits=15, format="fg", flag = "-"), fn.fixt.JM.linear.cen, ncolumns=1)
  
  # comp event
  fn.fixt.LM.comp <- paste(fn.base, k, "\\LM_comp_fixt.txt", sep='')
  fn.fixt.JM.spline.comp <- paste(fn.base, k, "\\JM_spline_comp_fixt.txt", sep='')
  fn.fixt.JM.linear.comp <- paste(fn.base, k, "\\JM_linear_comp_fixt.txt", sep='')
  write(formatC(unlist(vec.pe.LM.comp), digits=15, format="fg", flag = "-"), fn.fixt.LM.comp, ncolumns=1)
  write(formatC(unlist(vec.pe.JM.spline.comp), digits=15, format="fg", flag = "-"), fn.fixt.JM.spline.comp, ncolumns=1)
  write(formatC(unlist(vec.pe.JM.linear.comp), digits=15, format="fg", flag = "-"), fn.fixt.JM.linear.comp, ncolumns=1)
  
  
  print("fix t loop done")
  print(Sys.time())
 
  
  ###########################
  # Pred Window
  ###########################
  
  ################
  # Window 1
  ################
  
  # For fixed window w1=1 year work out PE for base time t=0-9 years (steps of 0.2 years)
  # For each base time t we must create a new LM data set and fit corresponding
  # LM models (using the training data)
  
  print("window 1 begin")
  
  vec.pe.LM.cen1 <-c()
  vec.pe.JM.spline.cen1 <-c()
  vec.pe.JM.linear.cen1 <-c()
  
  vec.pe.LM.comp1 <-c()
  vec.pe.JM.spline.comp1 <-c()
  vec.pe.JM.linear.comp1 <-c()
  
  for(i in 0:45){
    t.pbc <- i*0.2
    
    #create LM data set from training data at time t.pbc
    pbc.LM1 <- dataLM(train_data, t.pbc, respVar = "serBilir", timeVar = "year", 
                      idVar = "id", evTimeVar = "years", summary = "value")
    
    #censor event
    #Fit a standard Cox model to the landmark data set (for the censor event)
    Cox.LM.cen1 <- coxph(Surv(years, status2) ~ age + log(serBilir) + log(albumin) + 
                                                log(prothrombin), data = pbc.LM1)
    PE.LM.cen1 <- PE.AD.coxph(Cox.LM.cen1, newdata = test_data, Tstart = t.pbc, Thoriz = t.pbc + w1.pbc, 
                              idVar = "id", timeVar = "year", respVar = "serBilir",
                              evTimeVar = "years", lossFun = "square", summary = "value")
    PE.JM.spline.cen1 <- PE.AD.mvJM(JM.spline.cen, newdata = test_data, Tstart = t.pbc, 
                                   Thoriz = t.pbc + w1.pbc, lossFun = "square", idVar = "id", simulate=TRUE)
    PE.JM.linear.cen1 <- PE.AD.mvJM(JM.linear.cen, newdata = test_data, Tstart = t.pbc, 
                                   Thoriz = t.pbc + w1.pbc, lossFun = "square", idVar = "id", simulate=TRUE)
    
    #comp event
    #Fit a standard Cox model to the landmark data set (for the composite event)
    Cox.LM.comp1 <- coxph(Surv(years, status3) ~ age + log(serBilir) + log(albumin) + 
                                                 log(prothrombin), data = pbc.LM1)
    PE.LM.comp1 <- PE.AD.coxph(Cox.LM.comp1, newdata = test_data, Tstart = t.pbc, Thoriz = t.pbc + w1.pbc, 
                              idVar = "id", timeVar = "year", respVar = "serBilir",
                              evTimeVar = "years", lossFun = "square", summary = "value")
    PE.JM.spline.comp1 <- PE.AD.mvJM(JM.spline.comp, newdata = test_data, Tstart = t.pbc, 
                                    Thoriz = t.pbc + w1.pbc, lossFun = "square", idVar = "id", simulate=TRUE)
    PE.JM.linear.comp1 <- PE.AD.mvJM(JM.linear.comp, newdata = test_data, Tstart = t.pbc, 
                                    Thoriz = t.pbc + w1.pbc, lossFun = "square", idVar = "id", simulate=TRUE)
    
    #censor event
    vec.pe.LM.cen1 <-c(vec.pe.LM.cen1, PE.LM.cen1["prederr"])
    vec.pe.JM.spline.cen1 <-c(vec.pe.JM.spline.cen1, PE.JM.spline.cen1["prederr"])
    vec.pe.JM.linear.cen1 <-c(vec.pe.JM.linear.cen1, PE.JM.linear.cen1["prederr"])
    #comp event
    vec.pe.LM.comp1 <-c(vec.pe.LM.comp1, PE.LM.comp1["prederr"])
    vec.pe.JM.spline.comp1 <-c(vec.pe.JM.spline.comp1, PE.JM.spline.comp1["prederr"])
    vec.pe.JM.linear.comp1 <-c(vec.pe.JM.linear.comp1, PE.JM.linear.comp1["prederr"])
    
  }
  # save vector of prediction errors to files
  # censor event
  fn.w1.LM.cen1 <- paste(fn.base, k, "\\LM_cen_w1.txt", sep='')
  fn.w1.JM.spline.cen1 <- paste(fn.base, k, "\\JM_spline_cen_w1.txt", sep='')
  fn.w1.JM.linear.cen1 <- paste(fn.base, k, "\\JM_linear_cen_w1.txt", sep='')
  write(formatC(unlist(vec.pe.LM.cen1), digits=15, format="fg", flag = "-"), fn.w1.LM.cen1, ncolumns=1)
  write(formatC(unlist(vec.pe.JM.spline.cen1), digits=15, format="fg", flag = "-"), fn.w1.JM.spline.cen1, ncolumns=1)
  write(formatC(unlist(vec.pe.JM.linear.cen1), digits=15, format="fg", flag = "-"), fn.w1.JM.linear.cen1, ncolumns=1)
  
  # comp event
  fn.w1.LM.comp1 <- paste(fn.base, k, "\\LM_comp_w1.txt", sep='')
  fn.w1.JM.spline.comp1 <- paste(fn.base, k, "\\JM_spline_comp_w1.txt", sep='')
  fn.w1.JM.linear.comp1 <- paste(fn.base, k, "\\JM_linear_comp_w1.txt", sep='')
  write(formatC(unlist(vec.pe.LM.comp1), digits=15, format="fg", flag = "-"), fn.w1.LM.comp1, ncolumns=1)
  write(formatC(unlist(vec.pe.JM.spline.comp1), digits=15, format="fg", flag = "-"), fn.w1.JM.spline.comp1, ncolumns=1)
  write(formatC(unlist(vec.pe.JM.linear.comp1), digits=15, format="fg", flag = "-"), fn.w1.JM.linear.comp1, ncolumns=1)
  
  
  print("w1 loop done")
  print(Sys.time())
  
  ################
  # Window 2
  ################
  
  # For fixed window w2=2 years work out PE for base time t=0-8 years (steps of 0.2 years)
  # For each base time t we must create a new LM data set and fit corresponding
  # LM models (using the training data)
  
  
  print("window 2 begin")
  
  vec.pe.LM.cen2 <-c()
  vec.pe.JM.spline.cen2 <-c()
  vec.pe.JM.linear.cen2 <-c()
  
  vec.pe.LM.comp2 <-c()
  vec.pe.JM.spline.comp2 <-c()
  vec.pe.JM.linear.comp2 <-c()
  
  for(i in 0:40){
    t.pbc <- i*0.2
    
    #create LM data set from training data at time t.pbc
    pbc.LM2 <- dataLM(train_data, t.pbc, respVar = "serBilir", timeVar = "year", 
                      idVar = "id", evTimeVar = "years", summary = "value")
    
    #censor event
    #Fit a standard Cox model to the landmark data set (for the censor event)
    Cox.LM.cen2 <- coxph(Surv(years, status2) ~ age + log(serBilir) + log(albumin) + 
                                                log(prothrombin), data = pbc.LM2)
    PE.LM.cen2 <- PE.AD.coxph(Cox.LM.cen2, newdata = test_data, Tstart = t.pbc, Thoriz = t.pbc + w2.pbc, 
                              idVar = "id", timeVar = "year", respVar = "serBilir",
                              evTimeVar = "years", lossFun = "square", summary = "value")
    PE.JM.spline.cen2 <- PE.AD.mvJM(JM.spline.cen, newdata = test_data, Tstart = t.pbc, 
                                   Thoriz = t.pbc + w2.pbc, lossFun = "square", idVar = "id", simulate=TRUE)
    PE.JM.linear.cen2 <- PE.AD.mvJM(JM.linear.cen, newdata = test_data, Tstart = t.pbc, 
                                   Thoriz = t.pbc + w2.pbc, lossFun = "square", idVar = "id", simulate=TRUE)
    
    #comp event
    #Fit a standard Cox model to the landmark data set (for the composite event)
    Cox.LM.comp2 <- coxph(Surv(years, status3) ~ age + log(serBilir) + log(albumin) + 
                                                log(prothrombin), data = pbc.LM2)
    PE.LM.comp2 <- PE.AD.coxph(Cox.LM.comp2, newdata = test_data, Tstart = t.pbc, Thoriz = t.pbc + w2.pbc, 
                              idVar = "id", timeVar = "year", respVar = "serBilir",
                              evTimeVar = "years", lossFun = "square", summary = "value")
    PE.JM.spline.comp2 <- PE.AD.mvJM(JM.spline.comp, newdata = test_data, Tstart = t.pbc, 
                                    Thoriz = t.pbc + w2.pbc, lossFun = "square", idVar = "id", simulate=TRUE)
    PE.JM.linear.comp2 <- PE.AD.mvJM(JM.linear.comp, newdata = test_data, Tstart = t.pbc, 
                                    Thoriz = t.pbc + w2.pbc, lossFun = "square", idVar = "id", simulate=TRUE)
    
    #censor event
    vec.pe.LM.cen2 <-c(vec.pe.LM.cen2, PE.LM.cen2["prederr"])
    vec.pe.JM.spline.cen2 <-c(vec.pe.JM.spline.cen2, PE.JM.spline.cen2["prederr"])
    vec.pe.JM.linear.cen2 <-c(vec.pe.JM.linear.cen2, PE.JM.linear.cen2["prederr"])
    #comp event
    vec.pe.LM.comp2 <-c(vec.pe.LM.comp2, PE.LM.comp2["prederr"])
    vec.pe.JM.spline.comp2 <-c(vec.pe.JM.spline.comp2, PE.JM.spline.comp2["prederr"])
    vec.pe.JM.linear.comp2 <-c(vec.pe.JM.linear.comp2, PE.JM.linear.comp2["prederr"])
    
  }
  # save vector of prediction errors to files
  # censor event
  fn.w2.LM.cen2 <- paste(fn.base, k, "\\LM_cen_w2.txt", sep='')
  fn.w2.JM.spline.cen2 <- paste(fn.base, k, "\\JM_spline_cen_w2.txt", sep='')
  fn.w2.JM.linear.cen2 <- paste(fn.base, k, "\\JM_linear_cen_w2.txt", sep='')
  write(formatC(unlist(vec.pe.LM.cen2), digits=15, format="fg", flag = "-"), fn.w2.LM.cen2, ncolumns=1)
  write(formatC(unlist(vec.pe.JM.spline.cen2), digits=15, format="fg", flag = "-"), fn.w2.JM.spline.cen2, ncolumns=1)
  write(formatC(unlist(vec.pe.JM.linear.cen2), digits=15, format="fg", flag = "-"), fn.w2.JM.linear.cen2, ncolumns=1)
  
  # comp event
  fn.w2.LM.comp2 <- paste(fn.base, k, "\\LM_comp_w2.txt", sep='')
  fn.w2.JM.spline.comp2 <- paste(fn.base, k, "\\JM_spline_comp_w2.txt", sep='')
  fn.w2.JM.linear.comp2 <- paste(fn.base, k, "\\JM_linear_comp_w2.txt", sep='')
  write(formatC(unlist(vec.pe.LM.comp2), digits=15, format="fg", flag = "-"), fn.w2.LM.comp2, ncolumns=1)
  write(formatC(unlist(vec.pe.JM.spline.comp2), digits=15, format="fg", flag = "-"), fn.w2.JM.spline.comp2, ncolumns=1)
  write(formatC(unlist(vec.pe.JM.linear.comp2), digits=15, format="fg", flag = "-"), fn.w2.JM.linear.comp2, ncolumns=1)
  
  
  print("w2 loop done")
  print(Sys.time())
  
  ################
  # Window 3
  ################
  
  # For fixed window w3=3 years work out PE for base time t=0-7 years (steps of 0.2 years)
  # For each base time t we must create a new LM data set and fit corresponding
  # LM models (using the training data)
  
  print("window 3 begin")
  
  vec.pe.LM.cen3 <-c()
  vec.pe.JM.spline.cen3 <-c()
  vec.pe.JM.linear.cen3 <-c()
  
  vec.pe.LM.comp3 <-c()
  vec.pe.JM.spline.comp3 <-c()
  vec.pe.JM.linear.comp3 <-c()
  
  for(i in 0:35){
    t.pbc <- i*0.2
    
    #create LM data set from training data at time t.pbc
    pbc.LM3 <- dataLM(train_data, t.pbc, respVar = "serBilir", timeVar = "year", 
                      idVar = "id", evTimeVar = "years", summary = "value")
    
    #censor event
    #Fit a standard Cox model to the landmark data set (for the censor event)
    Cox.LM.cen3 <- coxph(Surv(years, status2) ~ age + log(serBilir) + log(albumin) + 
                                                log(prothrombin), data = pbc.LM3)
    PE.LM.cen3 <- PE.AD.coxph(Cox.LM.cen3, newdata = test_data, Tstart = t.pbc, Thoriz = t.pbc + w3.pbc, 
                              idVar = "id", timeVar = "year", respVar = "serBilir",
                              evTimeVar = "years", lossFun = "square", summary = "value")
    PE.JM.spline.cen3 <- PE.AD.mvJM(JM.spline.cen, newdata = test_data, Tstart = t.pbc, 
                                   Thoriz = t.pbc + w3.pbc, lossFun = "square", idVar = "id", simulate=TRUE)
    PE.JM.linear.cen3 <- PE.AD.mvJM(JM.linear.cen, newdata = test_data, Tstart = t.pbc, 
                                   Thoriz = t.pbc + w3.pbc, lossFun = "square", idVar = "id", simulate=TRUE)
    
    #comp event
    #Fit a standard Cox model to the landmark data set (for the composite event)
    Cox.LM.comp3 <- coxph(Surv(years, status3) ~ age + log(serBilir) + log(albumin) + 
                                                log(prothrombin), data = pbc.LM3)
    PE.LM.comp3 <- PE.AD.coxph(Cox.LM.comp3, newdata = test_data, Tstart = t.pbc, Thoriz = t.pbc + w3.pbc, 
                              idVar = "id", timeVar = "year", respVar = "serBilir",
                              evTimeVar = "years", lossFun = "square", summary = "value")
    PE.JM.spline.comp3 <- PE.AD.mvJM(JM.spline.comp, newdata = test_data, Tstart = t.pbc, 
                                    Thoriz = t.pbc + w3.pbc, lossFun = "square", idVar = "id", simulate=TRUE)
    PE.JM.linear.comp3 <- PE.AD.mvJM(JM.linear.comp, newdata = test_data, Tstart = t.pbc, 
                                    Thoriz = t.pbc + w3.pbc, lossFun = "square", idVar = "id", simulate=TRUE)
    
    #censor event
    vec.pe.LM.cen3 <-c(vec.pe.LM.cen3, PE.LM.cen3["prederr"])
    vec.pe.JM.spline.cen3 <-c(vec.pe.JM.spline.cen3, PE.JM.spline.cen3["prederr"])
    vec.pe.JM.linear.cen3 <-c(vec.pe.JM.linear.cen3, PE.JM.linear.cen3["prederr"])
    #comp event
    vec.pe.LM.comp3 <-c(vec.pe.LM.comp3, PE.LM.comp3["prederr"])
    vec.pe.JM.spline.comp3 <-c(vec.pe.JM.spline.comp3, PE.JM.spline.comp3["prederr"])
    vec.pe.JM.linear.comp3 <-c(vec.pe.JM.linear.comp3, PE.JM.linear.comp3["prederr"])
    
  }
  # save vector of prediction errors to files
  # censor event
  fn.w3.LM.cen3 <- paste(fn.base, k, "\\LM_cen_w3.txt", sep='')
  fn.w3.JM.spline.cen3 <- paste(fn.base, k, "\\JM_spline_cen_w3.txt", sep='')
  fn.w3.JM.linear.cen3 <- paste(fn.base, k, "\\JM_linear_cen_w3.txt", sep='')
  write(formatC(unlist(vec.pe.LM.cen3), digits=15, format="fg", flag = "-"), fn.w3.LM.cen3, ncolumns=1)
  write(formatC(unlist(vec.pe.JM.spline.cen3), digits=15, format="fg", flag = "-"), fn.w3.JM.spline.cen3, ncolumns=1)
  write(formatC(unlist(vec.pe.JM.linear.cen3), digits=15, format="fg", flag = "-"), fn.w3.JM.linear.cen3, ncolumns=1)
  
  # comp event
  fn.w3.LM.comp3 <- paste(fn.base, k, "\\LM_comp_w3.txt", sep='')
  fn.w3.JM.spline.comp3 <- paste(fn.base, k, "\\JM_spline_comp_w3.txt", sep='')
  fn.w3.JM.linear.comp3 <- paste(fn.base, k, "\\JM_linear_comp_w3.txt", sep='')
  write(formatC(unlist(vec.pe.LM.comp3), digits=15, format="fg", flag = "-"), fn.w3.LM.comp3, ncolumns=1)
  write(formatC(unlist(vec.pe.JM.spline.comp3), digits=15, format="fg", flag = "-"), fn.w3.JM.spline.comp3, ncolumns=1)
  write(formatC(unlist(vec.pe.JM.linear.comp3), digits=15, format="fg", flag = "-"), fn.w3.JM.linear.comp3, ncolumns=1)
  
  
  print("w3 loop done")
  print(Sys.time())
  
}