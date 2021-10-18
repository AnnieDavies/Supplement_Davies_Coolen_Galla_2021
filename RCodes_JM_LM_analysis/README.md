R Codes to perform Data Analysis of Joint Models and Landmarking Models in Davies, Coolen and Galla (2021)

The folder DataAnalysis contains three R codes:
1. AIDS_Split_Average.R 
2. Liver_Split_Average.R
3. PBC_Split_Average.R

Each code reads in the original data set (aids, prothro and pbc2 respectively) and performs the fixed base time and fixed prediction window analysis for 20 random splits of the data into a training data set and a test data set. 

At each iteration the training and test data sets are saved so they can be read into the C++ codes for the retarded kernel models.

Models are fitted to the training data while prediction error is calculated for the test data for each combination of base time, t, and prediction time, u.

The results for the fixed base time and fixed prediction window analyses are saved at each iteration. 
An average over the iterations (along with the retarded kernel results) was performed subsequently using the C++ code Average_Split.cpp. 

Details of the models fitted and the analysis performed is given in the manuscript 'Retarded kernels for longitudinal survival analysis and dynamic prediction', Davies, Coolen and Galla (2021). 


The folder Edited_JMbayes_functions contains the edited versions of the function prederrJM (from the JMbayes package) for coxph objects, JMbayes objects and mvJMbayes objects. 

The original codes were copied from prederrJM.coxph.R, prederr.JMbayes.R and prederr.mvJMbayes.R at https://github.com/drizopoulos/JMbayes/tree/master/R.

Edits were made so the calculation of prediction error exactly matches the prediction error equation (Eq. (4.26) in Davies, Galla and Coolen (2021) or, equivelently, the equation for prediction error on pg. 34 of Rizopoulos, D. (2016).  The R package JMbayes for fitting joint models for longitudinal andtime-to-event data using MCMC. Journal of Statistical Software 72(7), 1–46). 

All edits are described in the comments of each code and labelled in line. 
Further details can be found in the Supplementary Material for Davies, Galla and Coolen (2021).
