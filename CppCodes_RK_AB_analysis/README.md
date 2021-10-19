# C++ Codes to perform Data Analysis of Retarded Kernel models A and B in Davies, Coolen and Galla (2021)

## Model A

The folder ModelA contains three C++ codes
1. AIDS_RK_ModelA.cpp
2. LIVER_RK_ModelA.cpp
3. PBC_RK_ModelA.cpp

Each code reads in 20 training and test data sets created by the R code in the folder RCodes_JM_LM_analysis/DataAnalysis corresponding to the relevant data set (AIDS, Liver, PBC). At each iteration, the code performs the fixed base time and fixed prediction window analysis for retarded kernel Model A for these 20 random splits of the data.

Models are fitted to the training data while prediction error is calculated for the test data for each combination of base time, t, and prediction time, u.

The results for the fixed base time and fixed prediction window analyses are saved at each iteration. An average over the iterations (along with the retarded kernel results) was performed subsequently using the C++ code Average_Split.cpp.

Details of the models fitted and the analysis performed is given in the manuscript 'Retarded kernels for longitudinal survival analysis and dynamic prediction', Davies, Coolen and Galla (2021).

## Model B 

The folder ModelB contains three C++ codes
1. AIDS_RK_ModelB.cpp
2. LIVER_RK_ModelB.cpp
3. PBC_RK_ModelB.cpp

Each code reads in 20 training and test data sets created by the R code in the folder RCodes_JM_LM_analysis/DataAnalysis corresponding to the relevant data set (AIDS, Liver, PBC). At each iteration, the code performs the fixed base time and fixed prediction window analysis for retarded kernel Model B for these 20 random splits of the data.

Models are fitted to the training data while prediction error is calculated for the test data for each combination of base time, t, and prediction time, u.

The results for the fixed base time and fixed prediction window analyses are saved at each iteration. An average over the iterations (along with the retarded kernel results) was performed subsequently using the C++ code Average_Split.cpp.

Details of the models fitted and the analysis performed is given in the manuscript 'Retarded kernels for longitudinal survival analysis and dynamic prediction', Davies, Coolen and Galla (2021).

This folder also contains a pdf 'ModelB.pdf'. This contains detail about the Model B equation quoted in the supplementary material and the Model B equation used to write the codes. This is provided for the interested reader who wishes to better understand the codes - otherwise it can be ignored as it contains no new information.

## Average_Split.cpp

A simple C++ code that reads in the vectors of prediction error for each model (fixed base time and fixed prediction window) at each iteration and averages over the iterations.
