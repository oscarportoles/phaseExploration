Code supporting the paper: 
“Characterizing synchrony patterns across cognitive task stages of associative recognition memory”
Authors: Oscar Portoles, Jelmer P. Borst, Marieke van Vugt.
Contact: o.portoles.marin @ rug.nl
Date: July 2017

This document along with the paper provides all information to replicate the results of the paper. 

The data can be found in the paper “The discovery of processing stages: Analyzing EEG data with hidden semi-Markov models; J.P. Borst & J.R. Anderson; 2015.

The code and a manual to apply the hidden semi-Markov model multivariate pattern analysis (HSMM-MVPA) can be found in the paper: The discovery of processing stages: Extension of Sternberg’s method; J.R. Anderson et al.; 2016. From this code the variable “eventprobs” is needed. “eventprobs” contains the probability distributions of the location of stage’s onsets for all trials and subjects  concatenated together. Estimating “eventprobs” took about 14 hours in the the Peregrine cluster (http://www.rug.nl/society-business/centre-for-information-technology/research/services/hpc/facilities/peregrine-hpc-cluster?lang=en)

Once you have “eventprobs” you can run the Matlab file “runFiles.m”. This file will run all files in correct order to produce the figures of the paper from the original data. You may need to change the paths by default in all Matlab files to the path where you have the data. Each Matlab file will produce a data set that will be used by the next Matlab file. Running all files might several hours depending on your computer.

PlotPaperEJN.m will allow you to produce figures two and four.

