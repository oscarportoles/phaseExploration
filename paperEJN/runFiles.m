% This file runs all the files needed to reproduce all analysis in the
% paper when you already have the utput of the HSMM-MVPA. The HSMM-MVPA
% must be run before with the code probided along with the paper "The
% discovery of Prcessing Stages: Extension of Sternberg's Method, Anderson
% et al. 2016". It starts reading the events structure of the data
% provaided by Anderson & Borst 2015 in "The discovery of processing stages: 
% Analyzing EEG data with hidden semi-Markov models".

run TXTevent2EEGLABstruc.m
run downSampleFilterBands.m
run epochBaseLine.m
run epochPostStimuli.m
run powerUitpc.m
run dwPLIconnect.m
run visualizeConnetivityFigures.m