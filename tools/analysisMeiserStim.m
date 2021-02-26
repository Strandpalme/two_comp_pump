function [responseFeatures] = analysisMeiserStim(V1)
%% RF = analysisMeiserStim(V1) returns a struct array of response features
% This function analyzes response fetaures of the voltage trace V1 and
% returns a struct array of feature values.
% Input:
%       V1  
%           voltage trace of a neuronal response to the stimulus used in
%           Meiser2019.
%
% Output:
%       responseFetaures    
%           struct-array that holds a response feature value for the
%           corresponding field. 
%           Fields are: 
%           resting potential
%           input resistance
%           delta R
%           resting potential t20
%           input resistance t20
%           spontanious activity
%           spike count 1 nA
%           spike height 3rd
%           amp diff 2-4
%           voltage sag
%           spike count 0.5 nA
%           spike count 1.5 nA
%
% written by Kevin Sandbote
% march 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% defining regions of interest
