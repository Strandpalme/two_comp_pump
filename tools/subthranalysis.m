function [subthrFeatureVals] = subthranalysis(V)
%% subthrFeatureVals = subthranalysis(V1) returns a struct array of response features
% This function analyzes sub-threshold response fetaures of the voltage 
% trace V1 and returns a struct array of feature values.
% written by Kevin Sandbote
% march 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +++Input:
%       V  
%           voltage trace of a neuronal response to the stimulus used in
%           Meiser2019.
%
% +++Output:
%       responseFetaures    
%           array that holds a response feature value for following 
%           features: 
%           1. restpot    resting potential
%           2. inpres     input resistance
%           3. vlotsag    voltage sag
%           4. spontact   spontanious activity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% defining regions of interest
% give indices for time windows and determine response feature values

% indices of certain current stepps: IDX.currentStep = [start, stop]
IDXrest     = [ 250000,  300000]; % resting potential
IDXhyper    = [1100000, 1157000]; % hyperpolarization with -1 nA
IDXspont    = [2800000, 2850000]; % no stimulus period for spont. act. det.

% calc resting potential
restpot = median(V(IDXrest(1) : IDXrest(2)));

% calc input resistance
hyperpot = median(V(IDXhyper(1) : IDXhyper(2)));
inpres = restpot - hyperpot;

% calc voltage sag
voltsag = hyperpot - min(V(IDXhyper(1) : IDXhyper(2)));

% detecting sponaneous activity

spontact = length(findpeaks(V(IDXspont(1) : IDXspont(2)),...
           'MinPeakProminence',4));


%% assining values to output struct
subthrFeatureVals(1)   = restpot;
subthrFeatureVals(2)   = inpres;
subthrFeatureVals(3)   = voltsag;
subthrFeatureVals(4)   = spontact;
