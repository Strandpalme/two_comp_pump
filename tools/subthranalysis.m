function [subthrFeatureVals] = subthranalysis(V)
%% subthrFeatureVals = subthranalysis(V1) returns a struct array of response features
% This function analyzes sub-threshold response fetaures of the voltage 
% trace V1 and returns a struct array of feature values.
% Input:
%       V  
%           voltage trace of a neuronal response to the stimulus used in
%           Meiser2019.
%
% Output:
%       responseFetaures    
%           struct-array that holds a response feature value for the
%           corresponding field. 
%           Fields are: 
%           .restpot    resting potential
%           .inpres     input resistance
%           .vlotsag    voltage sag
%           .spontact   spontanious activity


% written by Kevin Sandbote
% march 2021
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
subthrFeatureVals.restpot   = restpot;
subthrFeatureVals.inpres    = inpres;
subthrFeatureVals.voltsag   = voltsag;
subthrFeatureVals.spontact  = spontact;
