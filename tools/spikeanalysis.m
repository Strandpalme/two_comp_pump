function [spikeFeatureVals] = spikeanalysis(V)
%% spikeFeatureVals = spikeanalysis(V1) returns a struct array of response features
% This function analyzes spike response fetaures of the voltage trace V1 
% and returns a struct array of feature values.
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
%           .sc1nA      spike count 1 nA
%           .amp3rd     amplitude of 3rd spike
%           .diff2to4   amp. diff. between spike 2 and 4
%           .sc05nA     spike count 0.5 nA
%           .sc15nA     spike count 1.5 nA

% written by Kevin Sandbote
% march 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% defining regions of interest
% give indices for time windows and determine response feature values

% indices of certain current stepps: IDX.currentStep = [start, stop]
IDXsc1nA    = [1900000, 1950000]; % depol. current step 1 nA
IDXsc05nA   = [1500000, 1550000]; % depol. current step 0.5 nA
IDXsc15nA   = [ 700000,  740000]; % depol. current step 1.5 nA

% assigning time windows 
curstep1    = V(IDXsc1nA(1)  : IDXsc1nA(2));
curstep05   = V(IDXsc05nA(1) : IDXsc05nA(2));
curstep15   = V(IDXsc15nA(1) : IDXsc15nA(2));


%% calculating response feature values
% detecting peak values in response signal to 1 nA delpol. current
peak1nA = findpeaks(curstep1,'MinPeakProminence',4);
% counting spikes
sc1nA = length(peak1nA);

% calculating amplitude of 3rd spike and difference between 2nd and 4th
% spike amplitude

if sc1nA >= 4 
    % case 1: more than 3 spikes. Spike features can be calculated.
    % Starting with detecting peak values of after hyperpolarization
    negpeak1nA =  -findpeaks(-curstep1,'MinPeakProminence',1);
    
    if length(negpeak1nA) >= 4
        % if there are more than 3 after hyperpolarizations, the amplitudes
        % can be calculated as follows
        amplitude1nA = peak1nA(1:4) - negpeak1nA(1:4);
        amp3rd = amplitude1nA(3);
        diff2to4 = amplitude1nA(2) - amplitude1nA(4);
    elseif length(negpeak1nA) == 3
        % if there are only 3 after hyperpolarizations, amplitude of the 
        % 3rd spike can be calculated but not the difference between the 
        % 2nd and 4th spike amplitude
        amplitude1nA = peak1nA(1:3) - negpeak1nA(1:3);
        amp3rd = amplitude1nA(3);
        diff2to4 = NaN;
    else
        % in case there are less than 3 after hyperpolarizatons, no ap
        % feature can be calculated
        amp3rd = NaN;
        diff2to4 = NaN;
    end
    
elseif nSpikes1nA == 3
    % when three spikes are detected, amplitude of the 3rd spike can be
    % calculated, but not the difference...
    negpeak1nA =  -findpeaks(-curstep1,'MinPeakProminence',1);
    
    if length(negpeak1nA) >= 3
        % three or more afterhyperpolarizations found
        amplitude1nA = peak1nA(1:3) - negpeak1nA(1:3);
        amp3rd = amplitude1nA(3);
        diff2to4 = NaN;
    else
        % less than three a.hps. found, nothing can be calculated
        amp3rd = NaN;
        diff2to4 = NaN;
    end
    
else
    % less than 3 spikes detected, no features can be calculated
    amp3rd = NaN;
    diff2to4 = NaN;
end    

% detcting peaks for 0.5 nA stimulation
peak05nA = findpeaks(curstep05,'MinPeakProminence',4);
% calculating spike count
sc05nA = length(peak05nA);

% detecting and counting spikes for 1.5 nA
peak15nA = findpeaks(curstep15,'MinPeakProminence',4);
sc15nA = length(peak15nA);


%% assining return values
spikeFeatureVals.sc1nA    = sc1nA;
spikeFeatureVals.amp3rd   = amp3rd;
spikeFeatureVals.diff2to4 = diff2to4;
spikeFeatureVals.sc05nA   = sc05nA;
spikeFeatureVals.sc15nA   = sc15nA;
