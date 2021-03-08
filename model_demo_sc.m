%% Demo script for two-comp-model
% calculation and analysis of a modeled T-cell response
%
% written by Kevin Sandbote
% march 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% adding the model to the matlab path
addpath(genpath('D:\Projekte\two_comp_pump'))


%% generating stimulus

dt = 0.01; % [ms] time steps

% original stimulus used in Meiser 2019 in 0.01 ms bins
load('D:\Projekte\two_comp_pump\data\originalstimulus_binsize_0p01ms.mat')
stimulus = stim_Tcharacteristics .* 1000; % [pA]
nTrials = 20;
Iinj = repmat(stimulus,nTrials);

%% pump parameters
% assigning pump parameters to a 1x3 array to test different parameter
% values
P(1) = 800e-6;          % [uA] max pump current (=800[pA])
P(2) = 10 * 0.06e-6;    % conversion factor from INa into [Na] [mM/pA.ms]
P(3) = 2 * 0.06e-6;     % conversion factor from Ipump into [Na] [mM/pA.ms]

%% features analyzed
% generating two cell arrays containing the name of features analyzed by
% "subthranalysis.m" and "spikeanalysis.m" in correct order and with units.

subthrFeatNames = {             
    'resting potential',    'mV'          
    'input resistance',     'Mohm'
    'voltage sag',          'mV'   
    'spont. activity',      'no. spikes'};
    
spikeFeatNames = {
    'Spike count 1 nA',                                     'no. spikes'
    'amplitude of 3rd spike',                               'mV'
    'difference between amplitudes of 2nd and 4th spike'    'mV'
    'spike count 0.5 nA',                                   'no. spikes'
    'spike count 1.5 nA',                                   'no. spikes'};

%% call model
tic
[V1, V2] = TcellDoublePumpFitting(Iinj, dt, P, 0);
toc


%% analyse response
% reorganizing the voltage trace into a trials x time matrix
datamat = reshape(V1,[],nTrials)';

% preallocating feature arrays
subthrFeat = zeros([nTrials,4]);
spikeFeat = zeros([nTrials,5]);
% analyzing trial by trial
for k = 1 : nTrials
    subthrFeat(k,:) = subthranalysis(datamat(k,:));
    spikeFeat(k,:) = spikeanalysis(datamat(k,:));
end


%% presenting results
figure(1)
for k = 1:size(subthrFeat,2)
    subplot(3,2,k)
    plot(subthrFeat(:,k))
    % labeling axes
    xlabel(subthrFeatNames(k,1))
    ylabel(subthrFeatNames(k,2))
end

figure
for k = 1:size(spikeFeat,2)
    subplot(3,2,k)
    plot(spikeFeat(:,k))
    % labeling axes
    xlabel(spikeFeatNames(k,1))
    ylabel(spikeFeatNames(k,2))
end
