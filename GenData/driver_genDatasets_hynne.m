%% Generate hi-Res datasets
clc; clear; close all;

dataDir = 'data';

nTotal = 1;

for k = 1:nTotal
    
    datasetNames{k} = wrapper_genOdeData_hynne(k,dataDir);
    
end

%% Generate my noisy datasets, over various nT and CoV combinations
nTList = [15 50];
covList = [0.05 0.15];
numSets = 3;

% a. Get the hi-res dataset file name
hiResDataFileName = ['data/hynne_k-01_hiRes.mat'];

% b. Loop nT values
for nT = nTList

    % c. Loop CoV Values
    for cov = covList
        % Generate my noisy datasets
        wrapper_genNoisyData(hiResDataFileName,nT,cov,numSets)
    end
end
    
