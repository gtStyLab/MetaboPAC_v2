clear
clc

% This script aims to sample 10,000 sets of RFs and their response factors
% for troubleshooting and penalty weight optimization purposes 

% Set up initial information about system
modelInfo_generic
load('generic_data.mat');
modelInfo.fixedFluxes = (modelInfo.vBounds(:,1)== modelInfo.vBounds(:,2));

numMetabs = size(modelInfo.S,1);
numFlux = size(modelInfo.S,2);

% Find Mass Action reactions
[row_MA col_MA] = find(modelInfo.S < 0);
MA_reactions = [row_MA col_MA];

% Load relative abundances
load('generic_trueRF');
trueRF = trueRF(1,:);

% Calculate relative abundances
relative_concMatrix = concMatrix.*trueRF;

% Calculate relative pooling fluxes by calculating change in relative
% abundances and dividing by change in time
for i = 2:size(relative_concMatrix,1)
    relative_Vpool(i-1,:) = (relative_concMatrix(i,:)-relative_concMatrix(i-1,:))./(timeVec(i)-timeVec(i-1));
end

%focus on optimization approach
knownMet = [];
RF_kinetics = trueRF(knownMet);
percKnownKinetics = 0;
knownKinetics = [];

maxRandVal = 1000;
minRandVal = 1;

gs_RF = trueRF(1,setdiff(1:numMetabs,knownMet));
%find penalty for true RF (gold standard) 
[penalty_gs,massbalance_penalty_gs,conc_penalty_gs,oneContMetCorr_penalty_gs,oneContMetCurveFit_penalty_gs,BST_penalty_gs,~,~] = calcPenalty_generic(gs_RF,modelInfo,relative_Vpool,numMetabs,...
        numFlux,MA_reactions,relative_concMatrix,timeVec,fluxTimeVec,RF_kinetics,...
        knownMet,knownKinetics);
penaltyTable_gs = table(trueRF,massbalance_penalty_gs,conc_penalty_gs,...
    oneContMetCorr_penalty_gs,oneContMetCurveFit_penalty_gs,BST_penalty_gs,penalty_gs);  
    
%find penalty for 10,000 randomly sampled RF
numSampling = 10000;
testRF_all = nan(numSampling,numMetabs-length(knownMet));
penalty_all = nan(numSampling,1);
mbPenalty_all = nan(numSampling,1);
conc_penalty_all = nan(numSampling,1);
oneCorrPen_all = nan(numSampling,1);
oneCurve_all = nan(numSampling,1);
BSTpen_all = nan(numSampling,1);

for rep = 1:10000
    rng(rep)
    testRF =  (maxRandVal-minRandVal).*rand(1,numMetabs-length(knownMet)) + minRandVal;

    [penalty,massbalance_penalty,conc_penalty,oneContMetCorr_penalty,oneContMetCurveFit_penalty,BST_penalty] ...
        = calcPenalty_generic(testRF,modelInfo,relative_Vpool,numMetabs,...
            numFlux,MA_reactions,relative_concMatrix,timeVec,fluxTimeVec,RF_kinetics,...
            knownMet,knownKinetics);

    testRF_all(rep,:) = testRF;
    penalty_all(rep) = penalty; 
    mbPenalty_all(rep) = massbalance_penalty;
    conc_penalty_all(rep) = conc_penalty;
    oneCorrPen_all(rep) = oneContMetCorr_penalty;
    oneCurve_all(rep) = oneContMetCurveFit_penalty;
    BSTpen_all(rep) = BST_penalty;

end

[sortedPenaltyAll,sortIdx] = sort(penalty_all,'ascend');
sortedTestRFAll = testRF_all(sortIdx,:);
sortedMBPen = mbPenalty_all(sortIdx);
sortedConcPen = conc_penalty_all(sortIdx);
sortedOneCorrPen = oneCorrPen_all(sortIdx);
sortedOneCurvePen = oneCurve_all(sortIdx);
sortedBstPen = BSTpen_all(sortIdx);

penaltyTable = table(sortedTestRFAll,sortedMBPen,sortedConcPen,sortedOneCorrPen,...
    sortedOneCurvePen,sortedBstPen,sortedPenaltyAll);

save('generic_sampling_final.mat');