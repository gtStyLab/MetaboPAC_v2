function [absolute_concMatrix predicted_responseFactors] = MetaboPAC_UDreg(approach_option,rep,percKnownKinetics,rand_idx)

% Set up initial information about system
modelInfo_UDreg
load('UDreg_k-01_hiRes.mat');

modelInfo.fixedFluxes = (modelInfo.vBounds(:,1)== modelInfo.vBounds(:,2));

numMetabs = size(modelInfo.S,1);
numFlux = size(modelInfo.S,2);

numExtraMets = numMetabs - size(concMatrix,2);
if numExtraMets > 0
    concMatrix(:,end+1:end+numExtraMets) = zeros(size(concMatrix,1),numExtraMets);
end

% Find Mass Action reactions
[row_MA col_MA] = find(modelInfo.S < 0);
MA_reactions = [row_MA col_MA];

% Load relative abundances
load('UDreg_trueRF');
trueRF = trueRF(rep,:);

% Calculate relative abundances
relative_concMatrix = concMatrix.*trueRF;

% Perform kinetic equations approach
knownKinetics_fixed = [4,5,7,8];
remKinetics = setdiff(1:numFlux,knownKinetics_fixed);
numKnownKinetics = round(percKnownKinetics/100*length(remKinetics));
rng(rand_idx)
knownKinetics_temp = remKinetics(randperm(length(remKinetics),numKnownKinetics));
knownKinetics = [knownKinetics_fixed,knownKinetics_temp];

switch approach_option
    case 1
        % KE
        approach_string = 'KE_';
        for i = 1:48
            [KEapproach_results(i,:) knownMet] = solveRF_KEapproach_UDreg(relative_concMatrix,timeVec,knownKinetics);
        end
        if knownMet ~= 0
            RF_kinetics = median(KEapproach_results);
        else
            RF_kinetics = [];
            knownMet = [];
        end
    case 2
        % opt only
        approach_string = '';
        RF_kinetics = [];
        knownMet = [];
end

% Set up information for optimization approach
maxRandVal = 1000;
minRandVal = 1;
lb = minRandVal*ones(1,numMetabs-length(knownMet));
ub = maxRandVal*ones(1,numMetabs-length(knownMet));
options.ConstraintTolerance	= 1e10;
options.MaxFunctionEvaluations = 5000;

% Calculate relative pooling fluxes by calculating change in relative
% abundances and dividing by change in time
for i = 2:size(relative_concMatrix,1)
    relative_Vpool(i-1,:) = (relative_concMatrix(i,:)-relative_concMatrix(i-1,:))./(timeVec(i)-timeVec(i-1));
end

% Perform optimization approach
fval = 1500*ones(48,1);
RF_opt = ones(48,numMetabs);
if length(knownMet) ~= numMetabs
    for i = 1:48
        % Set initial seed for optimizier
        rng(i)
        x0 = (maxRandVal-minRandVal).*rand(1,numMetabs-length(knownMet)) + minRandVal;
        ga_options = optimoptions('ga','MaxGenerations',100,'InitialPopulationMatrix',x0);

        % Perform genetic algorithm optimization
        [optimalRF fval(i,1)] = ga(@(testRF) calcPenalty(testRF,modelInfo,relative_Vpool,numMetabs,numFlux,MA_reactions,relative_concMatrix,timeVec,fluxTimeVec,RF_kinetics,knownMet,knownKinetics,rep),numMetabs-length(knownMet),[],[],[],[],lb,ub,[],ga_options);

        RF_temp = zeros(1,numMetabs);
        RF_temp(1,knownMet) = RF_kinetics;
        RF_temp(1,setdiff(1:numMetabs,knownMet)) = optimalRF;

        RF_opt(i,:) = RF_temp;
    end
else
    % If optimization approach not performed, use results from kinetic
    % equations approach
    RF_opt = KEapproach_results;
end

predicted_responseFactors = median(RF_opt);
absolute_concMatrix = relative_concMatrix./predicted_responseFactors;

save(sprintf('results/UDreg_MetaboPAC_%saddKinetics-%03d_rep-%03d_rand-%03d.mat',approach_string,percKnownKinetics,rep,rand_idx));
end



function penalty = calcPenalty(testRF,modelInfo,relative_Vpool,numMetabs,numFlux,MA_reactions,relative_concMatrix,timeVec,fluxTimeVec,RF_kinetics,knownMet,knownKinetics,rep)
    % Consolidate response factor values
    RF = zeros(1,numMetabs);
    RF(1,knownMet) = RF_kinetics;
    RF(1,setdiff(1:numMetabs,knownMet)) = testRF;

    % Calculate change in inferred absolute concentrations over time (i.e.
    % pooling fluxes, Vpool)
    test_Vpool(:,1:numMetabs) = relative_Vpool(:,1:numMetabs)./RF;
    
    % Calculate inferred absolute concentrations
    absolute_concMatrix(:,1:numMetabs) = relative_concMatrix(:,1:numMetabs)./RF;
    
    % Set up flux matrix before inferring fluxes
    newFluxMatrixTemp = calcFluxWithKinetics_UDReg(absolute_concMatrix,timeVec,knownKinetics);
    modelInfo.fixedFluxes = ~isnan(newFluxMatrixTemp(1,1:numFlux))';
    newFluxMatrixTemp(:,numFlux+1:end) = test_Vpool;
    
    % Infer fluxes using pinv
    Vcalc = calcFluxesViaPinv(newFluxMatrixTemp,modelInfo.S,modelInfo.fixedFluxes);    
    
    % Calculate penalties
    if sum(sum(Vcalc))~=0 && ~isnan(sum(sum(Vcalc)))
        % Calculate pooling flux penalty
        for i = 1:size(Vcalc,1)
            Vcalc_pool(i,:) = modelInfo.S*Vcalc(i,1:numFlux)';
        end

        % Calculate mass balance penalty
        massbalance_penalty = sqrt(sum(sum((Vcalc_pool(2:end,:) - Vcalc(2:end,numFlux+1:end)).^2)));
        %rmse = sqrt(mean((Vcalc_pool(2:end,:) - Vcalc(2:end,numFlux+1:end)).^2,'all'));
        %nrmse = sqrt(mean(((Vcalc_pool(2:end,:) - Vcalc(2:end,numFlux+1:end))./range(Vcalc(2:end,numFlux+1:end))).^2,'all'));
        
        % Calculate max. concentration penalty
        if any(any(absolute_concMatrix(:,1:numMetabs-2) > 50))
            conc_penalty = max(max(absolute_concMatrix(:,1:numMetabs-2)));
        else
            conc_penalty = 0;
        end
        

        % Calculate correlation penalty for reactions controlled by one
        % metabolite
        oneContMetCorr_penalty = pen_oneContMetCorr_UDreg(absolute_concMatrix,Vcalc,timeVec,fluxTimeVec);

        % Calculate curve fit penalty for reactions controlled by one
        % metabolite
        oneContMetCurveFit_penalty = pen_oneContMetCurveFit_UDreg(absolute_concMatrix,Vcalc,timeVec,fluxTimeVec);

        % Calculate BST fit penalty
        nT = size(absolute_concMatrix,1)-1;
        BST_penalty = pen_BSTfit_UDreg(nT,absolute_concMatrix,Vcalc);
        
        %Calculate steady state penalty
        ss_penalty = pen_devSSdist_UDreg(Vcalc(:,1:numFlux));
        
        % Calculate total penalty
        penalty = 1000 * abs(massbalance_penalty) + 10 * abs(conc_penalty) + 10 * abs(oneContMetCorr_penalty) + 10 * abs(oneContMetCurveFit_penalty) + 10 * abs(BST_penalty) + 10 * ss_penalty;

        if isnan(penalty)
            penalty = 1e7;
        end
    else
        penalty = sum(1./RF)*1e7;
    end
end
