function [penalty,massbalance_penalty,conc_penalty,oneContMetCorr_penalty,oneContMetCurveFit_penalty,BST_penalty,Vcalc_pool,Vcalc] = calcPenalty_Determined(testRF,modelInfo,relative_Vpool,numMetabs,numFlux,MA_reactions,relative_concMatrix,timeVec,fluxTimeVec,RF_kinetics,knownMet,knownKinetics,rep)
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
    newFluxMatrixTemp = calcFluxWithKinetics_Determined(absolute_concMatrix,timeVec,knownKinetics);
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
        oneContMetCorr_penalty = pen_oneContMetCorr_Determined(absolute_concMatrix,Vcalc,timeVec,fluxTimeVec);

        % Calculate curve fit penalty for reactions controlled by one
        % metabolite
        oneContMetCurveFit_penalty = pen_oneContMetCurveFit_Determined(absolute_concMatrix,Vcalc,timeVec,fluxTimeVec);

        % Calculate BST fit penalty
        nT = size(absolute_concMatrix,1)-1;
        BST_penalty = pen_BSTfit_Determined(nT,absolute_concMatrix,Vcalc,MA_reactions);
        
        % Calculate total penalty
        penalty = abs(massbalance_penalty) + abs(conc_penalty) + abs(oneContMetCorr_penalty) + abs(oneContMetCurveFit_penalty) + abs(BST_penalty);
        
        if isnan(penalty)
            penalty = 1e7;
        end
    else
        penalty = sum(1./RF)*1e7;
    end
end