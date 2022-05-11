function steadyStatePenalty = pen_devSSdist_UDreg(fluxMatrix)
    steadyStateProfiles = load('steadyStateProfiles_UDreg.mat');
    steadyStateFluxes = steadyStateProfiles.ssFluxProfile_UDreg;
    
    [fluxnT,numFlux] = size(fluxMatrix);
    
    steadyStateFluxMatrix = repmat(steadyStateFluxes,[fluxnT,1]);
    
    %gradually increased weighting towards end of time range
    fractionCount = 0.25;
    nT_count = round(fluxnT * fractionCount);
    nT_noCount = fluxnT - nT_count;
    stepSize = 1/(nT_count-1);
    weightVector = 0:stepSize:1;
    weightMatrix = repmat(weightVector',[1,numFlux]);
    weightedFluxDeviation = (weightMatrix .* ((steadyStateFluxMatrix(nT_noCount+1:end,1:numFlux)-fluxMatrix(nT_noCount+1:end,1:numFlux))./fluxMatrix(nT_noCount+1:end,1:numFlux)).^2);
    weightedFluxPenalty = sum(sum(weightedFluxDeviation));
  
    steadyStatePenalty = weightedFluxPenalty;






end