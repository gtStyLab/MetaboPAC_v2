function wrapper_genNoisyData(hiResDataFileName,nT,cov,numSets)

    % 1) Load in hiResDataFileName
    hiResData = load(hiResDataFileName);
    
    % Second '_' char should be the one at 'k-%02d_hiRes'
    idxEnd = strfind(hiResDataFileName,'_');
    loResDataFileNameStem = sprintf('%s_nT-%03d_cov-%02d_rep-',hiResDataFileName(1:idxEnd(2)-1),nT,cov*100);
    
    
    % 2) Interpolate to lo-res nT sampling
    timeVec = linspace(hiResData.tStart,hiResData.tEnd,nT+1)';
    fluxTimeVec = timeVec(1:end-1)+0.5*diff(timeVec(1:2))';
    %fluxTimeVec = timeVec;

    loResConcMatrix = interp1(hiResData.timeVec,hiResData.concMatrix,timeVec,'linear','extrap');
    loResFluxMatrix = interp1(hiResData.fluxTimeVec,hiResData.fluxMatrix,fluxTimeVec,'linear','extrap');
        
    
    loResData.tStart = hiResData.tStart;
    loResData.tEnd = hiResData.tEnd;
    loResData.nT = nT;
    loResData.paramsVec = hiResData.paramsVec;
    loResData.x0 = hiResData.x0;
    loResData.timeVec = timeVec;
    loResData.fluxTimeVec = fluxTimeVec;
    
    % 3) Loop through for numSets noisy datasets
    for k = 1:numSets
        
        % 3a) Add in noise: use noiseless data + random * cov 
        loResData.concMatrix = loResConcMatrix + loResConcMatrix .* (cov*randn(size(loResConcMatrix)));
        loResData.fluxMatrix = loResFluxMatrix + loResFluxMatrix .* (cov*randn(size(loResFluxMatrix)));
        
        % Guarantee we have no *negative* concentation values
        loResData.concMatrix(loResData.concMatrix<0) = 0;
        
        % First data point is always 'correct' and noiseless
        loResData.concMatrix(1,:) = loResData.x0;
        
        % v1 is fixed value, no noise added
        loResData.fluxMatrix(:,1) = loResFluxMatrix(:,1);
        
        % 3b) Save out this noisy dataset
        save(sprintf('%s%03d.mat',loResDataFileNameStem,k),'-struct','loResData')
        
    end
end