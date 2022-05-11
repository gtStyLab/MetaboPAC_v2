function [timeVec, concMatrix, fluxMatrix] = processToFittingData(odeData,nTFit,fluxScheme,missingX,fileOutName,fileOutSuffix)
% This function processes data which is to be fit by the various methods
% tested in the LK-DFBA manuscript and supplementary. The main DFE steps
% are performed in this function.
% 
% Written by R.A.Dromms 2016-06-02.

    % Grab ODE data metadata
    nT = length(odeData.timeVec)-1;
    rawTimeVec = odeData.timeVec(:);
    rawConcMatrix = odeData.concMatrix;
    rawFluxMatrix = odeData.fluxMatrix;
        
    tStart = rawTimeVec(1);
    tEnd = rawTimeVec(end);
    numMetabs = size(rawConcMatrix,2);
    
    % Set up flux data appropriately
    if isfield(odeData,'fluxTimeVec') && ~isempty(odeData.fluxTimeVec)
        rawFluxTimeVec = odeData.fluxTimeVec;
    else
        if size(rawFluxMatrix,1) == size(rawConcMatrix,1)
            rawFluxTimeVec = rawTimeVec;
        elseif size(rawFluxMatrix,1) == size(rawConcMatrix,1)
            rawFluxTimeVec = rawTimeVec(1:end-1)+0.5*diff(rawTimeVec(1:2));
        else
            % Automatically just use this case??
            rawFluxTimeVec = linspace(tStart,tEnd,size(rawFluxMatrix,1));
        end
    end

    timeVec = linspace(tStart,tEnd,nTFit+1)';
    fluxTimeVec = timeVec;
    deltaT = (timeVec(end) - timeVec(1))/nT;
            
    if exist('fileOutName','var') && exist('fileOutSuffix','var')
        % Okay, so we're creating a file for diagnostic visualization
        exportFiles = true;
        
        % We'll need this for later
        paramsVec = [];
    else
        exportFiles = false;
    end
    
    if ~exist('missingX','var')
        missingX = false(size(rawConcMatrix,2));
    else
        missingX = logical(missingX);
    end
    
    % We tried several ways of combining data smoothing and flux
    % estimation, some of which we abandoned from further analysis.
    switch fluxScheme
        
        % Scheme 1: Basic interpolation only
        case 1
            
            concMatrix = interp1(rawTimeVec,rawConcMatrix,timeVec,'linear','extrap');
            fluxMatrix = interp1(rawFluxTimeVec,rawFluxMatrix,fluxTimeVec,'linear','extrap');
            
            if exportFiles
                x0 = concMatrix(1,:);
                save(sprintf('%s%s_plotPostInterp.mat',fileOutName,fileOutSuffix),...
                    'concMatrix','fluxMatrix','fluxTimeVec','nT','paramsVec','tEnd','timeVec','tStart','x0');
                clear x0
            end
        
            
        % Scheme 2: Smoothing only
        case 2
            
            numMetabs = size(rawConcMatrix,2);
            impulseParams = nan(7,numMetabs);
            concMatrix = nan(length(timeVec),numMetabs);
            fluxMatrix = interp1(rawFluxTimeVec,rawFluxMatrix,fluxTimeVec,'linear','extrap');
            
            for k = 1:numMetabs
                if all(rawConcMatrix(:,k) == rawConcMatrix(1,k))
                    % Recalculate at timeVec
                    concMatrix(:,k) = rawConcMatrix(1,k)*ones(size(timeVec));
                    % Calculate Pooling Flux
                    fluxMatrix(:,(size(fluxMatrix,2)-numMetabs)+k) = zeros(size(fluxTimeVec));
                else
                    params.bF = 0.3;
                    % 1. Smooth data
                    [impulseParams(:,k), ~,~] = fitSmoothingFcn(rawTimeVec(:),rawConcMatrix(:,k),'impls',params);
                    % 2. Recalculate at timeVec
                    concMatrix(:,k) = evalSmoothingFcn(impulseParams(:,k),timeVec,'impls');
                    % Calculate Pooling Flux
                    fluxMatrix(:,(size(fluxMatrix,2)-numMetabs)+k) = evalSmoothingDeriv(impulseParams(:,k),fluxTimeVec,'impls');
                end
            end
            
            if exportFiles
                x0 = concMatrix(1,:);
                save(sprintf('%s%s_plotPostSmoothed.mat',fileOutName,fileOutSuffix),...
                    'concMatrix','fluxMatrix','fluxTimeVec','nT','paramsVec','tEnd','timeVec','tStart','x0');
                clear x0
            end
            
           
            
            
        % Scheme 3: Smooth + Backslash for fluxes
        case 3
            %impulseParams = nan(7,numMetabs);
            paramSolution = nan(5,numMetabs);
            concMatrix = nan(length(timeVec),numMetabs);
            fluxMatrixTemp = nan(length(fluxTimeVec),size(rawFluxMatrix,2));
            for k = 1:numMetabs
                if ~missingX(k)
                    % Get data smoothing parameters
                    % [impulseParams(:,k), ~,~] = fitSmoothingFcn(rawTimeVec,rawConcMatrix(:,k),'poly');
                    [paramSolution(:,k), ~,~] = fitSmoothingFcn(rawTimeVec,rawConcMatrix(:,k),'rat22');
                    % In principle, better, but only if nTFit = nTData
                    % If nTFit > nTData, then it's much more arguable.
                    concMatrix(:,k) = interp1(rawTimeVec,rawConcMatrix(:,k),timeVec,'linear','extrap');
                    %smoothConcMatrix(:,k) = evalSmoothingFcn(impulseParams(:,k),timeVec,'poly');
                    smoothConcMatrix(:,k) = evalSmoothingFcn(paramSolution(:,k),timeVec,'rat22');
                else
                    % Consistent with dx/dt = 0 assumption.
                    concMatrix(:,k) = odeData.concMatrix(1,k);
                end
            end
            disp(paramSolution)
            
            for k = 1:numMetabs
                if ~missingX(k)
                    % Calculate smoothing function from impulse derivative
                    fluxMatrixTemp(:,(size(fluxMatrixTemp,2)-numMetabs)+k) = evalSmoothingDeriv(paramSolution(:,k),fluxTimeVec,'rat22');
                    %fluxMatrixTemp(:,(size(fluxMatrixTemp,2)-numMetabs)+k) = evalSmoothingDeriv(impulseParams(:,k),fluxTimeVec,'poly');
                else
                    fluxMatrixTemp(:,(size(fluxMatrixTemp,2)-numMetabs)+k) = zeros();
                end
            end
            
            % Replace fixed fluxes with the appropriate values
            fluxMatrixTemp(:,odeData.fixedFluxes) = interp1(rawFluxTimeVec,rawFluxMatrix(:,find(odeData.fixedFluxes)),fluxTimeVec,'linear','extrap');
            
            if exportFiles
                x0 = concMatrix(1,:);
                fluxMatrix = fluxMatrixTemp;
                save(sprintf('%s%s_plotPostSmoothed.mat',fileOutName,fileOutSuffix),...
                    'concMatrix','fluxMatrix','fluxTimeVec','nT','paramsVec','tEnd','timeVec','tStart','x0');
                clear x0 fluxMatrix
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%
%             timePoint = size(concMatrix,1);
%             timeInterval = timeVec(2) - timeVec(1);
%             numMet = size(concMatrix,2);
%             poolingFluxMatrix = zeros(timePoint,5+numMet);
%             for m = 2:1:timePoint-1
%                 poolingFluxMatrix(m,numMet+1:end) = (concMatrix(m+1,:) - concMatrix(m-1,:))/(2*timeInterval);
%             end
%             fluxMatrixTemp = poolingFluxMatrix;
%             fluxMatrixTemp(:,odeData.fixedFluxes) = interp1(rawFluxTimeVec,rawFluxMatrix(:,find(odeData.fixedFluxes)),fluxTimeVec,'linear','extrap');
%             
            %%%%%%%%%%%%%%%%%%%%
            [fluxMatrix] = calcFluxesViaBackslash(fluxMatrixTemp,odeData.STM,odeData.fixedFluxes);
            
            if exportFiles
                x0 = concMatrix(1,:);
                save(sprintf('%s%s_plotPostSmoothedBackslash.mat',fileOutName,fileOutSuffix),...
                    'concMatrix','fluxMatrix','fluxTimeVec','smoothConcMatrix','nT','paramsVec','tEnd','timeVec','tStart','x0');
                clear x0
            end
            
        otherwise

            concMatrix = interp1(rawTimeVec,rawConcMatrix,timeVec,'linear','extrap');
            fluxMatrix = interp1(rawFluxTimeVec,rawFluxMatrix,fluxTimeVec,'linear','extrap');
            
    end

    
end