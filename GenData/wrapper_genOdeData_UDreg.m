function datasetFileName = wrapper_genOdeData_UDreg(paramSetNum,dataDir)
% The job of this function is to generate noiseless, hi-res data from the
% ODE model. For our noisy cases, we'll feed the results from this script
% into a second function to generate numReps noisy data at a given
% combination of CoV and lo-res nT values. 
    
    % Eh let's just hardcode the stem for this stuff here
    datasetFileName = sprintf('UDreg_k-%02d_hiRes.mat',paramSetNum);

    % We want hi-res version of our time-courses available:
    nT = 50;

    % For this function, we're going to just hard code stuff
    tStart = 0;
    tEnd = 10;
    numMetabs = 4;
    numFluxes = 12;
    
    % Tack on target directory, if we specified a destination
    if exist('dataDir','var')
        if ~exist(dataDir,'dir')
            mkdir(dataDir);
        end
        datasetFileName = sprintf('%s/%s',dataDir,datasetFileName);
    end

    % Define ODE Test Parameter sets -- we can add x0 values to this,
    % perhaps? But not doing so right now. Instead, we have:
    try
        paramsVec = getParamsVecNum_UDreg(paramSetNum);
        paramErrorCatch = false;
    catch
        paramErrorCatch = true;
    end


    % Check if we have a valid parameter set to work from before we move on
    if ~paramErrorCatch
    
        % Initialize our structure
        odeDataset.tStart = tStart;
        odeDataset.tEnd = tEnd;
        odeDataset.nT = nT;
        odeDataset.paramsVec = paramsVec(1:18);
        odeDataset.x0 = paramsVec(19:22);

        % Generate the actual ODE time course
% %         [odeDataset.timeVec, odeDataset.concMatrix, odeDataset.fluxMatrix] = solveODE(odeDataset.tStart,odeDataset.tEnd,odeDataset.nT,odeDataset.x0,odeDataset.paramsVec);
        [odeDataset.timeVec, odeDataset.concMatrix, odeDataset.fluxMatrix] = solveOdeMmUDreg(odeDataset.tStart,odeDataset.tEnd,odeDataset.nT,odeDataset.x0,odeDataset.paramsVec);
        odeDataset.fluxMatrix(end,:) = [];
        odeDataset.fluxTimeVec = odeDataset.timeVec(1:end-1)+0.5*diff(odeDataset.timeVec(1:2));

        % Let's actually export this now
        save(datasetFileName,'-struct','odeDataset');

        % For reference, in case it turns out I want the '-struct' part later
% %         save(sprintf('odeData_k%d.mat',k),'-struct','singleDataStruct')
        
    else
        warning('Error assigning parameters, due to invalid paramSetNum = %d\nSkipping this dataset.',paramSetNum)
    end

end