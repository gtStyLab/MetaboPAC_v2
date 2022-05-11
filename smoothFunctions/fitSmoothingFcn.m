function [paramSolution, rmsd, reportStruct] = fitSmoothingFcn(tList,yList,funcType,params)
% Fits a smoothing function to the functional form specified in FUNCTYPE.
% TLIST is the independent variable, and YLIST is the dependent variable.
% Returns the parameter list and the RMSD of the fit for that function.
    
    % Initialize the output to nothing useful
    rmsd = inf;
    paramSolution = NaN(7,1);
    reportStruct = [];
    
    % Set defaults
    ub = tList(end)+diff(tList(end-1:end));
    lb = tList(1)-diff(tList(1:2));
    hFactor = 0.1;
    betaFactor = 2;
    
    TolFun = 1e-6*min(abs(yList));
    TolX = 1e-6*min(abs(yList));
    MaxFunEvals = 1e5;
    MaxIter = 1e5;
    algorithmName = 'interior-point';
    maxVal = 1e5;
    numTries = 20;

    % Allow arbitrary UB and LB, or set a buffer of deltaT
    if exist('params','var')
        if isfield(params,'ub')
            ub = params.ub;
        end

        if isfield(params,'ub')
            lb = params.lb;
        end
        
        if isfield(params,'hF')
            hFactor = params.hF;
        end

        if isfield(params,'bF')
            betaFactor = params.bF/min(diff(tList));
        end

        if isfield(params,'maxVal')
            maxVal = params.maxVal;
        end

        if isfield(params,'TolFun')
            TolFun = params.TolFun;
        end

        if isfield(params,'TolX')
            TolX = params.TolX;
        end

        if isfield(params,'MaxFunEvals')
            MaxFunEvals = params.MaxFunEvals;
        end

        if isfield(params,'MaxIter')
            MaxIter = params.MaxIter;
        end

        if isfield(params,'numTries')
            numTries = params.numTries;
        end

        if isfield(params,'algorithm')
            if params.algorithm == 2
                algorithmName = 'active-set';
            elseif params.algorithm == 3
                algorithmName = 'sqp';
            elseif params.algorithm == 4
                algorithmName = 'trust-region-reflective';
            else
                algorithmName = 'interior-point';
            end
        end
    end
        
    solverOptions = optimset('Display','off','algorithm',algorithmName,'TolFun',TolFun,'TolX',TolX,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter);

    
    if strcmp(funcType(1:4),'poly')
        
        n = str2double(funcType(5));
        P = polyfit(tList,yList,n);
        paramSolution(1:length(P)) = P';
        rmsd = sqrt(sum((polyval(P,tList) - yList).^2)/(length(tList)-length(P)));
        
        
    elseif strcmp(funcType,'rat11')
    
        % Set up for the solver
        fitFcn = @(P,t) ((P(1)*t + P(2))./(t + P(3)));
        fMinFcn = @(P) sum((fitFcn(P,tList) - yList).^2);
        numParams = 3;
        
        success = 0;
        while ~success
            try
                [pInit,~,~,~,~] = nlinfit(tList,yList,fitFcn,rand(numParams,1));
                success = 1;
            end
        end

        tryK = 1;
        
        % We're going to run this multiple times to make sure we get
        % something reasonable
        for i = 1:numTries
            clear fOpt1 fOpt2 R1 R2

            % Fit to the data
            success = 0;
            while ~success
                try
                    pInit1 = pInit + exp(round(log10(abs(pInit))+2)).*rand(size(pInit));
                    pInit1(numParams) = -ub - exp(round(log10(abs(ub)))).*rand - 1;            
                    [pOpt(:,tryK), ~, flagList(tryK)] = fmincon(fMinFcn,pInit1,[],[],[],[],-maxVal*ones(numParams,1),[maxVal*ones(numParams-1,1); lb],[],solverOptions);
                    if (pOpt(numParams,tryK) < -ub || pOpt(numParams,tryK) > -lb) 
                        success = 1;
                    end
                end
            end
            
            % Evaluate the fit
            fOpt1 = evalSmoothingFcn(pOpt(:,1),tList,'rat11');
            R1 = fOpt1 - yList;
            rmsdList(tryK) = sqrt((R1'*R1)/(length(tList)-size(pOpt,1)));

            tryK = tryK+1;

            
            % Fit to the data
            success = 0;
            while ~success
                try
                    pInit2 = pInit + exp(round(log10(abs(pInit))+2)).*rand(size(pInit));
                    pInit2(numParams) = -lb + exp(round(log10(abs(lb)))).*rand + 1;
                    [pOpt(:,tryK), ~, flagList(tryK)] = fmincon(fMinFcn,pInit2,[],[],[],[],[-maxVal*ones(numParams-1,1); ub],maxVal*ones(numParams,1),[],solverOptions);
                    if (pOpt(numParams,tryK) < -ub || pOpt(numParams,tryK) > -lb) 
                        success = 1;
                    end
                end
            end
            
            % Evaluate the fit
            fOpt2 = evalSmoothingFcn(pOpt(:,tryK),tList,'rat11');
            R2 = fOpt2 - yList;
            rmsdList(tryK) = sqrt((R2'*R2)/(length(tList)-size(pOpt,1)));

            tryK = tryK+1;
        end
        
        reportStruct.params = pOpt;
        reportStruct.rmsd = rmsdList;
        reportStruct.flagList = flagList;
        
        % Pick the best fit
        [rmsd, idx] = min(rmsdList);
        paramSolution(1:size(pOpt,1)) = pOpt(:,idx);
        
        
        
    elseif strcmp(funcType,'rat31')
        
        % Set up for the solver
        fitFcn = @(P,t) ((P(1)*t.^3 + P(2)*t.^2 + P(3)*t + P(4))./(t + P(5)));
        fMinFcn = @(P) sum((fitFcn(P,tList) - yList).^2);
        numParams = 5;
        
        success = 0;
        while ~success
            try
                [pInit,~,~,~,~] = nlinfit(tList,yList,fitFcn,rand(numParams,1));
                success = 1;
            end
        end

        tryK = 1;
        
        % We're going to run this multiple times to make sure we get
        % something reasonable
        for i = 1:numTries
            clear fOpt1 fOpt2 R1 R2

            % Fit to the data
            success = 0;
            while ~success
                try
                    pInit1 = pInit + exp(round(log10(abs(pInit))+2)).*rand(size(pInit));
                    pInit1(numParams) = -ub - exp(round(log10(abs(ub)))).*rand - 1;            
                    [pOpt(:,tryK), ~, flagList(tryK)] = fmincon(fMinFcn,pInit1,[],[],[],[],-maxVal*ones(numParams,1),[maxVal*ones(numParams-1,1); lb],[],solverOptions);
                    if (pOpt(numParams,tryK) < -ub || pOpt(numParams,tryK) > -lb) 
                        success = 1;
                    end
                end
            end
            
            % Evaluate the fit
            fOpt1 = evalSmoothingFcn(pOpt(:,1),tList,'rat31');
            R1 = fOpt1 - yList;
            rmsdList(tryK) = sqrt((R1'*R1)/(length(tList)-size(pOpt,1)));

            tryK = tryK+1;

            
            % Fit to the data
            success = 0;
            while ~success
                try
                    pInit2 = pInit + exp(round(log10(abs(pInit))+2)).*rand(size(pInit));
                    pInit2(numParams) = -lb + exp(round(log10(abs(lb)))).*rand + 1;
                    [pOpt(:,tryK), ~, flagList(tryK)] = fmincon(fMinFcn,pInit2,[],[],[],[],[-maxVal*ones(numParams-1,1); ub],maxVal*ones(numParams,1),[],solverOptions);
                    if (pOpt(numParams,tryK) < -ub || pOpt(numParams,tryK) > -lb) 
                        success = 1;
                    end
                end
            end
            
            % Evaluate the fit
            fOpt2 = evalSmoothingFcn(pOpt(:,tryK),tList,'rat31');
            R2 = fOpt2 - yList;
            rmsdList(tryK) = sqrt((R2'*R2)/(length(tList)-size(pOpt,1)));

            tryK = tryK+1;
        end
        
        reportStruct.params = pOpt;
        reportStruct.rmsd = rmsdList;
        reportStruct.flagList = flagList;
        
        % Pick the best fit
        [rmsd, idx] = min(rmsdList);
        paramSolution(1:size(pOpt,1)) = pOpt(:,idx);
        
        
        
    elseif strcmp(funcType,'rat22')  

        % Some function Definitons
        fitFcn = @(P,t) ((P(1)*t.^2 + P(2)*t + P(3))./(t.^2 + P(4)*t + P(5)));
        fMinFcn = @(P) sum((fitFcn(P,tList) - yList).^2);
        
        % Set some parameters
        pOpt = nan(5,4);
        
        % Case 1) Straddling roots
        % P(5) < -(ub * P(4) + ub^2)
        % P(5) < -(lb * P(4) + lb^2)
        A1 = [1 0 0 0 0;
              0 1 0 0 0;
              0 0 1 0 0;
              0 0 0 lb 1;
              0 0 0 ub 1];
        b1 = [maxVal; maxVal; maxVal; -lb^2; -ub^2];
        
        % Case 2) Roots both > ub
            % P(4) < -2*ub
            % P(5) > -(ub*P(4)+ub^2)
        A2 = [1 0 0 0 0;
              0 1 0 0 0;
              0 0 1 0 0;
              0 0 0 1 0;
              0 0 0 -ub -1];
        b2 = [maxVal; maxVal; maxVal; -2*ub; ub^2];
        
        % Case 3) Roots both < lb
            % P(4) > -2*lb
            % P(5) > -(ub*P(4)+ub^2)
            A3 = [1 0 0 0 0;
                  0 1 0 0 0;
                  0 0 1 0 0;
                  0 0 0 -1 0;
                  0 0 0 -lb -1];
            b3 = [maxVal; maxVal; maxVal; 2*lb; lb^2];

        tryK = 1;

        for i = 1:numTries
            clear pInit1 pInit2 pInit3 pInit4 R1 R2 R3 R4
            
            % Case 1) Straddling roots
            success = 0;
            while ~success
                try
                    pInit1 = [10*rand-5; 10*rand-5; 10*rand-5;...
                              -(ub^2-lb^2)/(ub-lb) - lb;...
                              (ub^2*lb - lb^2*ub)/(ub-lb) - ub];
                    [pOpt(:,tryK), ~, flagList(tryK)] = fmincon(fMinFcn,pInit1,A1,b1,[],[],[],[],[],solverOptions);
                    
                    pRoots = roots([1 pOpt(4:5,tryK)']);
                    if ~isreal(pRoots) || all(~(pRoots>lb & pRoots<ub)) 
                        success = 1;
                    end
                catch
                    disp('caught case 1')
                end
            end

            R1 = evalSmoothingFcn(pOpt(:,tryK),tList,'rat22') - yList;
            rmsdList(tryK) = sqrt((R1'*R1)/(length(tList)-size(pOpt,1)));
           
            tryK = tryK+1;

            % Case 2) Roots both > ub
            success = 0;
            while ~success
                try
                    pInit2 = [10*rand-5; 10*rand-5; 10*rand-5; -2*ub-1; ub^2+ub+1];
                    [pOpt(:,tryK), ~, flagList(tryK)] = fmincon(fMinFcn,pInit2,A2,b2,[],[],[],[],[],solverOptions);
                    
                    pRoots = roots([1 pOpt(4:5,tryK)']);
                    if ~isreal(pRoots) || all(~(pRoots>lb & pRoots<ub))
                        success = 1;
                    end
                catch
                    disp('caught case 2')
                end
            end

            R2 = evalSmoothingFcn(pOpt(:,tryK),tList,'rat22') - yList;
            rmsdList(tryK) = sqrt((R2'*R2)/(length(tList)-size(pOpt,1)));
            
            tryK = tryK+1;


            % Case 3) Roots both < lb
            success = 0;
            while ~success
                try
                    pInit3 = [10*rand-5; 10*rand-5; 10*rand-5; -2*lb+1; lb^2-lb+1];
                    [pOpt(:,tryK), ~, flagList(tryK)] = fmincon(fMinFcn,pInit3,A3,b3,[],[],[],[],[],solverOptions);
                    
                    pRoots = roots([1 pOpt(4:5,tryK)']);
                    if ~isreal(pRoots) || all(~(pRoots>lb & pRoots<ub))
                        success = 1;
                    end
                catch
                    disp('caught case 3')
                end
            end

            R3 = evalSmoothingFcn(pOpt(:,tryK),tList,'rat22') - yList;
            rmsdList(tryK) = sqrt((R3'*R3)/(length(tList)-size(pOpt,1)));
            
            tryK = tryK+1;


            % Case 4) Imaginary roots in the center
            success = 0;
            while ~success
                try
                    pInit4 = [10*rand-5; 10*rand-5; 10*rand-5; 0; (ub-lb)/2];
                    [pOpt(:,tryK), ~, flagList(tryK)] = fmincon(fMinFcn,pInit4,[],[],[],[],[-inf; -inf; -inf; -2*ub; 0],[inf; inf; inf; -2*lb; inf],@nonlConFuncRat22,solverOptions);
                    
                    pRoots = roots([1 pOpt(4:5,tryK)']);
                    if  ~isreal(pRoots) || all(~(pRoots>lb & pRoots<ub))
                        success = 1;
                    end
                catch
                    disp('caught case 4')
                end
            end

            R4 = evalSmoothingFcn(pOpt(:,tryK),tList,'rat22') - yList;
            rmsdList(tryK) = sqrt((R4'*R4)/(length(tList)-size(pOpt,1)));
           
            tryK = tryK+1;

        end
        
        % Now that we've looked at all solution regions in parameter
        % space...
        rmsdList(rmsdList == 0) = inf;
        
        reportStruct.params = pOpt;
        reportStruct.rmsd = rmsdList;
        reportStruct.flagList = flagList;
        
        [rmsd, idxOpt] = min(rmsdList);
        paramSolution(1:size(pOpt,1)) = pOpt(:,idxOpt);
        
        
    % "Final" impulse code    
    elseif strcmp(funcType,'impls')
        % P = [h0, h1, h2, t1, t2,  b1, b2]';
        
        % Define our intermediate functions
        sigmoid = @(t,tm,hs,he,b) hs + (he-hs) ./ (1 + exp(-4*b*(t - tm)));
        fitFcn = @(P,t) 1/P(2) * sigmoid(t,P(4),P(1),P(2),P(6)) .* sigmoid(t,P(5),P(3),P(2),P(7)) ;
        nLinConFcn = @(P) fitFcn(P,tList) - yList; 
        
        % Just in case - don't permit h<0
        minVal = max([0, min(yList) - hFactor*range(yList)]);
        
        % Place bounds on acceptable parameter values
        pLB = [minVal
               minVal
               minVal
               min(tList)
               min(tList)
               -betaFactor; 
               -betaFactor;
                ];
            
        pUB = [max(yList) + hFactor*range(yList)
               max(yList) + hFactor*range(yList)
               max(yList) + hFactor*range(yList)
               max(tList)
               max(tList)
               betaFactor; 
               betaFactor;
                ];
        
        tryK = 1;

        % Seed the solver with multiple sets of random initial parameters
        for i = 1:numTries
            success = 0;
            while ~success
                
                % If we want to seed using a parameter set from outside
                if exist('params','var') && isfield(params,'pInit')
                    if tryK == 1
                        pInitK = params.pInit;
                    else
                        
                        pInitK = params.pInit + [1/2*range(yList)*rand - 1/2*range(yList)
                                                 1/2*range(yList)*rand - 1/2*range(yList)
                                                 1/2*range(yList)*rand - 1/2*range(yList)
                                                 1/2*range(tList)*rand - 1/2*range(tList)
                                                 1/2*range(tList)*rand - 1/2*range(tList)
                                                 betaFactor * randn
                                                 betaFactor * randn
                                                 ];
                    end
                    
                % Otherwise, let's just pick something in the general
                % vicinity of the data
                else
                    rng(i)
                    y0 = yList(1);
                    yM = yList(floor(1/2*length(tList)));
                    yF = yList(end);
                    t1 = tList(floor(1/3*length(tList)));
                    t2 = tList(ceil(2/3*length(tList)));
                    beta0 = betaFactor * randn;

                    pInitK = [y0 + 1/2*range(yList)*rand - 1/2*range(yList)
                              yM + 1/2*range(yList)*rand - 1/2*range(yList)
                              yF + 1/2*range(yList)*rand - 1/2*range(yList)
                              t1 + 1/2*range(tList)*rand - 1/2*range(tList)
                              t2 + 1/2*range(tList)*rand - 1/2*range(tList)
                              beta0
                              -beta0
                              ];

                end
                
                % Perform the actual fit
                [pOpt(:,tryK), ~, R, flagList(tryK)] = lsqnonlin(nLinConFcn,pInitK,pLB,pUB,optimset('Display','off','TolFun', TolFun,'TolX', TolX, 'MaxFunEvals', MaxFunEvals, 'MaxIter', MaxIter));

                % If the output is okay
                if any(flagList(tryK) == [0 1 2 3 4])
                    success = 1;
                end

            end
            
            % Calc and record the RMSD
            rmsdList(tryK) = sqrt((R'*R)/(length(tList)-size(pOpt,1)));

            % Enforce t2>t1
            if pOpt(4,tryK)>pOpt(5,tryK)
                pOpt(:,tryK) = pOpt([3,2,1,5,4,7,6],tryK);
            end
            
            % Get ready to move on
            tryK = tryK+1;
            
        end
        
        reportStruct.params = pOpt;
        reportStruct.rmsd = rmsdList;
        reportStruct.flagList = flagList;
        
        [rmsd, idxOpt] = min(rmsdList);
        paramSolution = pOpt(:,idxOpt);
        
        
    else
        error('Invalid fitting function, ''%s''',funcType)
    end
    
end

function [c, ceq] = nonlConFuncRat22(P)
    c = P(4).^2 - 4*P(5);
    ceq = 0;
end