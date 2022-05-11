function [RF knownMet res] = solveRF_KEapproach_UDreg(relative_concMatrix,timeVec,knownKinetics)
% Kinetic equations approach, assuming some or all kinetic equations for
% reactions are known in the system.

params.v1 = 1;
params.v2M = 0.6;
params.v2K = 0.5;
params.v3M = 0.8;
params.v3K = 0.4;
params.v3I = 0.7;
params.v4M = 0.3;
params.v4K = 0.5;
params.v5M = 0.6;
params.v5K = 0.4;
params.v6M = 0.4;
params.v6K = 0.3;
params.v7M = 0.5;
params.v7K = 0.3;
params.v8M = 0.7;
params.v8K = 0.4;
params.v8A = 0.6;
params.v8alpha = 0.3;
params.v8beta = 1.3;

for t = 2:size(relative_concMatrix,1)
    relative_Vpool(t-1,:) = (relative_concMatrix(t,:)-relative_concMatrix(t-1,:))./(timeVec(t)-timeVec(t-1));
end

knownMet = findKnownMet(knownKinetics);

minF = max(relative_concMatrix(:,knownMet))/50;

lb = minF;
ub = 1000*ones(length(knownMet),1);
x0 = (1000-1).*rand(1,length(knownMet)) + 1;

if ~isempty(knownMet)
    [RF,res] = lsqnonlin(@(F) dxdt_functions(F,relative_Vpool,relative_concMatrix,params,knownKinetics,knownMet),x0,lb,ub);
else
    RF = 0;
    knownMet = 0;
end

end


function func = dxdt_functions(F,relative_Vpool,relative_concMatrix,params,knownKinetics,knownMet)

nT = size(relative_Vpool,1);
fnum = zeros(size(relative_Vpool,2),1);
fnum(knownMet) = 1:length(knownMet);

func = [];
for t = 2:nT
    
    % dxdt1
    if all(ismember([1 2 3],knownKinetics))
        func(end+1) = ...
            abs(params.v1...
            -params.v2M * relative_concMatrix(t,1)/F(fnum(1)) / (params.v2K + relative_concMatrix(t,1)/F(fnum(1)))...
            -params.v3M * relative_concMatrix(t,1)/F(fnum(1)) / ((1 + relative_concMatrix(t,2)/F(fnum(2)) / params.v3I) * (params.v3K + relative_concMatrix(t,1)/F(fnum(1))))...
            -relative_Vpool(t,1)/F(fnum(1)));
    end
    
    % dxdt2 
    if all(ismember([2 4 8],knownKinetics))
        func(end+1) = ...
            abs(params.v2M * relative_concMatrix(t,1)/F(fnum(1)) / (params.v2K + relative_concMatrix(t,1)/F(fnum(1)))...
            -params.v4M * relative_concMatrix(t,2)/F(fnum(2)) / (params.v4K + relative_concMatrix(t,2)/F(fnum(2)))...
            -params.v8M * relative_concMatrix(t,2)/F(fnum(2)) / (params.v8K * ((1 + relative_concMatrix(t,4)/F(fnum(4)) / params.v8A) / (1 + params.v8beta * relative_concMatrix(t,4)/F(fnum(4)) / (params.v8alpha * params.v8A)))...
                + relative_concMatrix(t,2)/F(fnum(2)) * ((1 + relative_concMatrix(t,4)/F(fnum(4)) / (params.v8alpha * params.v8A)) / (1 + params.v8beta * relative_concMatrix(t,4)/F(fnum(4)) / (params.v8alpha * params.v8A))))...
            -relative_Vpool(t,2)/F(fnum(2)));
    end
        
    % dxdt3
    if all(ismember([3 5 6],knownKinetics))
        func(end+1) = ...
            abs(params.v3M * relative_concMatrix(t,1)/F(fnum(1)) / ((1 + relative_concMatrix(t,2)/F(fnum(2)) / params.v3I) * (params.v3K + relative_concMatrix(t,1)/F(fnum(1))))...
            -params.v5M * relative_concMatrix(t,3)/F(fnum(3)) / (params.v5K + relative_concMatrix(t,3)/F(fnum(3)))...
            -params.v6M * relative_concMatrix(t,3)/F(fnum(3)) / (params.v6K + relative_concMatrix(t,3)/F(fnum(3)))...
            -relative_Vpool(t,3)/F(fnum(3)));
    end
     
    % dxdt4    
    if all(ismember([6 7],knownKinetics))
        func(end+1) = ...
            abs(params.v6M * relative_concMatrix(t,3)/F(fnum(3)) / (params.v6K + relative_concMatrix(t,3)/F(fnum(3)))...
            -params.v7M * relative_concMatrix(t,4)/F(fnum(4)) / (params.v7K + relative_concMatrix(t,4)/F(fnum(4)))...
            -relative_Vpool(t,4)/F(fnum(4)));
    end
    
    if isempty(func)
        func = 0;
    end
end

end

function knownMet = findKnownMet(knownKinetics)
    knownMet = [];

    % dxdt1
    if all(ismember([1 2 3],knownKinetics))
        knownMet = [knownMet 1 2];
    end
    
    % dxdt2 
    if all(ismember([2 4 8],knownKinetics)) && max(knownKinetics)
        knownMet = [knownMet 1 2 4];
    end
        
    % dxdt3
    if all(ismember([3 5 6],knownKinetics)) && max(knownKinetics)
        knownMet = [knownMet 1 2 3];
    end
     
    % dxdt4    
    if all(ismember([6 7],knownKinetics)) && max(knownKinetics)
        knownMet = [knownMet 3 4];
    end

    knownMet = unique(knownMet);
end