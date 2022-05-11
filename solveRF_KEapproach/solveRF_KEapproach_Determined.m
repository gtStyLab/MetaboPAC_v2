function [RF knownMet res] = solveRF_KEapproach_Determined(relative_concMatrix,timeVec,knownKinetics)
% Kinetic equations approach, assuming some or all kinetic equations for
% reactions are known in the system.

params.v1 = 1;
params.v2M = 0.6;
params.v2K = 0.5;
params.v3M = 0.8;
params.v3K = 0.4;
params.v4M = 0.3;
params.v4K = 0.5;
params.v5M = 0.6;
params.v5K = 0.4;

for t = 2:size(relative_concMatrix,1)
    relative_Vpool(t-1,:) = (relative_concMatrix(t,:)-relative_concMatrix(t-1,:))./(timeVec(t)-timeVec(t-1));
end

knownMet = findKnownMet(knownKinetics);

minRF = max(relative_concMatrix(:,knownMet))/50;

lb = minRF;
ub = 1000*ones(length(knownMet),1);
x0 = (1000-1).*rand(1,length(knownMet)) + 1;

if ~isempty(knownMet)
    [RF,res] = lsqnonlin(@(RF) dxdt_functions(RF,relative_Vpool,relative_concMatrix,params,knownKinetics,knownMet),x0,lb,ub);
else
    RF = 0;
    knownMet = 0;
end

end


function func = dxdt_functions(RF,relative_Vpool,relative_concMatrix,params,knownKinetics,knownMet)

nT = size(relative_Vpool,1);
RFnum = zeros(size(relative_Vpool,2),1);
RFnum(knownMet) = 1:length(knownMet);

func = [];
for t = 2:nT
    
    % dxdt1
    if all(ismember([1 2 3],knownKinetics))
        func(end+1) = ...
            abs(params.v1...
            -params.v2M * relative_concMatrix(t,1)/RF(RFnum(1)) / (params.v2K + relative_concMatrix(t,1)/RF(RFnum(1)))...
            -params.v3M * relative_concMatrix(t,1)/RF(RFnum(1)) / (params.v3K + relative_concMatrix(t,1)/RF(RFnum(1)))...
            -relative_Vpool(t,1)/RF(RFnum(1)));
    end
    
    % dxdt2 
    if all(ismember([2 4],knownKinetics))
        func(end+1) = ...
            abs(params.v2M * relative_concMatrix(t,1)/RF(RFnum(1)) / (params.v2K + relative_concMatrix(t,1)/RF(RFnum(1)))...
            -params.v4M * relative_concMatrix(t,2)/RF(RFnum(2)) / (params.v4K + relative_concMatrix(t,2)/RF(RFnum(2)))...
            -relative_Vpool(t,2)/RF(RFnum(2)));
    end
        
    % dxdt3
    if all(ismember([3 5],knownKinetics))
        func(end+1) = ...
            abs(params.v3M * relative_concMatrix(t,1)/RF(RFnum(1)) / (params.v3K + relative_concMatrix(t,1)/RF(RFnum(1)))...
            -params.v5M * relative_concMatrix(t,3)/RF(RFnum(3)) / (params.v5K + relative_concMatrix(t,3)/RF(RFnum(3)))...
            -relative_Vpool(t,3)/RF(RFnum(3)));
    end
     
    % dxdt4    
    if all(ismember([4],knownKinetics))
        func(end+1) = ...
            abs(params.v4M * relative_concMatrix(t,2)/RF(RFnum(2)) / (params.v4K + relative_concMatrix(t,2)/RF(RFnum(2)))...
            -relative_Vpool(t,4)/RF(RFnum(4)));
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
        knownMet = [knownMet 1];
    end
    
    % dxdt2 
    if all(ismember([2 4],knownKinetics))
        knownMet = [knownMet 1 2];
    end
        
    % dxdt3
    if all(ismember([3 5],knownKinetics))
        knownMet = [knownMet 1 3];
    end
     
    % dxdt4    
    if all(ismember([4],knownKinetics))
        knownMet = [knownMet 2 4];
    end

    knownMet = unique(knownMet);
end