function [RF knownMet res] = solveRF_KEapproach_generic(relative_concMatrix,timeVec,knownKinetics)
% Kinetic equations approach, assuming some or all kinetic equations for
% reactions are known in the system.

% List parameters in kinetic equations
params.v1M = ;
params.v1K = ;

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

% Write out the exact kinetic equations for each dx/dt term that is fully
% known.
% Relative abundances are divided by the predicted response factor (RF)
% to convert these abundances into inferred absolute concentrations.
function func = dxdt_functions(RF,relative_Vpool,relative_concMatrix,params,knownKinetics,knownMet)

nT = size(relative_Vpool,1);
RFnum = zeros(size(relative_Vpool,2),1);
RFnum(knownMet) = 1:length(knownMet);

func = [];
for t = 2:nT
    
    % dxdt1
    if all(ismember([],knownKinetics))
        func(end+1) = ...
            % Example of sum of influx and efflux kinetic equations for
            % metabolite x1
            abs(params.v1...
            -params.v2M * relative_concMatrix(t,1)/RF(RFnum(1)) / (params.v2K + relative_concMatrix(t,1)/RF(RFnum(1)))...
            -params.v3M * relative_concMatrix(t,1)/RF(RFnum(1)) / (params.v3K + relative_concMatrix(t,1)/RF(RFnum(1)))...
            -relative_Vpool(t,1)/RF(RFnum(1)));
    end
    
    if isempty(func)
        func = 0;
    end
end

end

% List the known kinetic equations that are required to calculate each
% dx/dt term based on the stoichiometry of the system.
% This process could likely be "soft-coded" in the future.
function knownMet = findKnownMet(knownKinetics)
    knownMet = [];

    % dxdt1
    if all(ismember([],knownKinetics))
        knownMet = [knownMet 1];
    end

    knownMet = unique(knownMet);
end