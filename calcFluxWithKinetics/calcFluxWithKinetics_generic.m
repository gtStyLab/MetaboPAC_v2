function flux = calcFluxWithKinetics_generic(concMatrix,timeVec,knownKinetics)
% Calculate flux reaction rates in system of interest using known kinetic
% equations.

nT = size(concMatrix,1)-1;
deltaT = timeVec(2)-timeVec(1);

% List parameters in kinetic equations
params.v1M = ;
params.v1K = ;
    
% Set the size of the flux matrix, typically nT by # fluxes + # metabolites
flux = nan(nT,____);

% Calculate fluxes using known kinetic equations and parameters
for t = 1:nT
    
    % Flux reaction rate calculations
    if ismember(1,knownKinetics)
        flux(t,1) = params.v1M * concMatrix(t,1) / (params.v1K + concMatrix(t,1));
    end
    
    % Pooling fluxes (i.e. dx/dt) calculations
    flux(t,2) = (concMatrix(t+1,1) - concMatrix(t,1))/deltaT;
end
