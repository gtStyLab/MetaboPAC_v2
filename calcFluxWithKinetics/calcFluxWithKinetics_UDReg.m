function flux = calcFluxWithKinetics_UDReg(concMatrix,timeVec,knownKinetics)
% Calculate flux reaction rates in underdetermined system with regulation
% using known kinetic equations.

nT = size(concMatrix,1)-1;
deltaT = timeVec(2)-timeVec(1);

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
    
flux = nan(nT,12);
for t = 1:nT
    flux(t,1) = params.v1;
    
    if ismember(2,knownKinetics)
        flux(t,2) = params.v2M * concMatrix(t,1) / (params.v2K + concMatrix(t,1));
    end
    
    if ismember(3,knownKinetics)
        flux(t,3) = params.v3M * concMatrix(t,1) / ((1 + concMatrix(t,2) / params.v3I) * (params.v3K + concMatrix(t,1)));
    end
    
    if ismember(4,knownKinetics)
        flux(t,4) = params.v4M * concMatrix(t,2) / (params.v4K + concMatrix(t,2));
    end
    
    if ismember(5,knownKinetics)
        flux(t,5) = params.v5M * concMatrix(t,3) / (params.v5K + concMatrix(t,3));
    end
    
    if ismember(6,knownKinetics)
        flux(t,6) = params.v6M * concMatrix(t,3) / (params.v6K + concMatrix(t,3));
    end
    
    if ismember(7,knownKinetics)
        flux(t,7) = params.v7M * concMatrix(t,4) / (params.v7K + concMatrix(t,4));
    end
    
    if ismember(8,knownKinetics)
        flux(t,8) = params.v8M * concMatrix(t,2) / (params.v8K * ((1 + concMatrix(t,4) / params.v8A) / (1 + params.v8beta * concMatrix(t,4) / (params.v8alpha * params.v8A)))...
        + concMatrix(t,2) * ((1 + concMatrix(t,4) / (params.v8alpha * params.v8A)) / (1 + params.v8beta * concMatrix(t,4) / (params.v8alpha * params.v8A))));
    end
    
    flux(t,9) = (concMatrix(t+1,1) - concMatrix(t,1))/deltaT;
    flux(t,10) = (concMatrix(t+1,2) - concMatrix(t,2))/deltaT;
    flux(t,11) = (concMatrix(t+1,3) - concMatrix(t,3))/deltaT;
    flux(t,12) = (concMatrix(t+1,4) - concMatrix(t,4))/deltaT;
end
