function [fluxMatrixOut] = calcFluxesViaPinv(fluxMatrixIn,STM,fixedFluxes)
% This function calculates the system fluxes over time as a function of the
% pooling fluxes at each time point and the system stoichiometry using the
% pseudoinverse function.

    fluxMatrixOut = fluxMatrixIn;
    numSysFluxes = size(STM,2);
    poolingFluxMatrix = fluxMatrixIn(:,numSysFluxes+1:end);
    
    % Create two separate STMs, one each for fixed and non-fixed fluxes
    stmFixed = STM(:,fixedFluxes);
    stmVaried = STM(:,~fixedFluxes);
    fixedFluxMatrix = fluxMatrixIn(:,fixedFluxes);
    
    % For each t_k, use pseudoinverse on the system
    for k = 1:size(fluxMatrixIn)
        fluxMatrixOut(k,find(~fixedFluxes)) = (pinv(full(stmVaried))*(poolingFluxMatrix(k,:)' - stmFixed*fixedFluxMatrix(k,:)'))';
    end

end
