function [timeVec, concMatrix, fluxMatrix] = solveOdeMmUDreg(tStart,tEnd,nT,x0,params)
% This function is hard-coded to integrate the Underdetermined model with
% regulation as an ODE with Michaelis-Menten kinetic rate laws. 
%
% Written by R.A. Dromms 2015-07-29
% Edited JL10Aug20

    if exist('params','var')
        params = convertOdeParams(params);
    else
        params = setOdeParams;
    end
    
    [timeVec,concMatrix] = ode45(@(t,x)fRHS(t,x,params),linspace(tStart,tEnd,nT+1),x0');
    
    fluxMatrix = zeros(length(timeVec),10);
    for k = 1:length(timeVec)
        fluxMatrix(k,1:8) = calcFluxes(timeVec(k),concMatrix(k,:)',params)';
        fluxMatrix(k,9:12) = fRHS(timeVec(k),concMatrix(k,:)',params)';
    end
    
end

function xdot = fRHS(~,x,params)

    v = calcFluxes([],x,params);
    
    xdot = params.S*v; 

end

function v = calcFluxes(~,x,params)
    
    % Bad Things can happen with the math when we allow negative x
    % (Which isn't physically relevant, anyways)
    xMin = 1e-4;
    x(x<xMin) = xMin;

    v(1,1) = params.v0;
    v(2,1) = params.v2M * x(1) / (params.v2K + x(1)) ;
    v(3,1) = params.v3M * x(1) / ((1 + x(2) / params.v3I) * (params.v3K + x(1)));
    v(4,1) = params.v4M * x(2) / (params.v4K + x(2)) ;
    v(5,1) = params.v5M * x(3) / (params.v5K + x(3)) ;
    v(6,1) = params.v6M * x(3) / (params.v6K + x(3)) ;
    v(7,1) = params.v7M * x(4) / (params.v7K + x(4)) ;
    v(8,1) = params.v8M * x(2) / (params.v8K * ((1 + x(4) / params.v8A) / (1 + params.v8beta * x(4) / (params.v8alpha * params.v8A)))...
        + x(2) * ((1 + x(4) / (params.v8alpha * params.v8A)) / (1 + params.v8beta * x(4) / (params.v8alpha * params.v8A))));
    
end

function params = convertOdeParams(paramVec)
% ParamVec = [v2M; v2K; v3M; v3K; v3I; v4M; v4K; v5M; v5K; v6M; v6K; v7M; v7k; v8M; v8K; v8A; v8alpha; v8beta;]


    paramVec = paramVec(:);
    
    params.S = sparse([ 1 -1 -1  0  0  0  0  0;
                       0  1  0 -1  0  0  0 -1;
                       0  0  1  0 -1 -1  0  0;
                       0  0  0  0  0  1 -1  0;]);

    params.v0 = 1;
    
    params.v2M = paramVec(1);
    params.v2K = paramVec(2);
    params.v3M = paramVec(3);
    params.v3K = paramVec(4);
    params.v3I = paramVec(5);
    params.v4M = paramVec(6);
    params.v4K = paramVec(7);
    params.v5M = paramVec(8);
    params.v5K = paramVec(9);
    params.v6M = paramVec(10);
    params.v6K = paramVec(11);
    params.v7M = paramVec(12);
    params.v7K = paramVec(13);
    params.v8M = paramVec(14);
    params.v8K = paramVec(15);
    params.v8A = paramVec(16);
    params.v8alpha = paramVec(17);
    params.v8beta = paramVec(18);
    
end
