function [timeVec, concMatrix, fluxMatrix] = solveOdeMmDetermined(tStart,tEnd,nT,x0,params)
% This function is hard-coded to integrate the Determined model as an
% ODE with Michaelis-Menten kinetic rate laws. 
%
% Written by R.A. Dromms 2015-07-29
% Edited JL03Aug20

    if exist('params','var')
        params = convertOdeParams(params);
    else
        params = setOdeParams;
    end
    
    [timeVec,concMatrix] = ode45(@(t,x)fRHS(t,x,params),linspace(tStart,tEnd,nT+1),x0');
    
    fluxMatrix = zeros(length(timeVec),10);
    for k = 1:length(timeVec)
        fluxMatrix(k,1:5) = calcFluxes(timeVec(k),concMatrix(k,:)',params)';
        fluxMatrix(k,6:9) = fRHS(timeVec(k),concMatrix(k,:)',params)';
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
    v(3,1) = params.v3M * x(1) / (params.v3K + x(1)) ;
    v(4,1) = params.v4M * x(2) / (params.v4K + x(2)) ;
    v(5,1) = params.v5M * x(3) / (params.v5K + x(3)) ;
    
    
end

function params = convertOdeParams(paramVec)
% ParamVec = [bm3; bm4; v2M; v2K; v3M; v3K; v3A; v3alpha; v3beta; v4M; v4K; v4I; v5M; v5K;]


    paramVec = paramVec(:);
    
    params.S = sparse([ 1 -1 -1  0  0  ;
                       0  1  0 -1  0  ;
                       0  0  1  0 -1  ;
                       0  0  0  1  0  ;]);

    params.v0 = 1;
    
    params.v2M = paramVec(1);
    params.v2K = paramVec(2);
    params.v3M = paramVec(3);
    params.v3K = paramVec(4);
    params.v4M = paramVec(5);
    params.v4K = paramVec(6);
    params.v5M = paramVec(7);
    params.v5K = paramVec(8);
    
end
