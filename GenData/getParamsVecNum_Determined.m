function vecOut = getParamsVecNum_Determined(n)
% Just a hard-coded list of parameter vectors, by index

    % Vary ICs
        %  Params  = [ a2;  b21;     a3;  b31;    a4;  b42;    a5;  b53;  x0_1; x0_2; x0_3; x0_4;];
    
    paramsList{1}  = [0.6;  0.5;    0.8;  0.4;   0.3;  0.5;   0.6; 0.4;   0.4;  0.3;  0.7;  0.6; ];
    
    vecOut = paramsList{n};

end


