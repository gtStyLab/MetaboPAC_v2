function vecOut = getParamsVecNum_UDreg(n)
% Just a hard-coded list of parameter vectors, by index

    % Vary ICs
        %  Params  = [v2M;  v2K;   v3M;  v3K;   v3I;   v4M;  v4K;   v5M;  v5K;   v6M;  v6K;  v7M;  v7K;  v8M;  v8K;  v8A;  v8alpha;   v8beta;  x0_1; x0_2; x0_3; x0_4;];
    
    paramsList{1}  = [0.6;  0.5;   0.8;  0.4;   0.7;   0.3;  0.5;   0.6; 0.4;   0.4;   0.3;  0.5;  0.3;  0.7;  0.4;  0.6;      0.3;      1.3;   0.4;  0.3;  0.7;  0.6;];
    
    vecOut = paramsList{n};

end


