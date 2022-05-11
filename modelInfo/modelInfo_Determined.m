% Model/system info
modelInfo.S = sparse([ 1 -1 -1  0  0  ;
                       0  1  0 -1  0  ;
                       0  0  1  0 -1  ;
                       0  0  0  1  0  ;]);
                   
modelInfo.xBounds = [0,inf;
                     0,inf;
                     0,inf;
                     0,inf;
                     ]; 
                 
modelInfo.vBounds = [1,1;
                     0,inf;
                     0,inf;
                     0,inf;
                     0,inf;
                     ];

