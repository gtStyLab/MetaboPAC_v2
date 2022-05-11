function BST_penalty = pen_BSTfit_hynne(absolute_concMatrix,Vcalc)
% Calculate penalty for fitting BST equation to metabolite and flux data.
    nT = size(Vcalc,1); %#timepoints = #timepoints in flux profiles 
    numFlux = size(Vcalc,2) - size(absolute_concMatrix,2); 
    numMet = size(absolute_concMatrix,2);

kineticsMap = {[2,8,17]      ,2;
                [2,21]       ,3;
                [8,12]       ,4;
                 [12,20,21]  ,5;
                 [6,7,19]    ,6;
                 [6,7]       ,7;
                 [6,9,18,22] ,8;
                 [5,11,18,21],9;
                 [5,11]      ,10;
                 [14]        ,11;
                 [1,22]      ,12;
                 [4,15]      ,13;
                 [15]        ,14;
                 [7,9,22]    ,15;
                 [3,16]      ,16;
                 [16]        ,17;
                 [1,10]      ,18;
                 [10]        ,19;
                 [10,13]     ,20;
                 [13]        ,21;
                 [8,21]      ,22;
                 [21]        ,23;
                 [5,20,21]   ,24;};
     %% Estimate BST parameters 
     %Preassign BST parameters 
     BST_params = zeros(numFlux,numMet + 1);
     for i = 1:1:size(kineticsMap,1) %for each flux, calculate the parameters 
         
         controllerMets = kineticsMap{i,1};
         contFlux = kineticsMap{i,2};
         numCols = length(controllerMets);
         %if a certain flux contains a negative value, bypass it 
         if all(Vcalc(:,contFlux) > 0)
             %preassign A 
             A = nan(nT,numCols+1);
             for j = 1:1:numCols+1
                 if isequal(j,1)
                     A(:,j) = ones(nT,1);
                 else
                     contMet = controllerMets(j-1);
                     A(:,j) = log(absolute_concMatrix(1:nT,contMet));
                 end

             end
             %Finishes assigning A 
             b = log(Vcalc(:,contFlux));
             %Calculate parameters and put in BST_params
             params = A\b;
             %log(v) = log(a) + b1log(x_i) +...


             params(1) = exp(params(1));
             BST_params(contFlux,1) = params(1);
             for k = 1:1:numCols
                 contMet = controllerMets(k);
                 BST_params(contFlux,contMet+1) = params(k+1);
             end
         else
            continue
         end       
     end
     %Finish assigning BST_params 
     %% Calculating fitted fluxes 
     BST_flux = nan(nT,numFlux);
     for timeIdx = 1:nT
        v = FluxFormulation(absolute_concMatrix(timeIdx,:),BST_params);
        BST_flux(timeIdx,:) = v';
     end
     
     %% Calculate BST penalty 
     BST_error = zeros(1,numFlux);
     for fluxIdx = 2:1:numFlux
        if all(Vcalc(:,fluxIdx) > 0)
            % calculate RMSD for each flux
            BST_error(fluxIdx) = real(sqrt(sum((Vcalc(:,fluxIdx)-BST_flux(:,fluxIdx)).^2)/nT));
        end
     end
     
     BST_penalty = sum(BST_error);
     
end