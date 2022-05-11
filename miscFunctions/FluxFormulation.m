function v = FluxFormulation(x,powerLawKineticsParams,stoichMatrix)
    [numFlux,numMetp1] = size(powerLawKineticsParams);
    v = zeros(numFlux,1);
    for i = 1:1:numFlux
        v(i,1) = powerLawKineticsParams(i,1);
        for j = 2:1:numMetp1
            v(i,1) = v(i,1) * (x(j-1).^ powerLawKineticsParams(i,j));
        end
    end
    
%     for flux = 1:1:size(stoichMatrix,2)
%         mass_action_mets = find(stoichMatrix(:,flux)>0);
%         for met = 1:1:length(mass_action_mets)
%             if x(mass_action_mets(met)) < 1e-04
%                 v(flux,1) = 0;
%             end
%             
%         end
%     end
   
    
end