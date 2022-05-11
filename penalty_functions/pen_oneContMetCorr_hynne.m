function penalty = pen_oneContMetCorr_hynne(concMatrix,fluxMatrix,timeVec,fluxTimeVec)
% Calculate correlation penalty between metabolite and flux.

% First column is controller metabolite, second column is target flux
% that the controller metabolite interacts with.
oneContMet_interactions = {[17],1;
                           [14],11;
                           [15],14;
                           [10],19;
                           [16],17;
                           [13],21;
                           [21],23;};

penalty = nan(length(oneContMet_interactions),1);
for interaction_iter = 1:length(oneContMet_interactions)
    met = oneContMet_interactions{interaction_iter,1};
    flux = oneContMet_interactions{interaction_iter,2};
    
    x = concMatrix(1:end-1,met);
    v = fluxMatrix(:,flux);
    
    if flux == 1
        penalty(1) = abs(corr(x,v,'Type','Spearman') + 1); % Actually controlled by two metabolites, but one is constant
    elseif flux == 11
        penalty(2) = abs(corr(x,v,'Type','Spearman') - 1);
    elseif flux == 14
        penalty(3) = abs(corr(x,v,'Type','Spearman') - 1);
    elseif flux == 19
        penalty(4) = abs(corr(x,v,'Type','Spearman') - 1);
    elseif flux == 17
        penalty(5) = abs(corr(x,v,'Type','Spearman') - 1);
    elseif flux == 21
        penalty(6) = abs(corr(x,v,'Type','Spearman') + 1); % Actually controlled by two metabolites, but one is constant
    elseif flux == 23
        penalty(7) = abs(corr(x,v,'Type','Spearman') - 1);
    end

end
penalty = nansum(penalty);