function penalty = pen_oneContMetCorr_UDreg(concMatrix,fluxMatrix,timeVec,fluxTimeVec)
% Calculate correlation penalty between metabolite and flux.

% First column is controller metabolite, second column is target flux
% that the controller metabolite interacts with.
oneContMet_interactions = {[1],2;
                           [2],4;
                           [3],5;
                           [3],6;
                           [4],7;};

penalty = nan(length(oneContMet_interactions),1);
for interaction_iter = 1:length(oneContMet_interactions)
    met = oneContMet_interactions{interaction_iter,1};
    flux = oneContMet_interactions{interaction_iter,2};
    
    x = concMatrix(1:end-1,met);
    v = fluxMatrix(:,flux);
    
    if flux == 2
        penalty(1) = abs(corr(x,v,'Type','Spearman') - 1);
    elseif flux == 4
        penalty(2) = abs(corr(x,v,'Type','Spearman') - 1);
    elseif flux == 5
        penalty(3) = abs(corr(x,v,'Type','Spearman') - 1);
    elseif flux == 6
        penalty(4) = abs(corr(x,v,'Type','Spearman') - 1);
    elseif flux == 7
        penalty(5) = abs(corr(x,v,'Type','Spearman') - 1);
    end

end
penalty = nansum(penalty);