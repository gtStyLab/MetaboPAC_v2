function penalty = pen_oneContMetCurveFit_determined(concMatrix,fluxMatrix,timeVec,fluxTimeVec)
% Calculate curve fit penalty between metabolite and flux.

% First column is controller metabolite, second column is target flux
% that the controller metabolite interacts with.
oneContMet_interactions = {[1],2;
                           [1],3;
                           [2],4;
                           [3],5;};
                       
penalty = nan(length(oneContMet_interactions),1);
warning('off','all')
for interaction_iter = 1:length(oneContMet_interactions)
    met = oneContMet_interactions{interaction_iter,1};
    flux = oneContMet_interactions{interaction_iter,2};
    
    x = concMatrix(1:end-1,met);
    v = fluxMatrix(:,flux);
    
    [f gof] = fit(x,v,'poly2');
    penalty(interaction_iter) = 1 - gof.adjrsquare;
    if f.p1 > 0
        penalty(interaction_iter) = penalty(interaction_iter) + 10*f.p1;
    end
end
penalty = nansum(penalty);