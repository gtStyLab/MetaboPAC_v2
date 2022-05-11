function penalty = pen_oneContMetCurveFit_hynne(concMatrix,fluxMatrix,timeVec,fluxTimeVec)
% Calculate curve fit penalty between metabolite and flux.

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
warning('off','all')
for interaction_iter = 1:length(oneContMet_interactions)
    met = oneContMet_interactions{interaction_iter,1};
    flux = oneContMet_interactions{interaction_iter,2};
    
    x = concMatrix(1:end-1,met);
    v = fluxMatrix(:,flux);
    
    [f gof] = fit(x,v,'poly2');
    penalty(interaction_iter) = 1 - gof.adjrsquare;
end
penalty = nansum(penalty);