function penalty = pen_oneContMetCurveFit_chass(concMatrix,fluxMatrix,timeVec,fluxTimeVec)
% Calculate curve fit penalty between metabolite and flux.

% First column is controller metabolite, second column is target flux
% that the controller metabolite interacts with.
oneContMet_interactions = {[2],5;
                           [6],15;
                           [8],17;
                           [10],22;
                           [11],23;
                           [11],25;
                           [12],27;
                           [16],30;
                           [2],32;
                           [3],33;
                           [4],34;
                           [5],35;
                           [6],36;
                           [7],37;
                           [8],38;
                           [9],39;
                           [10],40;
                           [11],41;
                           [12],42;
                           [13],43;
                           [14],44;
                           [15],45;
                           [16],46;
                           [17],47;
                           [18],48};

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