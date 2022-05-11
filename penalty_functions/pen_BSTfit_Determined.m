function BST_penalty = pen_BSTfit_Determined(nT,absolute_concMatrix,Vcalc,MA_reactions)
% Calculate penalty for fitting BST equation to metabolite and flux data.

for i = 1:size(MA_reactions,1)
    current_conc = absolute_concMatrix(1:nT,MA_reactions(i,1));
    current_flux = Vcalc(1:nT,MA_reactions(i,2));

    A = [ones(nT,1) log(current_conc)];
    b = log(current_flux);
    x = A\b;
    alpha = exp(x(1));
    beta = x(2);

    BST_flux = alpha*current_conc.^beta;

    BST_error(i) = real(sqrt(sum(sum((current_flux - BST_flux).^2)))) + abs(imag(sqrt(sum(sum((current_flux - BST_flux).^2)))));
end
BST_penalty = nansum(BST_error);