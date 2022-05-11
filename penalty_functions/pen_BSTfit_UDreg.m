function BST_penalty = pen_BSTfit_UDreg(nT,absolute_concMatrix,Vcalc)
% Calculate penalty for fitting BST equation to metabolite and flux data.

% Set up linear equations
Av2 = [ones(nT,1), log(absolute_concMatrix(1:nT,1))];
Av3 = [ones(nT,1), log(absolute_concMatrix(1:nT,1)), log(absolute_concMatrix(1:nT,2))];
Av4 = [ones(nT,1), log(absolute_concMatrix(1:nT,2))];
Av5 = [ones(nT,1), log(absolute_concMatrix(1:nT,3))];
Av6 = [ones(nT,1), log(absolute_concMatrix(1:nT,3))];
Av7 = [ones(nT,1), log(absolute_concMatrix(1:nT,4))];
Av8 = [ones(nT,1), log(absolute_concMatrix(1:nT,2)), log(absolute_concMatrix(1:nT,4))];
A = blkdiag(Av2, Av3, Av4, Av5, Av6, Av7, Av8);
b = reshape(log(Vcalc(:,2:8)),[numel(Vcalc(:,2:8)),1]);

% Solve for linearized best fit
x = A\b;

% Switch a terms out of log-space
x([1 3 6 8 10 12 14]) = exp(x([1 3 6 8 10 12 14]));

% Calculate predicted BST flux
BST_flux(:,1) = x(1).*absolute_concMatrix(1:nT,1).^x(2);
BST_flux(:,2) = x(3).*absolute_concMatrix(1:nT,1).^x(4).*absolute_concMatrix(1:nT,2).^x(5);
BST_flux(:,3) = x(6).*absolute_concMatrix(1:nT,2).^x(7);
BST_flux(:,4) = x(8).*absolute_concMatrix(1:nT,3).^x(9);
BST_flux(:,5) = x(10).*absolute_concMatrix(1:nT,3).^x(11);
BST_flux(:,6) = x(12).*absolute_concMatrix(1:nT,4).^x(13);
BST_flux(:,7) = x(14).*absolute_concMatrix(1:nT,2).^x(15).*absolute_concMatrix(1:nT,4).^x(16);

% Calculate error between predicted flux using BST and inferred flux
BST_error(1) = real(sqrt(sum(sum((Vcalc(1:nT,2) - BST_flux(:,1)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,2) - BST_flux(:,1)).^2)))));
BST_error(2) = real(sqrt(sum(sum((Vcalc(1:nT,3) - BST_flux(:,2)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,3) - BST_flux(:,2)).^2)))));
BST_error(3) = real(sqrt(sum(sum((Vcalc(1:nT,4) - BST_flux(:,3)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,4) - BST_flux(:,3)).^2)))));
BST_error(4) = real(sqrt(sum(sum((Vcalc(1:nT,5) - BST_flux(:,4)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,5) - BST_flux(:,4)).^2)))));
BST_error(5) = real(sqrt(sum(sum((Vcalc(1:nT,6) - BST_flux(:,5)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,6) - BST_flux(:,5)).^2)))));
BST_error(6) = real(sqrt(sum(sum((Vcalc(1:nT,7) - BST_flux(:,6)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,7) - BST_flux(:,6)).^2)))));
BST_error(7) = real(sqrt(sum(sum((Vcalc(1:nT,8) - BST_flux(:,7)).^2)))) + abs(imag(sqrt(sum(sum((Vcalc(1:nT,8) - BST_flux(:,7)).^2)))));

BST_penalty = nansum(BST_error);