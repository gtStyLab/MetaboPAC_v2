function [dyList] = evalSmoothingDeriv(P,tList,functionType)
% Evaluates the derivative of the function listed in FUNCTIONTYPE using the
% parameters listed in P at times listed in TLIST. The dimensions of
% YLIST match those of TLIST.

    dyList = NaN(size(tList));

    if strcmp(functionType,'poly1')
        dyList = polyval(polyder(P(1:2)),tList);      
    elseif strcmp(functionType,'poly2')
        dyList = polyval(polyder(P(1:3)),tList);
    elseif strcmp(functionType,'poly3')
        dyList = polyval(polyder(P(1:4)),tList);
    elseif strcmp(functionType,'poly4')
        dyList = polyval(polyder(P(1:5)),tList);
    elseif strcmp(functionType,'poly5')
        dyList = polyval(polyder(P(1:6)),tList);
        
    elseif strcmp(functionType,'rat11') || strcmp(functionType,'rbt11') || strcmp(functionType,'rot11')
%         dyList = ((P(1)*tList + P(2))./(tList + P(3)));
        dyList = ((tList+P(3)).*P(1)-(P(1)*tList+P(2)))./(tList+P(3)).^2;
        
    elseif strcmp(functionType,'rat22') || strcmp(functionType,'rbt22') || strcmp(functionType,'rot22')
%         dyList = ((P(1)*tList.^2 + P(2)*tList + P(3))./(tList.^2 + P(4)*tList + P(5)));
        dyList = ((tList.^2+P(4)*tList+P(5)).*(2*P(1)*tList+P(2))-(P(1)*tList.^2+P(2)*tList+P(3)).*(2*tList+P(4)))./(tList.^2+P(4)*tList+P(5)).^2;
        
    elseif strcmp(functionType,'rat31') || strcmp(functionType,'rbt31') || strcmp(functionType,'rot31')
%         dyList = ((P(1)*tList.^3 + P(2)*tList.^2 + P(3)*tList + P(4))./(tList + P(5)));
        dyList = ((tList+P(5)).*(3*P(1)*tList.^2+2*P(2)*tList+P(3))-(P(1)*tList.^3+P(2)*tList.^2+P(3)*tList+P(4)))./(tList+P(5)).^2;
        
    elseif strcmp(functionType,'impls') || strcmp(functionType,'impl1') || strcmp(functionType,'impl2') || strcmp(functionType,'impl0') || strcmp(functionType,'implF')
        % P = [h0, h1, h2, t1, t2, b1, b2]';
        sigmoid = @(t,tm,hs,he,b) hs + (he-hs) ./ (1 + exp(-4*b*(t - tm)));
        ddtSigmoid = @(t,tm,hs,he,b) (he-hs) .* ((4*b*exp(-4*b*(t - tm))) ./ (1 + exp(-4*b*(t - tm))).^2);
        
        % So apparently the ddt denominator blows up sometimes and ends up with NaN values
        term1a = sigmoid(tList,P(4),P(1),P(2),P(6));
        term1b = ddtSigmoid(tList,P(5),P(3),P(2),P(7));
        term2a = ddtSigmoid(tList,P(4),P(1),P(2),P(6));
        term2b = sigmoid(tList,P(5),P(3),P(2),P(7));
        
        if any(isnan(term1b))
            term1b(isnan(term1b)) = zeros(1,nnz(isnan(term1b)));
        end

        if any(isnan(term2a))
            term2a(isnan(term2a)) = zeros(1,nnz(isnan(term2a)));
        end
        
        if any(isnan(term1a)) || any(isnan(term2b))
            error
        end
        
        dyList = 1/P(2) * (term1a .* term1b + term2a .* term2b) ;

        
    end
    
end