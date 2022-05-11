function [yList] = evalSmoothingFcn(P,tList,functionType)
% Evaluates the function listed in FUNCTIONTYPE using parameters listed in
% P at times listed in TLIST. The dimensions of YLIST match those of TLIST.

    yList = NaN(size(tList));

    if strcmp(functionType,'poly1')
        yList = polyval(P(1:2),tList);      
    elseif strcmp(functionType,'poly2')
        yList = polyval(P(1:3),tList);
    elseif strcmp(functionType,'poly3')
        yList = polyval(P(1:4),tList);
    elseif strcmp(functionType,'poly4')
        yList = polyval(P(1:5),tList);
    elseif strcmp(functionType,'poly5')
        yList = polyval(P(1:6),tList);
    elseif strcmp(functionType,'rat11') || strcmp(functionType,'rbt11') || strcmp(functionType,'rot11')
        yList = ((P(1)*tList + P(2))./(tList + P(3)));
    elseif strcmp(functionType,'rat22') || strcmp(functionType,'rbt22') || strcmp(functionType,'rot22')
        yList = ((P(1)*tList.^2 + P(2)*tList + P(3))./(tList.^2 + P(4)*tList + P(5)));
    elseif strcmp(functionType,'rat31') || strcmp(functionType,'rbt31') || strcmp(functionType,'rot31')
        yList = ((P(1)*tList.^3 + P(2)*tList.^2 + P(3)*tList + P(4))./(tList + P(5)));
    elseif strcmp(functionType,'impls') || strcmp(functionType,'impl1') || strcmp(functionType,'impl2') || strcmp(functionType,'impl0') || strcmp(functionType,'implF')
        % P = [h0, h1, h2, t1, t2, b1, b2]';
        sigmoid = @(t,tm,hs,he,b) hs + (he-hs) ./ (1 + exp(-4*b*(t - tm)));
        yList = 1/P(2) * sigmoid(tList,P(4),P(1),P(2),P(6)) .* sigmoid(tList,P(5),P(3),P(2),P(7)) ;
    end
    
end