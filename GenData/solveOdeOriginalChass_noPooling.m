function [timeVec, concMatrix, fluxMatrix] = solveOdeOriginalChass_noPooling(tStart,tEnd,nT,x0,params,stm)
% This function is hard-coded to integrate the branched pathway model as an
% ODE with Michaelis-Menten kinetic rate laws described in the
% Supplementary Methods. 
%
% Written by R.A. Dromms 2015-07-29
% Edited for Chassagnole 2002 model JL08Feb19
    %{
    if exist('params','var')
        params = convertOdeParams(params);
    else
        params = setOdeParams;
    end
    %}
    [timeVec,concMatrix] = ode15s(@(t,x)fRHS(t,x,params,stm),linspace(tStart,tEnd,nT+1),x0');
    
    % From solveOdeBst...
    %[timeVec,concMatrix] = ode15s(@(t,x)fRHS(t,x,params,stm),linspace(tStart,tEnd,nT+1),x0',odeset('NonNegative',ones(size(x0'))));
    
    fluxMatrix = zeros(length(timeVec),48);
    for k = 1:length(timeVec)
        fluxMatrix(k,1:48) = calcFluxes(timeVec(k),concMatrix(k,:)',params)';
%         fluxMatrix(k,49:66) = fRHS(timeVec(k),concMatrix(k,:)',params,stm)';
    end
    
end

function xdot = fRHS(~,x,params,stm)
    
    x(x<0) = 0;
    v = calcFluxes([],x,params);
    
    xdot = stm*v; 

end

function v = calcFluxes(~,x,params)
    
    % Bad Things can happen with the math when we allow negative x
    % (Which isn't physically relevant, anyways)
    xMin = 1e-4;
    x(x<xMin) = xMin;
    
    % List of reactions
    %Extracellular glucose kinetics
    v(1,1) = params(8)*(params(9)-x(1));    

    %Phosphotransferase system
    v(2,1) = params(10)*x(1)*(x(10)/x(11))/((params(11)+params(12)*(x(10)/x(11))+params(13)*x(1)+x(1)*(x(10)/x(11)))*(1+x(2)^params(14)/params(15)));

    %Glucose-6-phosphate isomerase
    v(3,1) = params(16)*(x(2)-x(3)/params(17))/(params(18)*(1+x(3)/(params(19)*(1+x(12)/params(20)))+x(12)/params(21))+x(2));

    %Phosphoglucomutase
    v(4,1) = params(22)*(x(2)-x(18)/params(23))/(params(24)*(1+x(18)/params(25))+x(2));

    %Glucose-6-phosphate dehydrogenase
    v(5,1) = params(26)*x(2)*params(4)/((x(2)+params(27))*(1+params(5)/params(28))*(params(29)*(1+params(5)/params(30))+params(4)));

    %Phosphofructokinase
    v(6,1) = params(31)*params(1)*x(3)/((params(1)+params(32)*(1+params(2)/params(33)))*(x(3)+params(34)*(1+x(10)/params(35)+params(2)/params(36)+params(3)/params(37))/(1+params(2)/params(38)+params(3)/params(39)))*(1+params(40)/(1+x(3)*(1+params(2)/params(38)+params(3)/params(39))/(params(34)*(1+x(10)/params(35)+params(2)/params(36)+params(3)/params(37))))^params(41)));

    %Transaldolase
    v(7,1) = params(42)*(x(5)*x(15)-x(17)*x(3)/params(43));

    %Transketolase a
    v(8,1) = params(44)*(x(16)*x(14)-x(15)*x(5)/params(45));

    %Transketolase b
    v(9,1) = params(46)*(x(14)*x(17)-x(3)*x(5)/params(47));

    %Mureine synthesis
    v(10,1) = params(48);

    %Aldolase
    v(11,1) = params(49)*(x(4)-x(5)*x(6)/params(50))/(params(51)+x(4)+params(52)*x(6)/(params(50)*params(53))+params(54)*x(5)/(params(50)*params(53))+x(4)*x(5)/params(55)+x(5)*x(6)/(params(53)*params(50)));

    %Glyceraldehyde-3-phosphate dehydrogenase
    v(12,1) = params(56)*(x(5)*params(6)-x(7)*params(7)/params(57))/((params(58)*(1+x(7)/params(59))+x(5))*(params(60)*(1+params(7)/params(61))+params(6)));
    
    %Triosephosphate isomerase
    v(13,1) = params(62)*(x(6)-x(5)/params(63))/(params(64)*(1+x(5)/params(65))+x(6));

    %Tryptophan synthesis
    v(14,1) = params(66);

    %Glycerol-3-phosphate dehydrogenase
    v(15,1) = params(67)*x(6)/(params(68)+x(6));

    %Phosphoglycerate kinase
    v(16,1) = params(69)*(params(2)*x(7)-params(1)*x(8)/params(70))/((params(71)*(1+params(1)/params(72))+params(2))*(params(73)*(1+x(8)/params(74))+x(7)));

    %Serine synthesis
    v(17,1) = params(75)*x(8)/(params(76)+x(8));

    %Phosphoglycerate Mutase
    v(18,1) = params(77)*(x(8)-x(9)/params(78))/(params(79)*(1+x(9)/params(80))+x(8));

    %Enolase
    v(19,1) = params(81)*(x(9)-x(10)/params(82))/(params(83)*(1+x(10)/params(84))+x(9));

    %Pyruvate kinase
    v(20,1) = params(85)*x(10)*(x(10)/params(86)+1)^(params(87)-1)*params(2)/(params(86)*(params(88)*((1+params(1)/params(89))/(x(4)/params(90)+params(3)/params(91)+1))^params(87)+(x(10)/params(86)+1)^params(87))*(params(2)+params(92)));

    %PEP carboxylase
    v(21,1) = params(93)*x(10)*(1+(x(4)/params(94))^params(95))/(params(96)+x(10));

    %Synthesis 1
    v(22,1) = params(97)*x(10)/(params(98)+x(10));

    %Synthesis 2
    v(23,1) = params(99)*x(11)/(params(100)+x(11));

    %DAHP synthesis
    v(24,1) = params(101)*x(17)^params(102)*x(10)^params(103)/((params(104)+x(17)^params(102))*(params(105)+x(10)^params(103)));

    %Pyruvate dehydrogenase
    v(25,1) = params(106)*x(11)^params(107)/(params(108)+x(11)^params(107));

    %Methionine synthesis
    v(26,1) = params(109);

    %6-Phosphogluconate dehydrogenase
    v(27,1) = params(110)*x(12)*params(4)/((x(12)+params(111))*(params(4)+params(112)*(1+params(5)/params(113))*(1+params(1)/params(114))));

    %Ribose-phosphate isomerase
    v(28,1) = params(115)*(x(13)-x(16)/params(116));

    %Ribulose-phosphate epimerase
    v(29,1) = params(117)*(x(13)-x(14)/params(118));

    %Ribose phosphate pyrophosphokinase
    v(30,1) = params(119)*x(16)/(params(120)+x(16));

    %Glucose-1-phosphate adenyltransferase
    v(31,1) = params(121)*x(18)*params(1)*(1+(x(4)/params(122))^params(123))/((params(124)+params(1))*(params(125)+x(18)));

    %G6P degradation
    v(32,1) = params(126)*x(2);

    %F6P degradation
    v(33,1) = params(126)*x(3);

    %FDP degradation
    v(34,1) = params(126)*x(4);

    %GAP degradation
    v(35,1) = params(126)*x(5);

    %DHAP degradation
    v(36,1) = params(126)*x(6);

    %PGP degradation
    v(37,1) = params(126)*x(7);

    %PG3 degradation
    v(38,1) = params(126)*x(8);

    %PG2 degradation
    v(39,1) = params(126)*x(9);

    %PEP degradation
    v(40,1) = params(126)*x(10);

    %Pyruvate dilution
    v(41,1) = params(126)*x(11);

    %PG dilution
    v(42,1) = params(126)*x(12);

    %Ribu5P dilution
    v(43,1) = params(126)*x(13);

    %XYL5P dilution
    v(44,1) = params(126)*x(14);

    %SED7P dilution
    v(45,1) = params(126)*x(15);

    %Rib5P dilution
    v(46,1) = params(126)*x(16);

    %E4P dilution
    v(47,1) = params(126)*x(17);

    %GLP dilution
    v(48,1) = params(126)*x(18);
    
end
%{
function params = convertOdeParams(paramVec)
% ParamVec = [bm3; bm4; v2M; v2K; v3M; v3K; v3A; v3alpha; v3beta; v4M; v4K; v4I; v5M; v5K; v53; v54]


    paramVec = paramVec(:);

    bm3 = paramVec(1);
    bm4 = paramVec(2);
    
    params.S = [ 1 -1  0 -1  0   ;
                 0  1 -1  0  0   ;
                 0  0  1  0 -bm3 ;
                 0  0  0  1 -bm4 ;
                 0  0  0  0  1   ];

    params.v0 = 1;
    
    params.v2M = paramVec(3);
    params.v2K = paramVec(4);
    params.v3M = paramVec(5);
    params.v3K = paramVec(6);
    params.v3A = paramVec(7);
    params.v3alpha = paramVec(8);
    params.v3beta = paramVec(9);
    params.v4M = paramVec(10);
    params.v4K = paramVec(11);
    params.v4I = paramVec(12);
    params.v5M = paramVec(13);
    params.v5K = paramVec(14);
    params.v53 = paramVec(15);
    params.v54 = paramVec(16);
    
end
%}