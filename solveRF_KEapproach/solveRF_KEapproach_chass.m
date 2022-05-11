function [RF knownMet res] = solveRF_KEapproach_chass(relative_concMatrix,timeVec,knownKinetics)
% Kinetic equations approach, assuming some or all kinetic equations for
% reactions are known in the system.
% Parameters and kinetics were extracted from Dynamic modeling of the
% central carbon metabolism of Escherichia coli, Chassagnole et al. 2002.

params(1) = 4.27; % catp
params(2) = 0.582; % cadp
params(3) = 0.954783; % camp
params(4) = 0.196759; % cnadp
params(5) = 0.062; % cnadph
params(6) = 1.4644; % cnad
params(7) = 0.0934; % cnadh
params(8) = 2.78e-5; % Dil
params(9) = 110.96; % cfeed
params(10) = 7829.78; % rmaxPTS
params(11) = 3082.3; % KPTSa1
params(12) = 0.01; % KPTSa2
params(13) = 245.3; % KPTSa3
params(14) = 3.66; % nPTSg6p
params(15) = 2.15; % KPTSg6p
params(16) = 650.988; % rmaxPGI
params(17) = 0.1725; % KPGIeq
params(18) = 2.9; % KPGIg6p
params(19) = 0.266; % KPGIf6p
params(20) = 0.2; % KPGIf6ppginh
params(21) = 0.2; % KPGlg6ppginh
params(22) = 0.839824; % rmaxPGM
params(23) = 0.196; % KPGMeq
params(24) = 1.038; % KPGMg6p
params(25) = 0.0136; % KPGMg1p
params(26) = 1.3802; % rmaxG6PDH
params(27) = 14.4; % KG6PDHg6p
params(28) = 6.43; % KG6PDHnadphg6pinh
params(29) = 0.0246; % KG6PDHnadp
params(30) = 0.01; % KG6PDHnadphnadpinh
params(31) = 1840.58; % rmaxPFK
params(32) = 0.123; % KPFKatps
params(33) = 4.14; % KPFKadpc
params(34) = 0.325; % KPFKf6ps
params(35) = 3.26; % KPFKpep
params(36) = 3.89; % KPFKadpb
params(37) = 3.2; % KPFKampb
params(38) = 128; % KPFKadpa
params(39) = 19.1; % KPFKampa
params(40) = 5.62907e6; % LPFK
params(41) = 11.1; % nPFK
params(42) = 10.8716; % rmaxTA
params(43) = 1.05; % KTAeq
params(44) = 9.47338; % rmaxTKa
params(45) = 1.2; % KTKaeq
params(46) = 86.5586; % rmaxTKb
params(47) = 10; % KTKbeq
params(48) = 0.00043711; % rmaxmurSynth
params(49) = 17.4146; % rmaxALDO
params(50) = 0.144; % kALDOeq
params(51) = 1.75; % kALDOfdp
params(52) = 0.088; % kALDOgap
params(53) = 2; % VALDOblf
params(54) = 0.088; % kALDOdhap
params(55) = 0.6; % kALDOgapinh
params(56) = 921.594; % rmaxGAPDH
params(57) = 0.63; % KGAPDHeq
params(58) = 0.683; % KGAPDHgap
params(59) = 1.04e-5; % KGAPDHpgp
params(60) = 0.252; % KGAPDHnad
params(61) = 1.09; % KGAPDHnadh
params(62) = 68.6747; % rmaxTIS
params(63) = 1.39; % kTISeq
params(64) = 2.8; % kTISdhap
params(65) = 0.3; % kTISgap
params(66) = 0.001037; % rmaxTrpSynth
params(67) = 0.0116204; % rmaxG3PDH
params(68) = 1; % KG3PDHdhap
params(69) = 3021.77; % rmaxPGK
params(70) = 1934.4; % KPGKeq
params(71) = 0.185; % KPGKadp
params(72) = 0.653; % KPGKatp
params(73) = 0.0468; % KPGKpgp
params(74) = 0.473; % KPGKpg3
params(75) = 0.0257121; % rmaxSerSynth
params(76) = 1; % KSerSynthpg3
params(77) = 89.0497; % rmaxPGlumu
params(78) = 0.188; % KPGlumueq
params(79) = 0.2; % KPGlumupg3
params(80) = 0.369; % KPGlumupg2
params(81) = 330.448; % rmaxENO
params(82) = 6.73; % KENOeq
params(83) = 0.1; % KENOpg2
params(84) = 0.135; % KENOpep
params(85) = 0.0611315; % rmaxPK
params(86) = 0.31; % KPKpep
params(87) = 4; % nPK
params(88) = 1000; % LPK
params(89) = 22.5; % KPKatp
params(90) = 0.19; % KPKfdp
params(91) = 0.2; % KPKamp
params(92) = 0.26; % KPKadp
params(93) = 0.107021; % rmaxpepCxylase
params(94) = 0.7; % KpepCxylasefdp
params(95) = 4.21; % npepCxylasefdp
params(96) = 4.07; % KpepCxylasepep
params(97) = 0.019539; % rmaxSynth1
params(98) = 1; % KSynth1pep
params(99) = 0.0736186; % rmaxSynth2
params(100) = 1; % KSynth2pyr
params(101) = 0.107953; % rmaxDAHPS
params(102) = 2.6; % nDAHPSe4p
params(103) = 2.2; % nDAHPSpep
params(104) = 0.035; % KDAHPSe4p
params(105) = 0.0053; % KDAHPSpep
params(106) = 6.05953; % rmaxPDH
params(107) = 3.68; % nPDH
params(108) = 1159; % KPDHpyr
params(109) = 0.0022627; % rmaxMetSynth
params(110) = 16.2324; % rmaxPGDH
params(111) = 37.5; % KPGDHpg
params(112) = 0.0506; % KPGDHnadp
params(113) = 0.0138; % KPGDHnadphinh
params(114) = 208; % KPGDHatpinh
params(115) = 4.83841; % rmaxR5PI
params(116) = 4; % KR5PIeq
params(117) = 6.73903; % rmaxRu5P
params(118) = 1.4; % KRu5Peq
params(119) = 0.0129005; % rmaxRPPK
params(120) = 0.1; % KRPPKrib5p
params(121) = 0.00752546; % rmaxG1PAT
params(122) = 0.119; % KG1PATfdp
params(123) = 1.2; % nG1PATfdp
params(124) = 4.42; % KG1PATatp
params(125) = 3.2; % KG1PATg1p
params(126) = 2.78e-5; % mu

for t = 2:size(relative_concMatrix,1)
    relative_Vpool(t-1,:) = (relative_concMatrix(t,:)-relative_concMatrix(t-1,:))./(timeVec(t)-timeVec(t-1));
end

knownMet = findKnownMet(knownKinetics);

minRF = max(relative_concMatrix(:,knownMet))/50;

lb = minRF;
ub = 1000*ones(length(knownMet),1);
x0 = (1000-1).*rand(1,length(knownMet)) + 1;

if ~isempty(knownMet)
    [RF,res] = lsqnonlin(@(RF) dxdt_functions(RF,relative_Vpool,relative_concMatrix,params,knownKinetics,knownMet),x0,lb,ub);
else
    RF = 0;
    knownMet = 0;
end

end


function func = dxdt_functions(RF,relative_Vpool,relative_concMatrix,params,knownKinetics,knownMet)

nT = size(relative_Vpool,1);
RFnum = zeros(size(relative_Vpool,2),1);
RFnum(knownMet) = 1:length(knownMet);

func = [];
for t = 2:nT
    
    % dxdt1
    if all(ismember([1 2],knownKinetics))
        func(end+1) = ...
            abs(params(8)*(params(9)-(relative_concMatrix(t,1)/RF(RFnum(1))))...
            -params(10)*(relative_concMatrix(t,1)/RF(RFnum(1)))*((relative_concMatrix(t,10)/RF(RFnum(10)))/(relative_concMatrix(t,11)/RF(RFnum(11))))/((params(11)+params(12)*((relative_concMatrix(t,10)/RF(RFnum(10)))/(relative_concMatrix(t,11)/RF(RFnum(11))))+params(13)*(relative_concMatrix(t,1)/RF(RFnum(1)))+(relative_concMatrix(t,1)/RF(RFnum(1)))*((relative_concMatrix(t,10)/RF(RFnum(10)))/(relative_concMatrix(t,11)/RF(RFnum(11)))))*(1+(relative_concMatrix(t,2)/RF(RFnum(2)))^params(14)/params(15)))...
            -relative_Vpool(t,1)/RF(RFnum(1)));
    end
    
    % dxdt2 
    if all(ismember([2 3 4 5],knownKinetics)) && max(knownKinetics) >= 32
        func(end+1) = ...
            abs(65 * (params(10)*(relative_concMatrix(t,1)/RF(RFnum(1)))*((relative_concMatrix(t,10)/RF(RFnum(10)))/(relative_concMatrix(t,11)/RF(RFnum(11))))/((params(11)+params(12)*((relative_concMatrix(t,10)/RF(RFnum(10)))/(relative_concMatrix(t,11)/RF(RFnum(11))))+params(13)*(relative_concMatrix(t,1)/RF(RFnum(1)))+(relative_concMatrix(t,1)/RF(RFnum(1)))*((relative_concMatrix(t,10)/RF(RFnum(10)))/(relative_concMatrix(t,11)/RF(RFnum(11)))))*(1+(relative_concMatrix(t,2)/RF(RFnum(2)))^params(14)/params(15))))...
            -1*(params(16)*((relative_concMatrix(t,2)/RF(RFnum(2)))-(relative_concMatrix(t,3)/RF(RFnum(3)))/params(17))/(params(18)*(1+(relative_concMatrix(t,3)/RF(RFnum(3)))/(params(19)*(1+(relative_concMatrix(t,12)/RF(RFnum(12)))/params(20)))+(relative_concMatrix(t,12)/RF(RFnum(12)))/params(21))+(relative_concMatrix(t,2)/RF(RFnum(2)))))...
            -1*(params(22)*((relative_concMatrix(t,2)/RF(RFnum(2)))-(relative_concMatrix(t,18)/RF(RFnum(18)))/params(23))/(params(24)*(1+(relative_concMatrix(t,18)/RF(RFnum(18)))/params(25))+(relative_concMatrix(t,2)/RF(RFnum(2)))))...
            -1*(params(26)*(relative_concMatrix(t,2)/RF(RFnum(2)))*params(4)/(((relative_concMatrix(t,2)/RF(RFnum(2)))+params(27))*(1+params(5)/params(28))*(params(29)*(1+params(5)/params(30))+params(4))))...
            -1*(params(126)*(relative_concMatrix(t,2)/RF(RFnum(2))))...
            -relative_Vpool(t,2)/RF(RFnum(2)));
    end
        
    % dxdt3
    if all(ismember([3 6 7 9 10],knownKinetics)) && max(knownKinetics) >= 32
        func(end+1) = ...
            abs((params(16)*((relative_concMatrix(t,2)/RF(RFnum(2)))-(relative_concMatrix(t,3)/RF(RFnum(3)))/params(17))/(params(18)*(1+(relative_concMatrix(t,3)/RF(RFnum(3)))/(params(19)*(1+(relative_concMatrix(t,12)/RF(RFnum(12)))/params(20)))+(relative_concMatrix(t,12)/RF(RFnum(12)))/params(21))+(relative_concMatrix(t,2)/RF(RFnum(2)))))...
            -1*(params(31)*params(1)*(relative_concMatrix(t,3)/RF(RFnum(3)))/((params(1)+params(32)*(1+params(2)/params(33)))*((relative_concMatrix(t,3)/RF(RFnum(3)))+params(34)*(1+(relative_concMatrix(t,10)/RF(RFnum(10)))/params(35)+params(2)/params(36)+params(3)/params(37))/(1+params(2)/params(38)+params(3)/params(39)))*(1+params(40)/(1+(relative_concMatrix(t,3)/RF(RFnum(3)))*(1+params(2)/params(38)+params(3)/params(39))/(params(34)*(1+(relative_concMatrix(t,10)/RF(RFnum(10)))/params(35)+params(2)/params(36)+params(3)/params(37))))^params(41))))...
            +1*(params(42)*((relative_concMatrix(t,5)/RF(RFnum(5)))*(relative_concMatrix(t,15)/RF(RFnum(15)))-(relative_concMatrix(t,17)/RF(RFnum(17)))*(relative_concMatrix(t,3)/RF(RFnum(3)))/params(43)))...
            +1*(params(46)*((relative_concMatrix(t,14)/RF(RFnum(14)))*(relative_concMatrix(t,17)/RF(RFnum(17)))-(relative_concMatrix(t,3)/RF(RFnum(3)))*(relative_concMatrix(t,5)/RF(RFnum(5)))/params(47)))...
            -2*(params(48))...
            -1*(params(126)*(relative_concMatrix(t,3)/RF(RFnum(3))))...
            -relative_Vpool(t,3)/RF(RFnum(3)));
    end
     
    % dxdt4    
    if all(ismember([6 11],knownKinetics)) && max(knownKinetics) >= 32
        func(end+1) = ...
            abs((params(31)*params(1)*(relative_concMatrix(t,3)/RF(RFnum(3)))/((params(1)+params(32)*(1+params(2)/params(33)))*((relative_concMatrix(t,3)/RF(RFnum(3)))+params(34)*(1+(relative_concMatrix(t,10)/RF(RFnum(10)))/params(35)+params(2)/params(36)+params(3)/params(37))/(1+params(2)/params(38)+params(3)/params(39)))*(1+params(40)/(1+(relative_concMatrix(t,3)/RF(RFnum(3)))*(1+params(2)/params(38)+params(3)/params(39))/(params(34)*(1+(relative_concMatrix(t,10)/RF(RFnum(10)))/params(35)+params(2)/params(36)+params(3)/params(37))))^params(41))))...
            -1*(params(49)*((relative_concMatrix(t,4)/RF(RFnum(4)))-(relative_concMatrix(t,5)/RF(RFnum(5)))*(relative_concMatrix(t,6)/RF(RFnum(6)))/params(50))/(params(51)+(relative_concMatrix(t,4)/RF(RFnum(4)))+params(52)*(relative_concMatrix(t,6)/RF(RFnum(6)))/(params(50)*params(53))+params(54)*(relative_concMatrix(t,5)/RF(RFnum(5)))/(params(50)*params(53))+(relative_concMatrix(t,4)/RF(RFnum(4)))*(relative_concMatrix(t,5)/RF(RFnum(5)))/params(55)+(relative_concMatrix(t,5)/RF(RFnum(5)))*(relative_concMatrix(t,6)/RF(RFnum(6)))/(params(53)*params(50))))...
            -1*(params(126)*(relative_concMatrix(t,4)/RF(RFnum(4))))...
            -relative_Vpool(t,4)/RF(RFnum(4)));
    end

    % dxdt5
    if all(ismember([7 8 9 11 12 13 14],knownKinetics)) && max(knownKinetics) >= 32
        func(end+1) = ...
            abs(-1*(params(42)*((relative_concMatrix(t,5)/RF(RFnum(5)))*(relative_concMatrix(t,15)/RF(RFnum(15)))-(relative_concMatrix(t,17)/RF(RFnum(17)))*(relative_concMatrix(t,3)/RF(RFnum(3)))/params(43)))...
            +1*(params(44)*((relative_concMatrix(t,16)/RF(RFnum(16)))*(relative_concMatrix(t,14)/RF(RFnum(14)))-(relative_concMatrix(t,15)/RF(RFnum(15)))*(relative_concMatrix(t,5)/RF(RFnum(5)))/params(45)))...
            +1*(params(46)*((relative_concMatrix(t,14)/RF(RFnum(14)))*(relative_concMatrix(t,17)/RF(RFnum(17)))-(relative_concMatrix(t,3)/RF(RFnum(3)))*(relative_concMatrix(t,5)/RF(RFnum(5)))/params(47)))...
            +1*(params(49)*((relative_concMatrix(t,4)/RF(RFnum(4)))-(relative_concMatrix(t,5)/RF(RFnum(5)))*(relative_concMatrix(t,6)/RF(RFnum(6)))/params(50))/(params(51)+(relative_concMatrix(t,4)/RF(RFnum(4)))+params(52)*(relative_concMatrix(t,6)/RF(RFnum(6)))/(params(50)*params(53))+params(54)*(relative_concMatrix(t,5)/RF(RFnum(5)))/(params(50)*params(53))+(relative_concMatrix(t,4)/RF(RFnum(4)))*(relative_concMatrix(t,5)/RF(RFnum(5)))/params(55)+(relative_concMatrix(t,5)/RF(RFnum(5)))*(relative_concMatrix(t,6)/RF(RFnum(6)))/(params(53)*params(50))))...
            -1*(params(56)*((relative_concMatrix(t,5)/RF(RFnum(5)))*params(6)-(relative_concMatrix(t,7)/RF(RFnum(7)))*params(7)/params(57))/((params(58)*(1+(relative_concMatrix(t,7)/RF(RFnum(7)))/params(59))+(relative_concMatrix(t,5)/RF(RFnum(5))))*(params(60)*(1+params(7)/params(61))+params(6))))...
            +1*(params(62)*((relative_concMatrix(t,6)/RF(RFnum(6)))-(relative_concMatrix(t,5)/RF(RFnum(5)))/params(63))/(params(64)*(1+(relative_concMatrix(t,5)/RF(RFnum(5)))/params(65))+(relative_concMatrix(t,6)/RF(RFnum(6)))))...
            +1*(params(66))...
            -1*(params(126)*(relative_concMatrix(t,5)/RF(RFnum(5))))...
            -relative_Vpool(t,5)/RF(RFnum(5)));
    end
    
    % dxdt6
    if all(ismember([11 13 15],knownKinetics)) && max(knownKinetics) >= 32
        func(end+1) = ...
            abs((params(49)*((relative_concMatrix(t,4)/RF(RFnum(4)))-(relative_concMatrix(t,5)/RF(RFnum(5)))*(relative_concMatrix(t,6)/RF(RFnum(6)))/params(50))/(params(51)+(relative_concMatrix(t,4)/RF(RFnum(4)))+params(52)*(relative_concMatrix(t,6)/RF(RFnum(6)))/(params(50)*params(53))+params(54)*(relative_concMatrix(t,5)/RF(RFnum(5)))/(params(50)*params(53))+(relative_concMatrix(t,4)/RF(RFnum(4)))*(relative_concMatrix(t,5)/RF(RFnum(5)))/params(55)+(relative_concMatrix(t,5)/RF(RFnum(5)))*(relative_concMatrix(t,6)/RF(RFnum(6)))/(params(53)*params(50))))...
            -1*(params(62)*((relative_concMatrix(t,6)/RF(RFnum(6)))-(relative_concMatrix(t,5)/RF(RFnum(5)))/params(63))/(params(64)*(1+(relative_concMatrix(t,5)/RF(RFnum(5)))/params(65))+(relative_concMatrix(t,6)/RF(RFnum(6)))))...
            -1*(params(67)*(relative_concMatrix(t,6)/RF(RFnum(6)))/(params(68)+(relative_concMatrix(t,6)/RF(RFnum(6)))))...
            -1*(params(126)*(relative_concMatrix(t,6)/RF(RFnum(6))))...
            -relative_Vpool(t,6)/RF(RFnum(6)));
    end

    % dxdt7
    if all(ismember([12 16],knownKinetics)) && max(knownKinetics) >= 32
        func(end+1) = ...
            abs((params(56)*((relative_concMatrix(t,5)/RF(RFnum(5)))*params(6)-(relative_concMatrix(t,7)/RF(RFnum(7)))*params(7)/params(57))/((params(58)*(1+(relative_concMatrix(t,7)/RF(RFnum(7)))/params(59))+(relative_concMatrix(t,5)/RF(RFnum(5))))*(params(60)*(1+params(7)/params(61))+params(6))))...
            -1*(params(69)*(params(2)*(relative_concMatrix(t,7)/RF(RFnum(7)))-params(1)*(relative_concMatrix(t,8)/RF(RFnum(8)))/params(70))/((params(71)*(1+params(1)/params(72))+params(2))*(params(73)*(1+(relative_concMatrix(t,8)/RF(RFnum(8)))/params(74))+(relative_concMatrix(t,7)/RF(RFnum(7))))))...
            -1*(params(126)*(relative_concMatrix(t,7)/RF(RFnum(7))))...
            -relative_Vpool(t,7)/RF(RFnum(7)));
        knownMet = [knownMet 5 7 8];
    end
    
    % dxdt8
    if all(ismember([16 17 18],knownKinetics)) && max(knownKinetics) >= 32
        func(end+1) = ...
            abs((params(69)*(params(2)*(relative_concMatrix(t,7)/RF(RFnum(7)))-params(1)*(relative_concMatrix(t,8)/RF(RFnum(8)))/params(70))/((params(71)*(1+params(1)/params(72))+params(2))*(params(73)*(1+(relative_concMatrix(t,8)/RF(RFnum(8)))/params(74))+(relative_concMatrix(t,7)/RF(RFnum(7))))))...
            -1*(params(75)*(relative_concMatrix(t,8)/RF(RFnum(8)))/(params(76)+(relative_concMatrix(t,8)/RF(RFnum(8)))))...
            -1*(params(77)*((relative_concMatrix(t,8)/RF(RFnum(8)))-(relative_concMatrix(t,9)/RF(RFnum(9)))/params(78))/(params(79)*(1+(relative_concMatrix(t,9)/RF(RFnum(9)))/params(80))+(relative_concMatrix(t,8)/RF(RFnum(8)))))...
            -1*(params(126)*(relative_concMatrix(t,9)/RF(RFnum(9))))...
            -relative_Vpool(t,8)/RF(RFnum(8)));
    end
    
    % dxdt9
    if all(ismember([18 19],knownKinetics)) && max(knownKinetics) >= 32
        func(end+1) = ...
            abs((params(77)*((relative_concMatrix(t,8)/RF(RFnum(8)))-(relative_concMatrix(t,9)/RF(RFnum(9)))/params(78))/(params(79)*(1+(relative_concMatrix(t,9)/RF(RFnum(9)))/params(80))+(relative_concMatrix(t,8)/RF(RFnum(8)))))...
            -1*(params(81)*((relative_concMatrix(t,9)/RF(RFnum(9)))-(relative_concMatrix(t,10)/RF(RFnum(10)))/params(82))/(params(83)*(1+(relative_concMatrix(t,10)/RF(RFnum(10)))/params(84))+(relative_concMatrix(t,9)/RF(RFnum(9)))))...
            -1*(params(126)*(relative_concMatrix(t,10)/RF(RFnum(10))))...
            -relative_Vpool(t,9)/RF(RFnum(9)));
    end
    
    % dxdt10
    if all(ismember([2 19 20 21 22 24],knownKinetics)) && max(knownKinetics) >= 32
        func(end+1) = ...
            abs(-65*(params(10)*(relative_concMatrix(t,1)/RF(RFnum(1)))*((relative_concMatrix(t,10)/RF(RFnum(10)))/(relative_concMatrix(t,11)/RF(RFnum(11))))/((params(11)+params(12)*((relative_concMatrix(t,10)/RF(RFnum(10)))/(relative_concMatrix(t,11)/RF(RFnum(11))))+params(13)*(relative_concMatrix(t,1)/RF(RFnum(1)))+(relative_concMatrix(t,1)/RF(RFnum(1)))*((relative_concMatrix(t,10)/RF(RFnum(10)))/(relative_concMatrix(t,11)/RF(RFnum(11)))))*(1+(relative_concMatrix(t,2)/RF(RFnum(2)))^params(14)/params(15))))...
            +1*(params(81)*((relative_concMatrix(t,9)/RF(RFnum(9)))-(relative_concMatrix(t,10)/RF(RFnum(10)))/params(82))/(params(83)*(1+(relative_concMatrix(t,10)/RF(RFnum(10)))/params(84))+(relative_concMatrix(t,9)/RF(RFnum(9)))))...
            -1*(params(85)*(relative_concMatrix(t,10)/RF(RFnum(10)))*((relative_concMatrix(t,10)/RF(RFnum(10)))/params(86)+1)^(params(87)-1)*params(2)/(params(86)*(params(88)*((1+params(1)/params(89))/((relative_concMatrix(t,4)/RF(RFnum(4)))/params(90)+params(3)/params(91)+1))^params(87)+((relative_concMatrix(t,10)/RF(RFnum(10)))/params(86)+1)^params(87))*(params(2)+params(92))))...
            -1*(params(93)*(relative_concMatrix(t,10)/RF(RFnum(10)))*(1+((relative_concMatrix(t,4)/RF(RFnum(4)))/params(94))^params(95))/(params(96)+(relative_concMatrix(t,10)/RF(RFnum(10)))))...
            -1*(params(97)*(relative_concMatrix(t,10)/RF(RFnum(10)))/(params(98)+(relative_concMatrix(t,10)/RF(RFnum(10)))))...
            -1*(params(101)*(relative_concMatrix(t,17)/RF(RFnum(17)))^params(102)*(relative_concMatrix(t,10)/RF(RFnum(10)))^params(103)/((params(104)+(relative_concMatrix(t,17)/RF(RFnum(17)))^params(102))*(params(105)+(relative_concMatrix(t,10)/RF(RFnum(10)))^params(103))))...
            -1*(params(126)*(relative_concMatrix(t,10)/RF(RFnum(10))))...
            -relative_Vpool(t,10)/RF(RFnum(10)));
    end
    
    % dxdt11
    if all(ismember([2 14 20 23 25 26],knownKinetics)) && max(knownKinetics) >= 32
        func(end+1) = ...
            abs(65*(params(10)*(relative_concMatrix(t,1)/RF(RFnum(1)))*((relative_concMatrix(t,10)/RF(RFnum(10)))/(relative_concMatrix(t,11)/RF(RFnum(11))))/((params(11)+params(12)*((relative_concMatrix(t,10)/RF(RFnum(10)))/(relative_concMatrix(t,11)/RF(RFnum(11))))+params(13)*(relative_concMatrix(t,1)/RF(RFnum(1)))+(relative_concMatrix(t,1)/RF(RFnum(1)))*((relative_concMatrix(t,10)/RF(RFnum(10)))/(relative_concMatrix(t,11)/RF(RFnum(11)))))*(1+(relative_concMatrix(t,2)/RF(RFnum(2)))^params(14)/params(15))))...
            +1*(params(66))...
            +1*(params(85)*(relative_concMatrix(t,10)/RF(RFnum(10)))*((relative_concMatrix(t,10)/RF(RFnum(10)))/params(86)+1)^(params(87)-1)*params(2)/(params(86)*(params(88)*((1+params(1)/params(89))/((relative_concMatrix(t,4)/RF(RFnum(4)))/params(90)+params(3)/params(91)+1))^params(87)+((relative_concMatrix(t,10)/RF(RFnum(10)))/params(86)+1)^params(87))*(params(2)+params(92))))...
            -1*(params(99)*(relative_concMatrix(t,11)/RF(RFnum(11)))/(params(100)+(relative_concMatrix(t,11)/RF(RFnum(11)))))...
            -1*(params(106)*(relative_concMatrix(t,11)/RF(RFnum(11)))^params(107)/(params(108)+(relative_concMatrix(t,11)/RF(RFnum(11)))^params(107)))...
            +1*(params(109))...
            -1*(params(126)*(relative_concMatrix(t,11)/RF(RFnum(11))))...
            -relative_Vpool(t,11)/RF(RFnum(11)));
    end
    
    % dxdt12
    if all(ismember([5 27],knownKinetics)) && max(knownKinetics) >= 32
        func(end+1) = ...
            abs((params(26)*(relative_concMatrix(t,2)/RF(RFnum(2)))*params(4)/(((relative_concMatrix(t,2)/RF(RFnum(2)))+params(27))*(1+params(5)/params(28))*(params(29)*(1+params(5)/params(30))+params(4))))...
            -1*(params(110)*(relative_concMatrix(t,12)/RF(RFnum(12)))*params(4)/(((relative_concMatrix(t,12)/RF(RFnum(12)))+params(111))*(params(4)+params(112)*(1+params(5)/params(113))*(1+params(1)/params(114)))))...
            -1*(params(126)*(relative_concMatrix(t,12)/RF(RFnum(12))))...
            -relative_Vpool(t,12)/RF(RFnum(12)));
    end
    
    % dxdt13
    if all(ismember([27 28 29],knownKinetics)) && max(knownKinetics) >= 32
        func(end+1) = ...
            abs((params(110)*(relative_concMatrix(t,12)/RF(RFnum(12)))*params(4)/(((relative_concMatrix(t,12)/RF(RFnum(12)))+params(111))*(params(4)+params(112)*(1+params(5)/params(113))*(1+params(1)/params(114)))))...
            -1*(params(115)*((relative_concMatrix(t,13)/RF(RFnum(13)))-(relative_concMatrix(t,16)/RF(RFnum(16)))/params(116)))...
            -1*(params(117)*((relative_concMatrix(t,13)/RF(RFnum(13)))-(relative_concMatrix(t,14)/RF(RFnum(14)))/params(118)))...
            -1*(params(126)*(relative_concMatrix(t,13)/RF(RFnum(13))))...
            -relative_Vpool(t,13)/RF(RFnum(13)));
    end
    
    % dxdt14
    if all(ismember([8 9 29],knownKinetics)) && max(knownKinetics) >= 32
        func(end+1) = ...
            abs(-1*(params(44)*((relative_concMatrix(t,16)/RF(RFnum(16)))*(relative_concMatrix(t,14)/RF(RFnum(14)))-(relative_concMatrix(t,15)/RF(RFnum(15)))*(relative_concMatrix(t,5)/RF(RFnum(5)))/params(45)))...
            -1*(params(46)*((relative_concMatrix(t,14)/RF(RFnum(14)))*(relative_concMatrix(t,17)/RF(RFnum(17)))-(relative_concMatrix(t,3)/RF(RFnum(3)))*(relative_concMatrix(t,5)/RF(RFnum(5)))/params(47)))...
            +1*(params(117)*((relative_concMatrix(t,13)/RF(RFnum(13)))-(relative_concMatrix(t,14)/RF(RFnum(14)))/params(118)))...
            -1*(params(126)*(relative_concMatrix(t,14)/RF(RFnum(14))))...
            -relative_Vpool(t,14)/RF(RFnum(14)));
    end
    
    % dxdt15
    if all(ismember([7 8],knownKinetics)) && max(knownKinetics) >= 32
        func(end+1) = ...
            abs(-1*(params(42)*((relative_concMatrix(t,5)/RF(RFnum(5)))*(relative_concMatrix(t,15)/RF(RFnum(15)))-(relative_concMatrix(t,17)/RF(RFnum(17)))*(relative_concMatrix(t,3)/RF(RFnum(3)))/params(43)))...
            +1*(params(44)*((relative_concMatrix(t,16)/RF(RFnum(16)))*(relative_concMatrix(t,14)/RF(RFnum(14)))-(relative_concMatrix(t,15)/RF(RFnum(15)))*(relative_concMatrix(t,5)/RF(RFnum(5)))/params(45)))...
            -1*(params(126)*(relative_concMatrix(t,15)/RF(RFnum(15))))...
            -relative_Vpool(t,15)/RF(RFnum(15)));
    end
    
    % dxdt16
    if all(ismember([8 28 30],knownKinetics)) && max(knownKinetics) >= 32
        func(end+1) = ...
            abs(-1*(params(44)*((relative_concMatrix(t,16)/RF(RFnum(16)))*(relative_concMatrix(t,14)/RF(RFnum(14)))-(relative_concMatrix(t,15)/RF(RFnum(15)))*(relative_concMatrix(t,5)/RF(RFnum(5)))/params(45)))...
            +1*(params(115)*((relative_concMatrix(t,13)/RF(RFnum(13)))-(relative_concMatrix(t,16)/RF(RFnum(16)))/params(116)))...
            -1*(params(119)*(relative_concMatrix(t,16)/RF(RFnum(16)))/(params(120)+(relative_concMatrix(t,16)/RF(RFnum(16)))))...
            -1*(params(126)*(relative_concMatrix(t,16)/RF(RFnum(16))))...
            -relative_Vpool(t,16)/RF(RFnum(16)));
    end
    
    % dxdt17
    if all(ismember([7 9 24],knownKinetics)) && max(knownKinetics) >= 32
        func(end+1) = ...
            abs((params(42)*((relative_concMatrix(t,5)/RF(RFnum(5)))*(relative_concMatrix(t,15)/RF(RFnum(15)))-(relative_concMatrix(t,17)/RF(RFnum(17)))*(relative_concMatrix(t,3)/RF(RFnum(3)))/params(43)))...
            -1*(params(46)*((relative_concMatrix(t,14)/RF(RFnum(14)))*(relative_concMatrix(t,17)/RF(RFnum(17)))-(relative_concMatrix(t,3)/RF(RFnum(3)))*(relative_concMatrix(t,5)/RF(RFnum(5)))/params(47)))...
            -1*(params(101)*(relative_concMatrix(t,17)/RF(RFnum(17)))^params(102)*(relative_concMatrix(t,10)/RF(RFnum(10)))^params(103)/((params(104)+(relative_concMatrix(t,17)/RF(RFnum(17)))^params(102))*(params(105)+(relative_concMatrix(t,10)/RF(RFnum(10)))^params(103))))...
            -1*(params(126)*(relative_concMatrix(t,17)/RF(RFnum(17))))...
            -relative_Vpool(t,17)/RF(RFnum(17)));
    end
    
    % dxdt18
    if all(ismember([4 31],knownKinetics)) && max(knownKinetics) >= 32
        func(end+1) = ...
            abs((params(22)*((relative_concMatrix(t,2)/RF(RFnum(2)))-(relative_concMatrix(t,18)/RF(RFnum(18)))/params(23))/(params(24)*(1+(relative_concMatrix(t,18)/RF(RFnum(18)))/params(25))+(relative_concMatrix(t,2)/RF(RFnum(2)))))...
            -1*(params(121)*(relative_concMatrix(t,18)/RF(RFnum(18)))*params(1)*(1+((relative_concMatrix(t,4)/RF(RFnum(4)))/params(122))^params(123))/((params(124)+params(1))*(params(125)+(relative_concMatrix(t,18)/RF(RFnum(18))))))...
            -1*(params(126)*(relative_concMatrix(t,18)/RF(RFnum(18))))...
            -relative_Vpool(t,18)/RF(RFnum(18)));
    end
    
    if isempty(func)
        func = 0;
    end
end

end

function knownMet = findKnownMet(knownKinetics)
    knownMet = [];

    % dxdt1
    if all(ismember([1 2],knownKinetics))
        knownMet = [knownMet 1 2 10 11];
    end
    
    % dxdt2 
    if all(ismember([2 3 4 5],knownKinetics)) && max(knownKinetics) >= 32
        knownMet = [knownMet 1 2 3 10 11 12 18];
    end
        
    % dxdt3
    if all(ismember([3 6 7 9 10],knownKinetics)) && max(knownKinetics) >= 32
        knownMet = [knownMet 2 3 5 10 12 14 15 17];
    end
     
    % dxdt4    
    if all(ismember([6 11],knownKinetics)) && max(knownKinetics) >= 32
        knownMet = [knownMet 3 4 5 6 10];
    end

    % dxdt5
    if all(ismember([7 8 9 11 12 13 14],knownKinetics)) && max(knownKinetics) >= 32
        knownMet = [knownMet 3 4 5 6 7 14 15 16 17];
    end
    
    % dxdt6
    if all(ismember([11 13 15],knownKinetics)) && max(knownKinetics) >= 32
        knownMet = [knownMet 4 5 6];
    end

    % dxdt7
    if all(ismember([12 16],knownKinetics)) && max(knownKinetics) >= 32
        knownMet = [knownMet 5 7 8];
    end
    
    % dxdt8
    if all(ismember([16 17 18],knownKinetics)) && max(knownKinetics) >= 32
        knownMet = [knownMet 7 8 9];
    end
    
    % dxdt9
    if all(ismember([18 19],knownKinetics)) && max(knownKinetics) >= 32
        knownMet = [knownMet 8 9 10];
    end
    
    % dxdt10
    if all(ismember([2 19 20 21 22 24],knownKinetics)) && max(knownKinetics) >= 32
        knownMet = [knownMet 1 4 9 10 11 17];
    end
    
    % dxdt11
    if all(ismember([2 14 20 23 25 26],knownKinetics)) && max(knownKinetics) >= 32
        knownMet = [knownMet 1 2 4 10 11];
    end
    
    % dxdt12
    if all(ismember([5 27],knownKinetics)) && max(knownKinetics) >= 32
        knownMet = [knownMet 2 12];
    end
    
    % dxdt13
    if all(ismember([27 28 29],knownKinetics)) && max(knownKinetics) >= 32
        knownMet = [knownMet 12 13 14 16];
    end
    
    % dxdt14
    if all(ismember([8 9 29],knownKinetics)) && max(knownKinetics) >= 32
        knownMet = [knownMet 3 5 13 14 15 16 17];
    end
    
    % dxdt15
    if all(ismember([7 8],knownKinetics)) && max(knownKinetics) >= 32
        knownMet = [knownMet 3 5 14 15 16 17];
    end
    
    % dxdt16
    if all(ismember([8 28 30],knownKinetics)) && max(knownKinetics) >= 32
        knownMet = [knownMet 5 13 14 15 16];
    end
    
    % dxdt17
    if all(ismember([7 9 24],knownKinetics)) && max(knownKinetics) >= 32
        knownMet = [knownMet 3 5 10 14 15 17];
    end
    
    % dxdt18
    if all(ismember([4 31],knownKinetics)) && max(knownKinetics) >= 32
        knownMet = [knownMet 2 4 18];
    end    
    knownMet = unique(knownMet);
end