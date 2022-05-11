function flux = calcFluxWithKinetics_chass(concMatrix,timeVec,knownKinetics)
% Calculate flux reaction rates in E. coli system using known kinetic
% equations.
% Parameters and kinetics were extracted from Dynamic modeling of the
% central carbon metabolism of Escherichia coli, Chassagnole et al. 2002.

nT = size(concMatrix,1)-1;
deltaT = timeVec(2)-timeVec(1);

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
    
flux = nan(nT,66);
for t = 1:nT% List of reactions
    %Extracellular glucose kinetics
    if ismember(1,knownKinetics)
        flux(t,1) = params(8)*(params(9)-concMatrix(t,1));    
    end
    
    %Phosphotransferase system
    if ismember(2,knownKinetics)
        flux(t,2) = params(10)*concMatrix(t,1)*(concMatrix(t,10)/concMatrix(t,11))/((params(11)+params(12)*(concMatrix(t,10)/concMatrix(t,11))+params(13)*concMatrix(t,1)+concMatrix(t,1)*(concMatrix(t,10)/concMatrix(t,11)))*(1+concMatrix(t,2)^params(14)/params(15)));
    end
    
    %Glucose-6-phosphate isomerase
    if ismember(3,knownKinetics)
        flux(t,3) = params(16)*(concMatrix(t,2)-concMatrix(t,3)/params(17))/(params(18)*(1+concMatrix(t,3)/(params(19)*(1+concMatrix(t,12)/params(20)))+concMatrix(t,12)/params(21))+concMatrix(t,2));
    end
    
    %Phosphoglucomutase
    if ismember(4,knownKinetics)
        flux(t,4) = params(22)*(concMatrix(t,2)-concMatrix(t,18)/params(23))/(params(24)*(1+concMatrix(t,18)/params(25))+concMatrix(t,2));
    end
    
    %Glucose-6-phosphate dehydrogenase
    if ismember(5,knownKinetics)
        flux(t,5) = params(26)*concMatrix(t,2)*params(4)/((concMatrix(t,2)+params(27))*(1+params(5)/params(28))*(params(29)*(1+params(5)/params(30))+params(4)));
    end
    
    %Phosphofructokinase
    if ismember(6,knownKinetics)
        flux(t,6) = params(31)*params(1)*concMatrix(t,3)/((params(1)+params(32)*(1+params(2)/params(33)))*(concMatrix(t,3)+params(34)*(1+concMatrix(t,10)/params(35)+params(2)/params(36)+params(3)/params(37))/(1+params(2)/params(38)+params(3)/params(39)))*(1+params(40)/(1+concMatrix(t,3)*(1+params(2)/params(38)+params(3)/params(39))/(params(34)*(1+concMatrix(t,10)/params(35)+params(2)/params(36)+params(3)/params(37))))^params(41)));
    end
    
    %Transaldolase
    if ismember(7,knownKinetics)
        flux(t,7) = params(42)*(concMatrix(t,5)*concMatrix(t,15)-concMatrix(t,17)*concMatrix(t,3)/params(43));
    end
    
    %Transketolase a
    if ismember(8,knownKinetics)
        flux(t,8) = params(44)*(concMatrix(t,16)*concMatrix(t,14)-concMatrix(t,15)*concMatrix(t,5)/params(45));
    end
    
    %Transketolase b
    if ismember(9,knownKinetics)
        flux(t,9) = params(46)*(concMatrix(t,14)*concMatrix(t,17)-concMatrix(t,3)*concMatrix(t,5)/params(47));
    end
    
    %Mureine synthesis
    if ismember(10,knownKinetics)
        flux(t,10) = params(48);
    end
    
    %Aldolase
    if ismember(11,knownKinetics)
        flux(t,11) = params(49)*(concMatrix(t,4)-concMatrix(t,5)*concMatrix(t,6)/params(50))/(params(51)+concMatrix(t,4)+params(52)*concMatrix(t,6)/(params(50)*params(53))+params(54)*concMatrix(t,5)/(params(50)*params(53))+concMatrix(t,4)*concMatrix(t,5)/params(55)+concMatrix(t,5)*concMatrix(t,6)/(params(53)*params(50)));
    end
    
    %Glyceraldehyde-3-phosphate dehydrogenase
    if ismember(12,knownKinetics)
        flux(t,12) = params(56)*(concMatrix(t,5)*params(6)-concMatrix(t,7)*params(7)/params(57))/((params(58)*(1+concMatrix(t,7)/params(59))+concMatrix(t,5))*(params(60)*(1+params(7)/params(61))+params(6)));
    end
    
    %Triosephosphate isomerase
    if ismember(13,knownKinetics)
        flux(t,13) = params(62)*(concMatrix(t,6)-concMatrix(t,5)/params(63))/(params(64)*(1+concMatrix(t,5)/params(65))+concMatrix(t,6));
    end
    
    %Tryptophan synthesis
    if ismember(14,knownKinetics)
        flux(t,14) = params(66);
    end
    
    %Glycerol-3-phosphate dehydrogenase
    if ismember(15,knownKinetics)
        flux(t,15) = params(67)*concMatrix(t,6)/(params(68)+concMatrix(t,6));
    end
    
    %Phosphoglycerate kinase
    if ismember(16,knownKinetics)
        flux(t,16) = params(69)*(params(2)*concMatrix(t,7)-params(1)*concMatrix(t,8)/params(70))/((params(71)*(1+params(1)/params(72))+params(2))*(params(73)*(1+concMatrix(t,8)/params(74))+concMatrix(t,7)));
    end
    
    %Serine synthesis
    if ismember(17,knownKinetics)
        flux(t,17) = params(75)*concMatrix(t,8)/(params(76)+concMatrix(t,8));
    end
    
    %Phosphoglycerate Mutase
    if ismember(18,knownKinetics)
        flux(t,18) = params(77)*(concMatrix(t,8)-concMatrix(t,9)/params(78))/(params(79)*(1+concMatrix(t,9)/params(80))+concMatrix(t,8));
    end
    
    %Enolase
    if ismember(19,knownKinetics)
        flux(t,19) = params(81)*(concMatrix(t,9)-concMatrix(t,10)/params(82))/(params(83)*(1+concMatrix(t,10)/params(84))+concMatrix(t,9));
    end
    
    %Pyruvate kinase
    if ismember(20,knownKinetics)
        flux(t,20) = params(85)*concMatrix(t,10)*(concMatrix(t,10)/params(86)+1)^(params(87)-1)*params(2)/(params(86)*(params(88)*((1+params(1)/params(89))/(concMatrix(t,4)/params(90)+params(3)/params(91)+1))^params(87)+(concMatrix(t,10)/params(86)+1)^params(87))*(params(2)+params(92)));
    end
    
    %PEP carboxylase
    if ismember(21,knownKinetics)
        flux(t,21) = params(93)*concMatrix(t,10)*(1+(concMatrix(t,4)/params(94))^params(95))/(params(96)+concMatrix(t,10));
    end
    
    %Synthesis 1
    if ismember(22,knownKinetics)
        flux(t,22) = params(97)*concMatrix(t,10)/(params(98)+concMatrix(t,10));
    end
    
    %Synthesis 2
    if ismember(23,knownKinetics)
        flux(t,23) = params(99)*concMatrix(t,11)/(params(100)+concMatrix(t,11));
    end
    
    %DAHP synthesis
    if ismember(24,knownKinetics)
        flux(t,24) = params(101)*concMatrix(t,17)^params(102)*concMatrix(t,10)^params(103)/((params(104)+concMatrix(t,17)^params(102))*(params(105)+concMatrix(t,10)^params(103)));
    end
    
    %Pyruvate dehydrogenase
    if ismember(25,knownKinetics)
        flux(t,25) = params(106)*concMatrix(t,11)^params(107)/(params(108)+concMatrix(t,11)^params(107));
    end
    
    %Methionine synthesis
    if ismember(26,knownKinetics)
        flux(t,26) = params(109);
    end
    
    %6-Phosphogluconate dehydrogenase
    if ismember(27,knownKinetics)
        flux(t,27) = params(110)*concMatrix(t,12)*params(4)/((concMatrix(t,12)+params(111))*(params(4)+params(112)*(1+params(5)/params(113))*(1+params(1)/params(114))));
    end
    
    %Ribose-phosphate isomerase
    if ismember(28,knownKinetics)
        flux(t,28) = params(115)*(concMatrix(t,13)-concMatrix(t,16)/params(116));
    end
    
    %Ribulose-phosphate epimerase
    if ismember(29,knownKinetics)
        flux(t,29) = params(117)*(concMatrix(t,13)-concMatrix(t,14)/params(118));
    end
    
    %Ribose phosphate pyrophosphokinase
    if ismember(30,knownKinetics)
        flux(t,30) = params(119)*concMatrix(t,16)/(params(120)+concMatrix(t,16));
    end
    
    %Glucose-1-phosphate adenyltransferase
    if ismember(31,knownKinetics)
        flux(t,31) = params(121)*concMatrix(t,18)*params(1)*(1+(concMatrix(t,4)/params(122))^params(123))/((params(124)+params(1))*(params(125)+concMatrix(t,18)));
    end
    
    %G6P degradation
    if ismember(32,knownKinetics)
        flux(t,32) = params(126)*concMatrix(t,2);
    end
    
    %F6P degradation
    if ismember(33,knownKinetics)
        flux(t,33) = params(126)*concMatrix(t,3);
    end
    
    %FDP degradation
    if ismember(34,knownKinetics)
        flux(t,34) = params(126)*concMatrix(t,4);
    end
    
    %GAP degradation
    if ismember(35,knownKinetics)
        flux(t,35) = params(126)*concMatrix(t,5);
    end
    
    %DHAP degradation
    if ismember(36,knownKinetics)
        flux(t,36) = params(126)*concMatrix(t,6);
    end
    
    %PGP degradation
    if ismember(37,knownKinetics)
        flux(t,37) = params(126)*concMatrix(t,7);
    end
    
    %PG3 degradation
    if ismember(38,knownKinetics)
        flux(t,38) = params(126)*concMatrix(t,8);
    end
    
    %PG2 degradation
    if ismember(39,knownKinetics)
        flux(t,39) = params(126)*concMatrix(t,9);
    end
    
    %PEP degradation
    if ismember(40,knownKinetics)
        flux(t,40) = params(126)*concMatrix(t,10);
    end
    
    %Pyruvate dilution
    if ismember(41,knownKinetics)
        flux(t,41) = params(126)*concMatrix(t,11);
    end
    
    %PG dilution
    if ismember(42,knownKinetics)
        flux(t,42) = params(126)*concMatrix(t,12);
    end
    
    %Ribu5P dilution
    if ismember(43,knownKinetics)
        flux(t,43) = params(126)*concMatrix(t,13);
    end
    
    %XYL5P dilution
    if ismember(44,knownKinetics)
        flux(t,44) = params(126)*concMatrix(t,14);
    end
    
    %SED7P dilution
    if ismember(45,knownKinetics)
        flux(t,45) = params(126)*concMatrix(t,15);
    end
    
    %Rib5P dilution
    if ismember(46,knownKinetics)
        flux(t,46) = params(126)*concMatrix(t,16);
    end
    
    %E4P dilution
    if ismember(47,knownKinetics)
        flux(t,47) = params(126)*concMatrix(t,17);
    end
    
    %GLP dilution
    if ismember(48,knownKinetics)
        flux(t,48) = params(126)*concMatrix(t,18);
    end
    
    %Pooling fluxes
    flux(t,49) = (concMatrix(t+1,1) - concMatrix(t,1))/deltaT;
    flux(t,50) = (concMatrix(t+1,2) - concMatrix(t,2))/deltaT;
    flux(t,51) = (concMatrix(t+1,3) - concMatrix(t,3))/deltaT;
    flux(t,52) = (concMatrix(t+1,4) - concMatrix(t,4))/deltaT;
    flux(t,53) = (concMatrix(t+1,5) - concMatrix(t,5))/deltaT;
    flux(t,54) = (concMatrix(t+1,6) - concMatrix(t,6))/deltaT;
    flux(t,55) = (concMatrix(t+1,7) - concMatrix(t,7))/deltaT;
    flux(t,56) = (concMatrix(t+1,8) - concMatrix(t,8))/deltaT;
    flux(t,57) = (concMatrix(t+1,9) - concMatrix(t,9))/deltaT;
    flux(t,58) = (concMatrix(t+1,10) - concMatrix(t,10))/deltaT;
    flux(t,59) = (concMatrix(t+1,11) - concMatrix(t,11))/deltaT;
    flux(t,60) = (concMatrix(t+1,12) - concMatrix(t,12))/deltaT;
    flux(t,61) = (concMatrix(t+1,13) - concMatrix(t,13))/deltaT;
    flux(t,62) = (concMatrix(t+1,14) - concMatrix(t,14))/deltaT;
    flux(t,63) = (concMatrix(t+1,15) - concMatrix(t,15))/deltaT;
    flux(t,64) = (concMatrix(t+1,16) - concMatrix(t,16))/deltaT;
    flux(t,65) = (concMatrix(t+1,17) - concMatrix(t,17))/deltaT;
    flux(t,66) = (concMatrix(t+1,18) - concMatrix(t,18))/deltaT;
end
