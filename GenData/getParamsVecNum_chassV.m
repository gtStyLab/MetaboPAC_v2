function vecOut = getParamsVecNum_chassV(n)
% Just a hard-coded list of parameter vectors, by index

    % Vary ICs
    paramsList{1}  = [
        4.27; % catp
        0.582; % cadp
        0.954783; % camp
        0.196759; % cnadp
        0.062; % cnadph
        1.4644; % cnad
        0.0934; % cnadh
        2.78e-5; % Dil
        110.96; % cfeed
        7829.78; % rmaxPTS
        3082.3; % KPTSa1
        0.01; % KPTSa2
        245.3; % KPTSa3
        3.66; % nPTSg6p
        2.15; % KPTSg6p
        650.988; % rmaxPGI
        0.1725; % KPGIeq
        2.9; % KPGIg6p
        0.266; % KPGIf6p
        0.2; % KPGIf6ppginh
        0.2; % KPGlg6ppginh
        0.839824; % rmaxPGM
        0.196; % KPGMeq
        1.038; % KPGMg6p
        0.0136; % KPGMg1p
        1.3802; % rmaxG6PDH
        14.4; % KG6PDHg6p
        6.43; % KG6PDHnadphg6pinh
        0.0246; % KG6PDHnadp
        0.01; % KG6PDHnadphnadpinh
        1840.58; % rmaxPFK
        0.123; % KPFKatps
        4.14; % KPFKadpc
        0.325; % KPFKf6ps
        3.26; % KPFKpep
        3.89; % KPFKadpb
        3.2; % KPFKampb
        128; % KPFKadpa
        19.1; % KPFKampa
        5.62907e6; % LPFK
        11.1; % nPFK
        10.8716; % rmaxTA
        1.05; % KTAeq
        9.47338; % rmaxTKa
        1.2; % KTKaeq
        86.5586; % rmaxTKb
        10; % KTKbeq
        0.00043711; % rmaxmurSynth
        17.4146; % rmaxALDO
        0.144; % kALDOeq
        1.75; % kALDOfdp
        0.088; % kALDOgap
        2; % VALDOblf
        0.088; % kALDOdhap
        0.6; % kALDOgapinh
        921.594; % rmaxGAPDH
        0.63; % KGAPDHeq
        0.683; % KGAPDHgap
        1.04e-5; % KGAPDHpgp
        0.252; % KGAPDHnad
        1.09; % KGAPDHnadh
        68.6747; % rmaxTIS
        1.39; % kTISeq
        2.8; % kTISdhap
        0.3; % kTISgap
        0.001037; % rmaxTrpSynth
        0.0116204; % rmaxG3PDH
        1; % KG3PDHdhap
        3021.77; % rmaxPGK
        1934.4; % KPGKeq
        0.185; % KPGKadp
        0.653; % KPGKatp
        0.0468; % KPGKpgp
        0.473; % KPGKpg3
        0.0257121; % rmaxSerSynth
        1; % KSerSynthpg3
        89.0497; % rmaxPGlumu
        0.188; % KPGlumueq
        0.2; % KPGlumupg3
        0.369; % KPGlumupg2
        330.448; % rmaxENO
        6.73; % KENOeq
        0.1; % KENOpg2
        0.135; % KENOpep
        0.0611315; % rmaxPK
        0.31; % KPKpep
        4; % nPK
        1000; % LPK
        22.5; % KPKatp
        0.19; % KPKfdp
        0.2; % KPKamp
        0.26; % KPKadp
        0.107021; % rmaxpepCxylase
        0.7; % KpepCxylasefdp
        4.21; % npepCxylasefdp
        4.07; % KpepCxylasepep
        0.019539; % rmaxSynth1
        1; % KSynth1pep
        0.0736186; % rmaxSynth2
        1; % KSynth2pyr
        0.107953; % rmaxDAHPS
        2.6; % nDAHPSe4p
        2.2; % nDAHPSpep
        0.035; % KDAHPSe4p
        0.0053; % KDAHPSpep
        6.05953; % rmaxPDH
        3.68; % nPDH
        1159; % KPDHpyr
        0.0022627; % rmaxMetSynth
        16.2324; % rmaxPGDH
        37.5; % KPGDHpg
        0.0506; % KPGDHnadp
        0.0138; % KPGDHnadphinh
        208; % KPGDHatpinh
        4.83841; % rmaxR5PI
        4; % KR5PIeq
        6.73903; % rmaxRu5P
        1.4; % KRu5Peq
        0.0129005; % rmaxRPPK
        0.1; % KRPPKrib5p
        0.00752546; % rmaxG1PAT
        0.119; % KG1PATfdp
        1.2; % nG1PATfdp
        4.42; % KG1PATatp
        3.2; % KG1PATg1p
        2.78e-5; % mu
        
        1.67; % x0_1 Extracellular Glucose
        3.481128719; % x0_2 Glucose-6-Phosphate
        0.6001920503; % x0_3 Fructose-6-Phosphate
        0.2870670232; % x0_4 Fructose-1,6-bisphosphate
        0.2243062372; % x0_5 Glyceraldehyde-3-Phosphate
        0.1717249964; % x0_6 Dihydroxyacetonephosphate
        0.008194516317; % x0_7 1,3-diphosphosphoglycerate
        2.135151754; % x0_8 3-Phosphoglycerate
        0.399768912; % x0_9 2-Phosphoglycerate
        2.675065316; % x0_10 Phosphoenol pyruvate
        2.67400793; % x0_11 Pyruvate
        0.8068521166; % x0_12 6-Phosphogluconate
        0.1119250357; % x0_13 Ribulose-5-phosphate
        0.1392324284; % x0_14 Xylulose-5-phosphate
        0.2732399134; % x0_15 Sedoheptulose-7-phosphate
        0.4014264584; % x0_16 Ribose-5-phosphate
        0.09988007306; % x0_17 Erythrose-4-phosphate
        0.6513824412;]; % x0_18 Glucose-1-Phosphate
    

    vecOut = paramsList{n};

end


