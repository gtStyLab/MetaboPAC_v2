function flux = calcFluxWithKinetics_hynne(concMatrix,timeVec,knownKinetics)
% Calculate flux reaction rates in yeast system using known kinetic
% equations.
% Parameters and kinetics were extracted from Full-scale model of
% glycolysis in Saccharomyces cerevisiae, Hynne and Sorensen 2001.

nT = size(concMatrix,1)-1;
deltaT = timeVec(2)-timeVec(1);

params(1) =  1014.96; % V2f
params(2) =  1.7; % K2Glc
params(3) =  1; % P2
params(4) =  1.2; % K2IG6P
params(5) =  7.2; % K2IIG6P
params(6) =  1014.96; % V2r
params(7) =  51.7547; % V3m
params(8) =  0.37; % K3DGlc
params(9) =  0.1; % K3ATP
params(10) =  0; % K3Glc
params(11) =  496.042; % V4f
params(12) =  0.8; % K4G6P
params(13) =  0.15; % K4F6P
params(14) =  496.042; % V4r
params(15) =  0.13; % K4eq
params(16) =  45.4327; % V5m
params(17) =  0.021; % K5
params(18) =  0.15; % kappa5
params(19) =  2207.82; % V6f
params(20) =  0.3; %K6FBP
params(21) =  2; % K6DHAP
params(22) =  0.081; % K6eq
params(23) =  5; % ratio6
params(24) =  4; % K6GAP
params(25) =  10; % K6IGAP
params(26) =  116.365; % V7f
params(27) =  1.23; % K7DHAP
params(28) =  1.27; % K7GAP
params(29) =  116.365; % V7r
params(30) =  0.055; % K7eq
params(31) =  833.858; % V8f
params(32) =  0.6; % K8GAP
params(33) =  0.1; % K8NAD
params(34) =  0.01; % K8BPG
params(35) =  0.06; % K8NADH
params(36) =  833.858; % V8r
params(37) =  0.0055; % K8eq
params(38) =  443866; % PEPsynth_kf
params(39) =  1528.62; % PEPsynth_kr
params(40) =  343.096; % V10m
params(41) =  0.2; % K10PEP
params(42) =  0.17; % K10ADP
params(43) =  53.1328; % V11m
params(44) =  0.3; % K11
params(45) =  89.8023; % V12m
params(46) =  0.1; % K12NADH
params(47) =  0.71; % K12ACA
params(48) =  16.72; % k13
params(49) =  81.4797; % V15m
params(50) =  25; % K15DHAP
params(51) =  0.034; % K15INADH
params(52) =  0.13; % K15INAD
params(53) =  0.13; % K15NADH
params(54) =  1.9; % k16
params(55) =  24.7; % k18
params(56) =  0.00283828; % k20
params(57) =  2.25932; % k22
params(58) =  3.2076; % k23
params(59) =  432.9; % Adenylate_kinase_kf
params(60) =  133.333; % Adenylate_kinase_kr

params(61) =  59; % Yvol
params(62) =  0.048; % k0

params(63) =  8.5; % GlcX0 (default 18.5)
params(64) =  19; % CNX0


flux = nan(nT,46);
for t = 1:nT% List of reactions
    %Glucose Mixed flow to extracellular medium
    if ismember(1,knownKinetics)
        flux(t,1) = params(61)*params(62)*(params(63)-concMatrix(t,17));
    end

    %Glucose uptake
    if ismember(2,knownKinetics)
    flux(t,2) = (params(1))*(concMatrix(t,17)/params(2))/(1 + concMatrix(t,17)/params(2) + ((params(3)*concMatrix(t,17)/params(2) + 1)/(params(3)*concMatrix(t,2)/params(2) + 1)) * (1 + concMatrix(t,2)/params(2) + concMatrix(t,8)/params(4) + (concMatrix(t,2)*concMatrix(t,8))/(params(2)*params(5))))...
        - (params(6))*(concMatrix(t,2)/params(2))/(1 + concMatrix(t,2)/params(2) + ((params(3)*concMatrix(t,2)/params(2) + 1)/(params(3)*concMatrix(t,17)/params(2) + 1)) * (1 + concMatrix(t,2)/params(2) + concMatrix(t,8)/params(4) + (concMatrix(t,2)*concMatrix(t,8))/(params(2)*params(5))));
    end
    
    %Hexokinase
    if ismember(3,knownKinetics)
        flux(t,3) = (params(7)*concMatrix(t,21)*concMatrix(t,2))/(params(8)*params(9) + params(10)*concMatrix(t,21) + params(9)*concMatrix(t,2) + concMatrix(t,2)*concMatrix(t,21));
    end
    
    %Phosphoglucoisomerase
    if ismember(4,knownKinetics)
        flux(t,4) = (params(11)*concMatrix(t,8))/(params(12) + concMatrix(t,8) + (params(12)/params(13))*concMatrix(t,12))...
            - (params(14)*concMatrix(t,12)/params(15))/(params(12) + concMatrix(t,8) + (params(12)/params(13))*concMatrix(t,12));
    end
    
    %Phosphofructokinase
    if ismember(5,knownKinetics)
        flux(t,5) = (params(16)*concMatrix(t,12)^2)/(params(17)*(1 + params(18)*(concMatrix(t,21)/concMatrix(t,20))*(concMatrix(t,21)/concMatrix(t,20))) + concMatrix(t,12)^2);
    end
    
    %Aldolase
    if ismember(6,knownKinetics)
        flux(t,6) = (params(19)*concMatrix(t,19))/(params(20) + concMatrix(t,19) + (concMatrix(t,6)*params(21)*params(19))/(params(22)*params(19)*params(23)) + (concMatrix(t,7)*params(24)*params(19))/(params(22)*params(19)*params(23)) + concMatrix(t,19)*concMatrix(t,6)/params(25) + (concMatrix(t,6)*concMatrix(t,7)*params(19))/(params(22)*params(19)*params(23)))...
            - ((params(19)*concMatrix(t,6)*concMatrix(t,7))/params(22))/(params(20) + concMatrix(t,19) + (concMatrix(t,6)*params(21)*params(19))/(params(22)*params(19)*params(23)) + (concMatrix(t,7)*params(24)*params(19))/(params(22)*params(19)*params(23)) + concMatrix(t,19)*concMatrix(t,6)/params(25) + (concMatrix(t,6)*concMatrix(t,7)*params(19))/(params(22)*params(19)*params(23)));
    end
    
    %Triosephosphate isomerase
    if ismember(7,knownKinetics)
        flux(t,7) = (params(26)*concMatrix(t,7))/(params(27) + concMatrix(t,7) + (params(27)/params(28))*concMatrix(t,6)) - (params(29)*concMatrix(t,6)/params(30))/(params(27) + concMatrix(t,7) + (params(27)/params(28))*concMatrix(t,6));
    end
    
    %Glyceraldehyde 3-phosphate dehydrogenase
    if ismember(8,knownKinetics)
        flux(t,8) = ((params(31)*concMatrix(t,6)*concMatrix(t,9))/params(32)/params(33))/((1 + concMatrix(t,6)/params(32) + concMatrix(t,18)/params(34))*(1 + concMatrix(t,9)/params(33) + concMatrix(t,22)/params(35)))...
            - ((params(36)*concMatrix(t,18)*concMatrix(t,22))/params(37)/params(32)/params(33))/((1 + concMatrix(t,6)/params(32) + concMatrix(t,18)/params(34))*(1 + concMatrix(t,9)/params(33) + concMatrix(t,22)/params(35)));
    end
    
    %Phosphoenolpyruvate synthesis***
    if ismember(9,knownKinetics)
        flux(t,9) = params(38)*concMatrix(t,18)*concMatrix(t,5) - params(39)*concMatrix(t,11)*concMatrix(t,21);
    end
    
    %Pyruvate kinase
    if ismember(10,knownKinetics)
        flux(t,10) = (params(40)*concMatrix(t,5)*concMatrix(t,11))/((params(41) + concMatrix(t,11))*(params(42) + concMatrix(t,5)));
    end
    
    %Pyruvate decarboxylase
    if ismember(11,knownKinetics)
        flux(t,11) = (params(43)*concMatrix(t,14))/(params(44) + concMatrix(t,14));
    end
    
    %Alcohol dehydrogenase
    if ismember(12,knownKinetics)
        flux(t,12) = (params(45)*concMatrix(t,1)*concMatrix(t,22))/((params(46) + concMatrix(t,22))*(params(47) + concMatrix(t,1)));
    end
    
    %Ethanol out
    if ismember(13,knownKinetics)
        flux(t,13) = (params(48))*(concMatrix(t,4) - concMatrix(t,15));
    end
    
    %Ethanol flow
    if ismember(14,knownKinetics)
        flux(t,14) = params(61) * params(62)*concMatrix(t,15);
    end
    
    %Glycerol synthesis
    if ismember(15,knownKinetics)
        flux(t,15) = (params(49)*concMatrix(t,7))/(params(50)*(1 + (params(51)/concMatrix(t,22))*(1 + concMatrix(t,9)/params(52))) + concMatrix(t,7)*(1 + (params(53)/concMatrix(t,22))*(1 + concMatrix(t,9)/params(52))));
    end
    
    %Glycerol out
    if ismember(16,knownKinetics)
        flux(t,16) = (params(54))*(concMatrix(t,3) - concMatrix(t,16));
    end
    
    %Glycerol flow
    if ismember(17,knownKinetics)
        flux(t,17) = params(61) * params(62)*concMatrix(t,16);
    end
    
    %Acetaldehyde out
    if ismember(18,knownKinetics)
        flux(t,18) = (params(55))*(concMatrix(t,1) - concMatrix(t,10));
    end
    
    %Acetaldehyde flow
    if ismember(19,knownKinetics)
        flux(t,19) = params(61) * params(62)*concMatrix(t,10);
    end
    
    %Cyanide-Acetaldehyde flow
    if ismember(20,knownKinetics)
        flux(t,20) = params(61) * params(56)*concMatrix(t,10)*concMatrix(t,13);
    end
    
    %Cyanide flow
    if ismember(21,knownKinetics)
        flux(t,21) = params(61) * params(62)*(params(64) - concMatrix(t,13));
    end
    
    %Storage
    if ismember(22,knownKinetics)
        flux(t,22) = params(57)*concMatrix(t,21)*concMatrix(t,8);
    end
    
    %ATP consumption
    if ismember(23,knownKinetics)
        flux(t,23) = params(58)*concMatrix(t,21);
    end
    
    %Adenylate kinase***
    if ismember(24,knownKinetics)
        flux(t,24) = params(59)*concMatrix(t,21)*concMatrix(t,20) - params(60)*concMatrix(t,5)*concMatrix(t,5);
    end
    
    %Pooling fluxes
    flux(t,25) = (concMatrix(t+1,1) - concMatrix(t,1))/deltaT;
    flux(t,26) = (concMatrix(t+1,2) - concMatrix(t,2))/deltaT;
    flux(t,27) = (concMatrix(t+1,3) - concMatrix(t,3))/deltaT;
    flux(t,28) = (concMatrix(t+1,4) - concMatrix(t,4))/deltaT;
    flux(t,29) = (concMatrix(t+1,5) - concMatrix(t,5))/deltaT;
    flux(t,30) = (concMatrix(t+1,6) - concMatrix(t,6))/deltaT;
    flux(t,31) = (concMatrix(t+1,7) - concMatrix(t,7))/deltaT;
    flux(t,32) = (concMatrix(t+1,8) - concMatrix(t,8))/deltaT;
    flux(t,33) = (concMatrix(t+1,9) - concMatrix(t,9))/deltaT;
    flux(t,34) = (concMatrix(t+1,10) - concMatrix(t,10))/deltaT;
    flux(t,35) = (concMatrix(t+1,11) - concMatrix(t,11))/deltaT;
    flux(t,36) = (concMatrix(t+1,12) - concMatrix(t,12))/deltaT;
    flux(t,37) = (concMatrix(t+1,13) - concMatrix(t,13))/deltaT;
    flux(t,38) = (concMatrix(t+1,14) - concMatrix(t,14))/deltaT;
    flux(t,39) = (concMatrix(t+1,15) - concMatrix(t,15))/deltaT;
    flux(t,40) = (concMatrix(t+1,16) - concMatrix(t,16))/deltaT;
    flux(t,41) = (concMatrix(t+1,17) - concMatrix(t,17))/deltaT;
    flux(t,42) = (concMatrix(t+1,18) - concMatrix(t,18))/deltaT;
    flux(t,43) = (concMatrix(t+1,19) - concMatrix(t,19))/deltaT;
    flux(t,44) = (concMatrix(t+1,20) - concMatrix(t,20))/deltaT;
    flux(t,45) = (concMatrix(t+1,21) - concMatrix(t,21))/deltaT;
    flux(t,46) = (concMatrix(t+1,22) - concMatrix(t,22))/deltaT;
end
