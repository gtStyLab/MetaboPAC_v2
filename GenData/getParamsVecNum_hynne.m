function vecOut = getParamsVecNum_hynne(n)
% Just a hard-coded list of parameter vectors, by index

    % Vary ICs
    paramsList{1}  = [
        1014.96; % V2f
        1.7; % K2Glc
        1; % P2
        1.2; % K2IG6P
        7.2; % K2IIG6P
        1014.96; % V2r
        51.7547; % V3m
        0.37; % K3DGlc
        0.1; % K3ATP
        0; % K3Glc
        496.042; % V4f
        0.8; % K4G6P
        0.15; % K4F6P
        496.042; % V4r
        0.13; % K4eq
        45.4327; % V5m
        0.021; % K5
        0.15; % kappa5
        2207.82; % V6f
        0.3; %K6FBP
        2; % K6DHAP
        0.081; % K6eq
        5; % ratio6
        4; % K6GAP
        10; % K6IGAP
        116.365; % V7f
        1.23; % K7DHAP
        1.27; % K7GAP
        116.365; % V7r
        0.055; % K7eq
        833.858; % V8f
        0.6; % K8GAP
        0.1; % K8NAD
        0.01; % K8BPG
        0.06; % K8NADH
        833.858; % V8r
        0.0055; % K8eq
        443866; % PEPsynth_kf
        1528.62; % PEPsynth_kr
        343.096; % V10m
        0.2; % K10PEP
        0.17; % K10ADP
        53.1328; % V11m
        0.3; % K11
        89.8023; % V12m
        0.1; % K12NADH
        0.71; % K12ACA
        16.72; % k13
        81.4797; % V15m
        25; % K15DHAP
        0.034; % K15INADH
        0.13; % K15INAD
        0.13; % K15NADH
        1.9; % k16
        24.7; % k18
        0.00283828; % k20
        2.25932; % k22
        3.2076; % k23
        432.9; % Adenylate_kinase_kf
        133.333; % Adenylate_kinase_kr
        
        59; % Yvol
        0.048; % k0
        
        8.5; % GlcX0 (default 18.5)
        19; % CNX0
        
        1.48153; % x0_1 Acetaldehyde
        0.573074; % x0_2 Cytsosolic glucose
        4.196; % x0_3 Glycerol
        19.2379; % x0_4 EtOH
        1.5; % x0_5 ADP
        0.115; % x0_6 Glyceraldehyde 3-phosphate
        2.95; % x0_7 Dihydroxyacetone phosphate
        4.2; % x0_8 Glucose-6-Phosphate
        0.65; % x0_9 NAD
        1.28836; % x0_10 Extracellular acetaldehyde
        0.04; % x0_11 Phosphoenolpyruvate
        0.49; % x0_12 Fructose-6-Phosphate
        5.20358; % x0_13 Extracellular cyanide
        8.7; % x0_14 Pyruvate
        16.4514; % x0_15 Extracellular ethanol
        1.68478; % x0_16 Extracellular glycerol
        1.55307; % x0_17 Extracellular Glucose (default 1.55307)
        0.00027; % x0_18 1,3-Bisphosphoglycerate
        4.64; % x0_19 Fructose 1,6-bisphosphate
        0.33; % x0_20 AMP
        2.1; % x0_21 ATP
        0.33;]; % x0_22 NADH
    
    vecOut = paramsList{n};

end


