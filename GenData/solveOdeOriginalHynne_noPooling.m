function [timeVec, concMatrix, fluxMatrix] = solveOdeOriginalHynne_noPooling(tStart,tEnd,nT,x0,params,stm)
% This function is hard-coded to integrate the branched pathway model as an
% ODE with Michaelis-Menten kinetic rate laws described in the
% Supplementary Methods. 
%
% Written by R.A. Dromms 2015-07-29
% Edited for Hynne S. cerevisiae model JL17May20
    %{
    if exist('params','var')
        params = convertOdeParams(params);
    else
        params = setOdeParams;
    end
    %}
%     [timeVec,concMatrix] = ode45(@(t,x)fRHS(t,x,params,stm),linspace(tStart,tEnd,nT+1),x0');
    
    % From solveOdeBst...
    [timeVec,concMatrix] = ode15s(@(t,x)fRHS(t,x,params,stm),linspace(tStart,tEnd,nT+1),x0,odeset('NonNegative',1:length(x0)));
    
    fluxMatrix = zeros(length(timeVec),24);
    for k = 1:length(timeVec)
        fluxMatrix(k,1:24) = calcFluxes(timeVec(k),concMatrix(k,:)',params)';
        %fluxMatrix(k,25:42) = fRHS(timeVec(k),concMatrix(k,:)',params,stm)';
    end
    
end

function xdot = fRHS(~,x,params,stm)
    
    x(x<0) = 0;
    v = calcFluxes([],x,params);
    
    xdot = stm*v; 

end

function v = calcFluxes(~,x,params)
    
%     % Bad Things can happen with the math when we allow negative x
%     % (Which isn't physically relevant, anyways)
%     xMin = 1e-4;
%     x(x<xMin) = xMin;
    
    % List of reactions
    %Glucose Mixed flow to extracellular medium
    v(1,1) = params(61) * params(62)*(params(63)-x(17));

    %Glucose uptake
    v(2,1) = (params(1))*(x(17)/params(2))/(1 + x(17)/params(2) + ((params(3)*x(17)/params(2) + 1)/(params(3)*x(2)/params(2) + 1)) * (1 + x(2)/params(2) + x(8)/params(4) + (x(2)*x(8))/(params(2)*params(5))))...
        - (params(6))*(x(2)/params(2))/(1 + x(2)/params(2) + ((params(3)*x(2)/params(2) + 1)/(params(3)*x(17)/params(2) + 1)) * (1 + x(2)/params(2) + x(8)/params(4) + (x(2)*x(8))/(params(2)*params(5))));
    
    %Hexokinase
    v(3,1) = (params(7)*x(21)*x(2))/(params(8)*params(9) + params(10)*x(21) + params(9)*x(2) + x(2)*x(21));

    %Phosphoglucoisomerase
    v(4,1) = (params(11)*x(8))/(params(12) + x(8) + (params(12)/params(13))*x(12))...
        - (params(14)*x(12)/params(15))/(params(12) + x(8) + (params(12)/params(13))*x(12));
    
    %Phosphofructokinase
    v(5,1) = (params(16)*x(12)^2)/(params(17)*(1 + params(18)*(x(21)/x(20))*(x(21)/x(20))) + x(12)^2);
    
    %Aldolase
    v(6,1) = (params(19)*x(19))/(params(20) + x(19) + (x(6)*params(21)*params(19))/(params(22)*params(19)*params(23)) + (x(7)*params(24)*params(19))/(params(22)*params(19)*params(23)) + x(19)*x(6)/params(25) + (x(6)*x(7)*params(19))/(params(22)*params(19)*params(23)))...
        - ((params(19)*x(6)*x(7))/params(22))/(params(20) + x(19) + (x(6)*params(21)*params(19))/(params(22)*params(19)*params(23)) + (x(7)*params(24)*params(19))/(params(22)*params(19)*params(23)) + x(19)*x(6)/params(25) + (x(6)*x(7)*params(19))/(params(22)*params(19)*params(23)));

    %Triosephosphate isomerase
    v(7,1) = (params(26)*x(7))/(params(27) + x(7) + (params(27)/params(28))*x(6)) - (params(29)*x(6)/params(30))/(params(27) + x(7) + (params(27)/params(28))*x(6));

    %Glyceraldehyde 3-phosphate dehydrogenase
    v(8,1) = ((params(31)*x(6)*x(9))/params(32)/params(33))/((1 + x(6)/params(32) + x(18)/params(34))*(1 + x(9)/params(33) + x(22)/params(35)))...
        - ((params(36)*x(18)*x(22))/params(37)/params(32)/params(33))/((1 + x(6)/params(32) + x(18)/params(34))*(1 + x(9)/params(33) + x(22)/params(35)));
    
    %Phosphoenolpyruvate synthesis***
    v(9,1) = params(38)*x(18)*x(5) - params(39)*x(11)*x(21);
    
    %Pyruvate kinase
    v(10,1) = (params(40)*x(5)*x(11))/((params(41) + x(11))*(params(42) + x(5)));

    %Pyruvate decarboxylase
    v(11,1) = (params(43)*x(14))/(params(44) + x(14));

    %Alcohol dehydrogenase
    v(12,1) = (params(45)*x(1)*x(22))/((params(46) + x(22))*(params(47) + x(1)));
    
    %Ethanol out
    v(13,1) = params(48)*(x(4) - x(15));

    %Ethanol flow
    v(14,1) = params(61) * params(62)*x(15);

    %Glycerol synthesis
    v(15,1) = (params(49)*x(7))/(params(50)*(1 + (params(51)/x(22))*(1 + x(9)/params(52))) + x(7)*(1 + (params(53)/x(22))*(1 + x(9)/params(52))));

    %Glycerol out
    v(16,1) = params(54) *(x(3) - x(16));

    %Glycerol flow
    v(17,1) = params(61) * params(62)*x(16);

    %Acetaldehyde out
    v(18,1) = params(55)*(x(1) - x(10));

    %Acetaldehyde flow
    v(19,1) = params(61) * params(62)*x(10);

    %Cyanide-Acetaldehyde flow
    v(20,1) = params(61) * params(56)*x(10)*x(13);

    %Cyanide flow
    v(21,1) = params(61) * params(62)*(params(64) - x(13));

    %Storage
    v(22,1) = params(57)*x(21)*x(8);

    %ATP consumption
    v(23,1) = params(58)*x(21);

    %Adenylate kinase***
    v(24,1) = params(59)*x(21)*x(20) - params(60)*x(5)*x(5);
    
end
%{
    % List of reactions
    %Glucose Mixed flow to extracellular medium
    v(1,1) = k0*(GlcX0-GlcX);

    %Glucose uptake
    v(2,1) = (V2f/Yvol)*([GlcX]/K2Glc)/(1 + [GlcX]/K2Glc + ((P2*[GlcX]/K2Glc + 1)/(P2*[Glc]/K2Glc + 1)) * (1 + [Glc]/K2Glc + [G6P]/K2IG6P + ([Glc]*[G6P])/(K2Glc*K2IIG6P)))...
        - (V2r/Yvol)*([Glc]/K2Glc)/(1 + [Glc]/K2Glc + ((P2*[Glc]/K2Glc + 1)/(P2*[GlcX]/K2Glc + 1)) * (1 + [Glc]/K2Glc) + [G6P]/K2IG6P + ([Glc]*[G6P])/(K2Glc*K2IIG6P));
    
    %Hexokinase
    v(3,1) = (V3m*[ATP]*[Glc])/(K3DGlc*K3ATP + K3Glc*[ATP] + K3ATP*[Glc] + [Glc]*[ATP]);

    %Phosphoglucoisomerase
    v(4,1) = (V4f*[G6P])/(K4G6P + [G6P] + (K4G6P/K4F6P)*[F6P])...
        - (V4r*[F6P]/K4eq)/(K4G6P + [G6P] + (K4G6P/K4F6P)*[F6P]);
    
    %Phosphofructokinase
    v(5,1) = (V5m*[F6P]^2)/(K5*(1 + kappa5*([ATP]/[AMP])*([ATP]/[AMP])) + [F6P]^2);
    
    %Aldolase
    v(6,1) = (V6f*[FBP])/(K6FBP + [FBP] + ([GAP]*K6DHAP*V6f)/(K6eq*V6f*ratio6) + ([DHAP]*K6GAP*V6f)/(K6eq*V6f*ratio6) + [FBP]*[GAP]/K6IGAP + ([GAP]*[DHAP]*V6f)/(K6eq*V6f*ratio6))...
        - ((V6f*[GAP]*[DHAP])/K6eq)/(K6FBP + [FBP] + ([GAP]*K6DHAP*V6f)/(K6eq*V6f*ratio6) + ([DHAP]*K6GAP*V6f)/(K6eq*V6f*ratio6) + [FBP]*[GAP]/K6IGAP + ([GAP]*[DHAP]*V6f)/(K6eq*V6f*ratio6));
    
    %Triosephosphate isomerase
    v(7,1) = (V7f*[DHAP])/(K7DHAP + [DHAP] + (K7DHAP/K7GAP)*[GAP]) - (V7r*[GAP]/K7eq)/(K7DHAP + [DHAP] + (K7DHAP/K7GAP)*[GAP]);

    %Glyceraldehyde 3-phosphate dehydrogenase
    v(8,1) = ((V8f*[GAP]*[NAD])/K8GAP/K8NAD)/((1 + [GAP]/K8GAP + [BPG]/K8BPG)*(1 + [NAD]/K8NAD + [NADH]/K8NADH))...
        - ((V8r*[BPG]*[NADH])/K8eq/K8GAP/K8NAD)/((1 + [GAP]/K8GAP + [BPG]/K8BPG)*(1 + [NAD]/K8NAD + [NADH]/K8NADH));
    
    %Phosphoenolpyruvate synthesis***
    v(9,1) = PEPsynth_kf*[13bpg]*[ADP] - PEPsynth_kr*[PEP]*[ATP];
    
    %Pyruvate kinase
    v(10,1) = (V10m*[ADP]*[PEP])/((K10PEP + [PEP])*(K10ADP + [ADP]));

    %Pyruvate decarboxylase
    v(11,1) = (V11m*[Pyr])/(K11 + [Pyr]);

    %Alcohol dehydrogenase
    v(12,1) = (V12m*[ACA]*[NADH])/((K12NADH + [NADH])*(K12ACA + [ACA]));
    
    %Ethanol out
    v(13,1) = (k13/Yvol)*([EtOH] - [EtOHX]);

    %Ethanol flow
    v(14,1) = k0*[EtOHX];

    %Glycerol synthesis
    v(15,1) = (V15m*[DHAP])/(K15DHAP*(1 + (K15INADH/[NADH])*(1 + [NAD]/K15INAD)) + [DHAP]*(1 + (K15NADH/[NADH])*(1 + [NAD]/K15INAD)));

    %Glycerol out
    v(16,1) = (k16/Yvol)*([Glyc] - [GlycX]);

    %Glycerol flow
    v(17,1) = k0*[GlycX];

    %Acetaldehyde out
    v(18,1) = (k18/Yvol)*([ACA] - [ACAX]);

    %Acetaldehyde flow
    v(19,1) = k0*[ACAX];

    %Cyanide-Acetaldehyde flow
    v(20,1) = k20*[ACAX]*[CNX];

    %Cyanide flow
    v(21,1) = k0*([CNX0] - [CNX]);

    %Storage
    v(22,1) = k22*[ATP]*[G6P];

    %ATP consumption
    v(23,1) = k23*[ATP];

    %Adenylate kinase***
    v(24,1) = Adenylate_kinase_kf*[ATP]*[AMP] - Adenylate_kinase_kr*[ADP]*[ADP];
%}