%% Generate true response factors for each model

%% Determined model

clear;clc;

maxRandVal = 1000;
minRandVal = 1;
numMetabs = 4;
numSets = 100;

trueRF = (maxRandVal - minRandVal).*rand(numSets,numMetabs) + minRandVal;

save('Determined_trueRF','trueRF');

%% Determined model (log)

clear;clc;

trueRF_lower_mag = 0;
trueRF_upper_mag = 3;
numMetabs = 4;
numSets = 100;

trueRF = 10.^(trueRF_lower_mag + (trueRF_upper_mag-trueRF_lower_mag)*rand(numSets,numMetabs));

save('Determined_trueRF_log','trueRF');

%% Underdetermined model with regulation

clear;clc;

maxRandVal = 1000;
minRandVal = 1;
numMetabs = 4;
numSets = 100;

trueRF = (maxRandVal - minRandVal).*rand(numSets,numMetabs) + minRandVal;

save('UDreg_trueRF','trueRF');

%% Underdetermined model with regulation (log)

clear;clc;

trueRF_lower_mag = 0;
trueRF_upper_mag = 3;
numMetabs = 4;
numSets = 100;

trueRF = 10.^(trueRF_lower_mag + (trueRF_upper_mag-trueRF_lower_mag)*rand(numSets,numMetabs));

save('UDreg_trueRF_log','trueRF');

%% Hynne model

clear;clc;

maxRandVal = 1000;
minRandVal = 1;
numMetabs = 22;
numSets = 20;

trueRF = (maxRandVal - minRandVal).*rand(numSets,numMetabs) + minRandVal;

save('hynne_trueRF','trueRF');

%% Hynne model (log)

clear;clc;

trueRF_lower_mag = 0;
trueRF_upper_mag = 3;
numMetabs = 22;
numSets = 20;

trueRF = 10.^(trueRF_lower_mag + (trueRF_upper_mag-trueRF_lower_mag)*rand(numSets,numMetabs));

save('hynne_trueRF_log','trueRF');

%% Chass model

clear;clc;

maxRandVal = 1000;
minRandVal = 1;
numMetabs = 18;
numSets = 20;

trueRF = (maxRandVal - minRandVal).*rand(numSets,numMetabs) + minRandVal;

save('chass_trueRF','trueRF');

%% Chass model (log)

clear;clc;

trueRF_lower_mag = 0;
trueRF_upper_mag = 3;
numMetabs = 18;
numSets = 20;

trueRF = 10.^(trueRF_lower_mag + (trueRF_upper_mag-trueRF_lower_mag)*rand(numSets,numMetabs));

save('chass_trueRF_log','trueRF');