%This script serves to determine the penalty weight specific to each system
%studied. The idea is to balance the magnitude of penalties to ensure that
%mass balance penalty is dominant and that no other one penalty has a large
%numerical value that it overthrows other penalties 

system_name = '';%{Determined,UDreg,ecoli,yeast}

%Load sampling results
results = load(sprintf('%s_sampling_final.mat',system_name));
penaltyTable = results.penaltyTable; 
penaltyTable_gs = results.penaltyTable_gs;

% use the minimum value in each penalty as the reference point and
% bring those to the same magnitude
MBpenalty_values = penaltyTable.sortedMBPen; 
conc_values = penaltyTable.sortedConcPen;
oneCorr_values = penaltyTable.sortedOneCorrPen;
oneCurve_values = penaltyTable.sortedOneCurvePen;
BST_values = penaltyTable.sortedBstPen;
if ismember('sortedssPen',penaltyTable.Properties.VariableNames)
    ss_values = penaltyTable.sortedssPen;
end

all_penalties = [MBpenalty_values,conc_values,oneCorr_values,oneCurve_values,BST_values,ss_values];

min_all_penalties = min(all_penalties);
reference_penalty_weight = zeros(size(min_all_penalties)) - round(log10(min_all_penalties));
reference_penalty_weight(reference_penalty_weight>1) = 1; %We do not wish to over-penalize any one penalty term 

%Obtain the penalty weight vector for the system. The weight for mass
%balance is multiplied by 100. 
penalty_weight_vec = 10.^(reference_penalty_weight);
penalty_weight_vec(1) = penalty_weight_vec(1) * 100;


