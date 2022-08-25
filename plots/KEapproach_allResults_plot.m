%Plot results of all systems when 100% kinetics are known 
%Figure 2, Figure S2

%error threshold
error_threshold_list = [log2(1.1),log2(1.3),log2(1.5)];

%trueRF fileName format
ecoli_trueRF_file = load('trueRF/chass_trueRF.mat');
ecoli_trueRF = ecoli_trueRF_file.trueRF;
yeast_trueRF_file = load('trueRF/hynne_trueRF.mat');
yeast_trueRF = yeast_trueRF_file.trueRF;
UDreg_trueRF_file = load('trueRF/UDreg_trueRF.mat');
UDreg_trueRF = UDreg_trueRF_file.trueRF;
Determined_trueRF_file = load('trueRF/Determined_trueRF.mat');
Determined_trueRF = Determined_trueRF_file.trueRF;

%Load results files and analyze 
count_ecoli = nan(50,3);
count_yeast = nan(50,3);
count_UDreg = nan(50,3);
count_Determined = nan(50,3); 
%results fileName format
log2error_ecoli = nan(50,18);
log2error_yeast = nan(50,22);
log2error_UDreg = nan(50,4);
log2error_Determined = nan(50,4); 
for rep = 1:1:5
    for rand_idx = 1:1:10
        idx = (rep - 1) * 10 + rand_idx; 
        ecoli_fileName = sprintf('results/ecoli_MetaboPAC_KE_addKinetics-100_rep-%03d_rand-%03d',rep,rand_idx);
        yeast_fileName = sprintf('results/yeast_MetaboPAC_k-01_addKinetics-100_rep-%03d_rand-%03d',rep,rand_idx);
        ecoli_results_file = load(ecoli_fileName); 
        predicted_RF_ecoli = ecoli_results_file.predicted_responseFactors;
        yeast_results_file = load(yeast_fileName); 
        predicted_RF_yeast = yeast_results_file.predicted_responseFactors;
    
        log2error_ecoli(idx,:) = abs(log2(predicted_RF_ecoli) - log2(ecoli_trueRF(rep,:)));

        log2error_yeast(idx,:) = abs(log2(predicted_RF_yeast) - log2(yeast_trueRF(rep,:)));

        UDreg_fileName = sprintf('results/UDreg_MetaboPAC_KE_addKinetics-100_rep-%03d_rand-%03d.mat',rep,rand_idx);
        UDreg_results_file = load(UDreg_fileName);
        predicted_RF_UDreg = UDreg_results_file.predicted_responseFactors;
        log2error_UDreg(idx,:) = abs(log2(predicted_RF_UDreg) - log2(UDreg_trueRF(rep,:)));

        Determined_fileName = sprintf('results/Determined_MetaboPAC_KE_percKnownKinetics-100_rep-%03d_rand-%03d.mat',rep,rand_idx);
        Determined_results_file = load(Determined_fileName);
        predicted_RF_Determined = Determined_results_file.predicted_responseFactors;
        log2error_Determined(idx,:) = abs(log2(predicted_RF_Determined) - log2(Determined_trueRF(rep,:)));

        for error_threshold_idx = 1:1:length(error_threshold_list)
            error_threshold = error_threshold_list(error_threshold_idx);
            count_Determined(idx,error_threshold_idx) = sum(log2error_Determined(idx,:) < error_threshold);
            count_ecoli(idx,error_threshold_idx) = sum(log2error_ecoli(idx,:) < error_threshold);
            count_yeast(idx,error_threshold_idx) = sum(log2error_yeast(idx,:) < error_threshold);
            count_UDreg(idx,error_threshold_idx) = sum(log2error_UDreg(idx,:) < error_threshold);
        end
    end
end
percentage_ecoli = 100 .* count_ecoli./18;
percentage_UDreg = 100 .*count_UDreg./4;
percentage_Determined = 100 .* count_Determined/4;
percentage_yeast = 100 .* count_yeast./22;
avg_count = [mean(percentage_Determined);mean(percentage_UDreg,1);mean(percentage_ecoli);mean(percentage_yeast);]';
sd_count = [std(percentage_Determined,0,1);std(percentage_UDreg,0,1);std(percentage_ecoli,0,1);std(percentage_yeast,0,1)]';
standard_error = sd_count ./ avg_count; 



X = categorical({'<log2(1.1)','<log2(1.3)','<log2(1.5)'});
X = reordercats(X,{'<log2(1.1)','<log2(1.3)','<log2(1.5)'}); 
b = bar(X,avg_count,'grouped');
hold on 
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(avg_count);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',avg_count,standard_error,'k','linestyle','none');
hold off

legend('Determined','Underdetermined with Regulation','E. coli','S. cerevisiae','Location','bestoutside')
hold on 
ylim([0 120]) 
ylabel({'percentage of response factors','predicted within error range'})
% title({'100% known kinetics - KE approach'})
set(gca,'FontSize',14)












