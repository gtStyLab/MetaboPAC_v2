% Plot results for noiseless synthetic systems.

clear;clc;
% approach_strings = {'','KE_'};
% Figure 3
approach_string = '';
% Figure 5 , Figure S3, and Figure S4
% approach_string = 'KE_';


% for string_idx = 1:1:length(approach_strings) 
%     approach_string = approach_strings{string_idx};
for log2_error = [1.1 1.3 1.5]
    percKnown_count = 1;
    for percKnownFlux = [0 20 40 60 80 100]
        run_count = 1;
        for run = 1:5
            for rand_idx = 1:10
                run_count = 10 * (run-1) + rand_idx;
                fileName = sprintf('results/Determined_MetaboPAC_%spercKnownKinetics-%03d_rep-%03d_rand-%03d.mat',approach_string,percKnownFlux,run,rand_idx);

                if isfile(fileName)
                    load(fileName,'RF_opt','trueRF','numMetabs','maxRandVal','minRandVal','knownMet');

                    med_F = median(RF_opt(:,1:numMetabs));
                    rand_F = (maxRandVal-minRandVal).*rand(1,numMetabs) + minRandVal;
                    %rand_F = 10.^(trueRF_lower_mag + (trueRF_upper_mag-trueRF_lower_mag)*rand(1,numMetabs));
                    F500 = 500*ones(1,numMetabs);
                    num_knownMet(run_count,percKnown_count) = length(knownMet);
                    
                    abs_log2_error{run_count,percKnown_count} = abs(log2(med_F) - log2(trueRF));
                    abs_log2_error_kinetic{run_count,percKnown_count} = abs(log2(med_F(:,knownMet)) - log2(trueRF(knownMet)));
                    abs_log2_error_optimization{run_count,percKnown_count} = abs(log2(med_F(:,setdiff(1:numMetabs,knownMet))) - log2(trueRF(setdiff(1:numMetabs,knownMet))));

                    abs_log2_error_random{run_count,percKnown_count} = abs(log2(rand_F) - log2(trueRF));

                    abs_log2_error_F500{run_count,percKnown_count} = abs(log2(F500) - log2(trueRF));

                end
            end
        end
        percKnown_count = percKnown_count + 1;
    end
end

% Log2 error plot
f = figure(1);
f.Position = [500 500 1600 800];
annotation('line', [0.05 0.95], [0.5, 0.5])
a = annotation('textbox', [0.05 0.97 0 0],'String',{'A.'},'EdgeColor','none');
a.FontSize = 24;
a.FontWeight = 'bold';
b = annotation('textbox', [0.05 0.5 0 0],'String',{'B.'},'EdgeColor','none');
b.FontSize = 24;
b.FontWeight = 'bold';
subplot_count = 1;
for log2_error = [1.1 1.3 1.5]
    
    percKnown_count = 1;
    mean_framework_results = [];
    mean_random_results = [];
    mean_F500_results = [];
    std_framework_results = [];
    std_random_results = [];
    std_F500_results = [];
    numRuns = nan(1,size(abs_log2_error,2));
    for percKnownFlux = [0 20 40 60 80 100]
        for run = 1:size(abs_log2_error,1)
            framework_results(run,percKnown_count) = mean(abs_log2_error{run,percKnown_count} < log2(log2_error));
            random_results(run,percKnown_count) = mean(abs_log2_error_random{run,percKnown_count} < log2(log2_error));
            F500_results(run,percKnown_count) = mean(abs_log2_error_F500{run,percKnown_count} < log2(log2_error));
        end
        % manage the dummy 10 replicates of different subsets of known
        % kinetic equations for 0% and 100% known kinetics 
        % remove those and fill in nan values
        if isequal(percKnownFlux,0) || isequal(percKnownFlux,100)
            framework_results_v2 = nan(size(framework_results(:,percKnown_count)));
            framework_results_v2(1:10:end) = framework_results(1:10:end,percKnown_count);
            
            random_results_v2 = nan(size(random_results(:,percKnown_count)));
            random_results_v2(1:10:end) = random_results(1:10:end,percKnown_count);
            
            F500_results_v2 = nan(size(F500_results(:,percKnown_count)));
            F500_results_v2(1:10:end) = F500_results(1:10:end,percKnown_count);
        else
            framework_results_v2 = framework_results(:,percKnown_count);
            random_results_v2 = random_results(:,percKnown_count);
            F500_results_v2 = F500_results(:,percKnown_count);
        end
        
        numRuns(percKnown_count) = sum(~isnan(framework_results_v2));
        mean_framework_results(1,percKnown_count) = nanmean(framework_results_v2);
        mean_random_results(1,percKnown_count) = nanmean(random_results_v2);
        mean_F500_results(1,percKnown_count) = nanmean(F500_results_v2);
        std_framework_results(1,percKnown_count) = nanstd(framework_results_v2);
        std_random_results(1,percKnown_count) = nanstd(random_results_v2);
        std_F500_results(1,percKnown_count) = nanstd(F500_results_v2);
        
        significance1 = ttest2(framework_results_v2(~isnan(framework_results_v2)),random_results_v2(~isnan(random_results_v2)));
        significance2 = ttest2(framework_results_v2(~isnan(framework_results_v2)),F500_results_v2(~isnan(F500_results_v2)));
        if significance1 == 1 && significance2 == 1 && mean_framework_results(1,percKnown_count) > mean_random_results(1,percKnown_count) && mean_framework_results(1,percKnown_count) > mean_F500_results(1,percKnown_count)
            sig_plot(percKnown_count) = 103;
        else
            sig_plot(percKnown_count) = -30;
        end

        percKnown_count = percKnown_count + 1;
    end
    
    subplot(2,3,subplot_count);
    hold on
    %Solid line for combined approach, dashed line for optimization
    %approach 
%     if strcmp(approach_string,'') 
%         lineStyle = ':';
%         errorbar([0 20 40 60 80 100],100*mean_framework_results,100*std_framework_results./sqrt(numRuns),'LineWidth',1.5,'LineStyle',lineStyle,'Color',[0, 0.4470, 0.7410]);
%     elseif strcmp(approach_string,'KE_')
        lineStyle = '-';
        errorbar([0 20 40 60 80 100],100*mean_framework_results,100*std_framework_results./sqrt(numRuns),'LineWidth',1.5,'LineStyle',lineStyle,'Color',[0, 0.4470, 0.7410]);
        errorbar([0 20 40 60 80 100],100*mean_random_results,100*std_random_results./sqrt(numRuns),'LineWidth',1.5,'LineStyle',lineStyle,'Color',[0.8500, 0.3250, 0.0980]);
        errorbar([0 20 40 60 80 100],100*mean_F500_results,100*std_F500_results./sqrt(numRuns),'LineWidth',1.5,'LineStyle',lineStyle,'Color',[0.9290, 0.6940, 0.1250]);
        plot([0 20 40 60 80 100],sig_plot,'k*')
%     end
    
    

    xlim([-10 110])
    ylim([0 105])


    if subplot_count == 1
        legend('MetaboPAC','Random','500','Location','Northwest');
    end

    title({'Determined',sprintf('< log_2(%.1f) error',log2_error)})
    xlabel('% of kinetic equations known');
    ylabel({'% of predicted response','factors within error range'});
    set(gca,'FontSize',14)
    subplot_count = subplot_count + 1;
end
% end
% Kinetic equation approach vs optimization approach plot
f = figure(2);
f.Position = [500 500 1600 800];
annotation('line', [0.05 0.95], [0.5, 0.5])
a = annotation('textbox', [0.05 0.97 0 0],'String',{'A.'},'EdgeColor','none');
a.FontSize = 24;
a.FontWeight = 'bold';
b = annotation('textbox', [0.05 0.5 0 0],'String',{'B.'},'EdgeColor','none');
b.FontSize = 24;
b.FontWeight = 'bold';
subplot_count = 1;
for log2_error = [1.1 1.3 1.5]
    
    percKnown_count = 1;
    mean_kinetic_results = [];
    mean_optimization_results = [];
    std_kinetic_results = [];
    std_optimization_results = [];
    numRuns = nan(1,size(abs_log2_error,2));
    for percKnownFlux = [0 20 40 60 80 100]
        abs_log2_kinetic_error_results = [abs_log2_error_kinetic{:,percKnown_count}];
        abs_log2_optimization_error_results = [abs_log2_error_optimization{:,percKnown_count}];
        % manage the dummy 10 replicates of different subsets of known
        % kinetic equations for 0% and 100% known kinetics 
        % remove those and fill in nan values
        if isequal(percKnownFlux,0) || isequal(percKnownFlux,100)
            abs_log2_kinetic_error_results = abs_log2_kinetic_error_results(1:10:end);
            abs_log2_optimization_error_results = abs_log2_optimization_error_results(1:10:end);
        end
        numRuns(percKnown_count) = (length(abs_log2_kinetic_error_results) + length(abs_log2_optimization_error_results))/numMetabs;
        mean_kinetic_results = [mean_kinetic_results mean(abs_log2_kinetic_error_results < log2(log2_error),'all')];
        mean_optimization_results = [mean_optimization_results mean(abs_log2_optimization_error_results < log2(log2_error),'all')];
        std_kinetic_results = [std_kinetic_results std(abs_log2_kinetic_error_results < log2(log2_error),[],'all')];
        std_optimization_results = [std_optimization_results std(abs_log2_optimization_error_results < log2(log2_error),[],'all')];
        
        percKnown_count = percKnown_count + 1;
    end
   
    subplot(2,3,subplot_count);
    hold on
    errorbar([0 20 40 60 80 100],100*mean_kinetic_results,100*std_kinetic_results./sqrt(numRuns),'LineWidth',1.5);
    errorbar([0 20 40 60 80 100],100*mean_optimization_results,100*std_optimization_results./sqrt(numRuns),'LineWidth',1.5);

    xlim([-10 110])
    ylim([0 105])


    if subplot_count == 1
        legend('Kinetic Equation','Optimization','Location','Northwest');
    end

    title({'Determined',sprintf('< log_2(%.1f) error',log2_error)})
    %ylabel({'Percent of predicted conversion';'factors within \chi% error'});
    xlabel('% of kinetic equations known');
    ylabel({'% of predicted response','factors within error range'});
    set(gca,'FontSize',14)
    subplot_count = subplot_count + 1;
end

% Percent of response factors found using kinetic equations approach
f = figure(3);
f.Position = [500 500 1400 600];
annotation('line', [0.5 0.5], [0.05, 0.95])
a = annotation('textbox', [0.05 0.99 0 0],'String',{'A.'},'EdgeColor','none');
a.FontSize = 24;
a.FontWeight = 'bold';
b = annotation('textbox', [0.55 0.99 0 0],'String',{'B.'},'EdgeColor','none');
b.FontSize = 24;
b.FontWeight = 'bold';
subplot(1,2,1)
hold on;
errorbar([0 20 40 60 80 100],100*(mean(num_knownMet)/numMetabs),100*(std(num_knownMet)/numMetabs)./sqrt(numRuns),'LineWidth',1.5)
errorbar([0 20 40 60 80 100],100*(mean(numMetabs-num_knownMet)/numMetabs),100*(std(numMetabs-num_knownMet)/numMetabs)./sqrt(numRuns),'LineWidth',1.5)
title('Determined')
ylabel('Percent of response factors found in each approach')
xlabel('% of kinetic equations known')
legend('Kinetic Equations','Optimization')
xlim([-10 110])
ylim([-5 105])
set(gca,'FontSize',14)


% for string_idx = 1:1:length(approach_strings) 
%     approach_string = approach_strings{string_idx};
for log2_error = [1.1 1.3 1.5]
    percKnown_count = 1;
    for percKnownFlux = [0 20 40 60 80 100]
        run_count = 1;
        for run = 1:5
            for rand_idx = 1:10
                run_count = 10 * (run-1) + rand_idx;
                fileName = sprintf('results/UDreg_MetaboPAC_%saddKinetics-%03d_rep-%03d_rand-%03d.mat',approach_string,percKnownFlux,run,rand_idx);

                if isfile(fileName)
                    load(fileName,'RF_opt','trueRF','numMetabs','maxRandVal','minRandVal','knownMet');

                    med_F = median(RF_opt(:,1:numMetabs));
                    rand_F = (maxRandVal-minRandVal).*rand(1,numMetabs) + minRandVal;
                    F500 = 500*ones(1,numMetabs);

                    num_knownMet(run_count,percKnown_count) = length(knownMet);

                    abs_log2_error{run_count,percKnown_count} = abs(log2(med_F) - log2(trueRF));
                    abs_log2_error_kinetic{run_count,percKnown_count} = abs(log2(med_F(:,knownMet)) - log2(trueRF(knownMet)));
                    abs_log2_error_optimization{run_count,percKnown_count} = abs(log2(med_F(:,setdiff(1:numMetabs,knownMet))) - log2(trueRF(setdiff(1:numMetabs,knownMet))));

                    abs_log2_error_random{run_count,percKnown_count} = abs(log2(rand_F) - log2(trueRF));

                    abs_log2_error_F500{run_count,percKnown_count} = abs(log2(F500) - log2(trueRF));

                end
            end
        end
        percKnown_count = percKnown_count + 1;
    end
end

% Log2 error plot
figure(1);
subplot_count = 4;
for log2_error = [1.1 1.3 1.5]
    
    percKnown_count = 1;
    mean_framework_results = [];
    mean_random_results = [];
    mean_F500_results = [];
    std_framework_results = [];
    std_random_results = [];
    std_F500_results = [];
    numRuns = nan(1,size(abs_log2_error,2));
    for percKnownFlux = [0 20 40 60 80 100]
        for run = 1:size(abs_log2_error,1)
            framework_results(run,percKnown_count) = mean(abs_log2_error{run,percKnown_count} < log2(log2_error));
            random_results(run,percKnown_count) = mean(abs_log2_error_random{run,percKnown_count} < log2(log2_error));
            F500_results(run,percKnown_count) = mean(abs_log2_error_F500{run,percKnown_count} < log2(log2_error));
            
            
        end
        % manage the dummy 10 replicates of different subsets of known
        % kinetic equations for 0% and 100% known kinetics 
        % remove those and fill in nan values
        if isequal(percKnownFlux,0) || isequal(percKnownFlux,100)
            framework_results_v2 = nan(size(framework_results(:,percKnown_count)));
            framework_results_v2(1:10:end) = framework_results(1:10:end,percKnown_count);
            
            random_results_v2 = nan(size(random_results(:,percKnown_count)));
            random_results_v2(1:10:end) = random_results(1:10:end,percKnown_count);
            
            F500_results_v2 = nan(size(F500_results(:,percKnown_count)));
            F500_results_v2(1:10:end) = F500_results(1:10:end,percKnown_count);
        else
            framework_results_v2 = framework_results(:,percKnown_count);
            random_results_v2 = random_results(:,percKnown_count);
            F500_results_v2 = F500_results(:,percKnown_count);
        end
        
        numRuns(percKnown_count) = sum(~isnan(framework_results_v2));
        mean_framework_results(1,percKnown_count) = nanmean(framework_results_v2);
        mean_random_results(1,percKnown_count) = nanmean(random_results_v2);
        mean_F500_results(1,percKnown_count) = nanmean(F500_results_v2);
        std_framework_results(1,percKnown_count) = nanstd(framework_results_v2);
        std_random_results(1,percKnown_count) = nanstd(random_results_v2);
        std_F500_results(1,percKnown_count) = nanstd(F500_results_v2);
        
        significance1 = ttest2(framework_results_v2(~isnan(framework_results_v2)),random_results_v2(~isnan(random_results_v2)));
        significance2 = ttest2(framework_results_v2(~isnan(framework_results_v2)),F500_results_v2(~isnan(F500_results_v2)));
        if significance1 == 1 && significance2 == 1 && mean_framework_results(1,percKnown_count) > mean_random_results(1,percKnown_count) && mean_framework_results(1,percKnown_count) > mean_F500_results(1,percKnown_count)
            sig_plot(percKnown_count) = 103;
        else
            sig_plot(percKnown_count) = -30;
        end

        percKnown_count = percKnown_count + 1;
    end
        
    
    subplot(2,3,subplot_count);
    hold on
    %Solid line for combined approach, dashed line for optimization
    %approach 
%     if strcmp(approach_string,'') 
%         lineStyle = ':';
%         errorbar([0 20 40 60 80 100],100*mean_framework_results,100*std_framework_results./sqrt(numRuns),'LineWidth',1.5,'LineStyle',lineStyle,'Color',[0, 0.4470, 0.7410]);
%     elseif strcmp(approach_string,'KE_')
        lineStyle = '-';
        errorbar([0 20 40 60 80 100],100*mean_framework_results,100*std_framework_results./sqrt(numRuns),'LineWidth',1.5,'LineStyle',lineStyle,'Color',[0, 0.4470, 0.7410]);
        errorbar([0 20 40 60 80 100],100*mean_random_results,100*std_random_results./sqrt(numRuns),'LineWidth',1.5,'LineStyle',lineStyle,'Color',[0.8500, 0.3250, 0.0980]);
        errorbar([0 20 40 60 80 100],100*mean_F500_results,100*std_F500_results./sqrt(numRuns),'LineWidth',1.5,'LineStyle',lineStyle,'Color',[0.9290, 0.6940, 0.1250]);
        plot([0 20 40 60 80 100],sig_plot,'k*')
%     end

    xlim([-10 110])
    ylim([0 105])

    if subplot_count == 1
        legend('MetaboPAC','Random','500','Location','Northwest');
    end

    title({'Underdetermined with Regulation',sprintf('< log_2(%.1f) error',log2_error)})
    xlabel('% of additional kinetic equations known');
    ylabel({'% of predicted response','factors within error range'});
    set(gca,'FontSize',14)
    subplot_count = subplot_count + 1;
end
% end
% Kinetic equation approach vs optimization approach plot
figure(2);
subplot_count = 4;
for log2_error = [1.1 1.3 1.5]
    
    percKnown_count = 1;
    mean_kinetic_results = [];
    mean_optimization_results = [];
    std_kinetic_results = [];
    std_optimization_results = [];
    numRuns = nan(1,size(abs_log2_error,2));
    for percKnownFlux = [0 20 40 60 80 100]
        abs_log2_kinetic_error_results = [abs_log2_error_kinetic{:,percKnown_count}];
        abs_log2_optimization_error_results = [abs_log2_error_optimization{:,percKnown_count}];
        % manage the dummy 10 replicates of different subsets of known
        % kinetic equations for 0% and 100% known kinetics 
        % remove those and fill in nan values
        if isequal(percKnownFlux,0) || isequal(percKnownFlux,100)
            abs_log2_kinetic_error_results = abs_log2_kinetic_error_results(1:10:end);
            abs_log2_optimization_error_results = abs_log2_optimization_error_results(1:10:end);
        end
        numRuns(percKnown_count) = (length(abs_log2_kinetic_error_results) + length(abs_log2_optimization_error_results))/numMetabs;
        
        mean_kinetic_results = [mean_kinetic_results mean(abs_log2_kinetic_error_results < log2(log2_error),'all')];
        mean_optimization_results = [mean_optimization_results mean(abs_log2_optimization_error_results < log2(log2_error),'all')];
        std_kinetic_results = [std_kinetic_results std(abs_log2_kinetic_error_results < log2(log2_error),[],'all')];
        std_optimization_results = [std_optimization_results std(abs_log2_optimization_error_results < log2(log2_error),[],'all')];
        
        percKnown_count = percKnown_count + 1;
    end
        
    numRuns = size(abs_log2_error,1);
    subplot(2,3,subplot_count);
    hold on
    
    errorbar([0 20 40 60 80 100],100*mean_kinetic_results,100*std_kinetic_results./sqrt(numRuns),'LineWidth',1.5);
    errorbar([0 20 40 60 80 100],100*mean_optimization_results,100*std_optimization_results./sqrt(numRuns),'LineWidth',1.5);

    xlim([-10 110])
    ylim([0 105])


    if subplot_count == 1
        legend('Kinetic Equation','Optimization','Location','Northwest');
    end

    title({'Underdetermined with Regulation',sprintf('< log_2(%.1f) error',log2_error)})
    %ylabel({'Percent of predicted conversion';'factors within \chi% error'});
    xlabel('% of additional kinetic equations known');
    ylabel({'% of predicted response','factors within error range'});
    set(gca,'FontSize',14)
    subplot_count = subplot_count + 1;
end

% Percent of response factors found using kinetic equations approach
figure(3);
subplot(1,2,2)
hold on;
errorbar([0 20 40 60 80 100],100*(mean(num_knownMet)/numMetabs),100*(std(num_knownMet)/numMetabs)/sqrt(numRuns),'LineWidth',1.5)
errorbar([0 20 40 60 80 100],100*(mean(numMetabs-num_knownMet)/numMetabs),100*(std(numMetabs-num_knownMet)/numMetabs)/sqrt(numRuns),'LineWidth',1.5)
title('Underdetermined with Regulation')
ylabel('Percent of response factors found in each approach')
xlabel('% of additional kinetic equations known')
xlim([-10 110])
ylim([-5 105])
set(gca,'FontSize',14)