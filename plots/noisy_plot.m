% Plot results for noisy synthetic and biological systems 
% Figure 7 & 8
clear;clc;

% Figure 7 
% start_idx = 1;
% end_idx = 2; 
% Figure 8
start_idx = 1;
end_idx = 2;


model_name_list = {'Determined','UDreg','ecoli','yeast'};
model_full_name_list = {'Determined','Underdetermined with Regulation','E. coli','S. cerevisiae'};
subplot_count1 = 1;
subplot_count2 = 1;
f = figure(1);
f.Position = [100 100 800 800];
for model_idx = start_idx:end_idx
    model_name = model_name_list{model_idx};
    model_full_name = model_full_name_list{model_idx};
    subplot_count1 = subplot_count2;
    for nt = [50 15]
        for cov = [5 15]
            percKnown_count = 1;
            for percKnownFlux = [0 20 40 60 80 100]
                run = 1;
                for RFrep =  1:1:3
                    for noise_rep = 1:1:3
                        for rand_idx = 1:1:3
                                run_count = (run - 1) * 27 + (RFrep - 1) * 9 + (noise_rep - 1) * 3 + rand_idx;
                                if strcmp(model_name,'Determined')
                                    string_for_fileName = 'percKnownKinetics';
                                    string_for_xlabel = '% of kinetic equations known';
                                else
                                    string_for_fileName = 'addKinetics';
                                    string_for_xlabel = '% of additional kinetic equations known';
                                end
                                
                                fileName = sprintf('results/%s_MetaboPAC_nT-%03d_cov-%02d_%s-%03d_RFrep-%03d_noiserep-%03d_rand-%03d.mat',model_name,nt,cov,string_for_fileName,percKnownFlux,RFrep,noise_rep,rand_idx);

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
                end
                run = run + 1;

                percKnown_count = percKnown_count + 1;
            end

            % Log2 error plot
            figure(1);
            Position = [500 500 800 800];
            log2_error = 1.3;

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
                        framework_results_v2(1:3:end) = framework_results(1:3:end,percKnown_count);

                        random_results_v2 = nan(size(random_results(:,percKnown_count)));
                        random_results_v2(1:3:end) = random_results(1:3:end,percKnown_count);

                        F500_results_v2 = nan(size(F500_results(:,percKnown_count)));
                        F500_results_v2(1:3:end) = F500_results(1:3:end,percKnown_count);
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

                s = subplot(4,2,subplot_count1);
                hold on

                errorbar([0 20 40 60 80 100],100*mean_framework_results,100*std_framework_results./sqrt(numRuns),'LineWidth',1.5);
                errorbar([0 20 40 60 80 100],100*mean_random_results,100*std_random_results./sqrt(numRuns),'LineWidth',1.5);
                errorbar([0 20 40 60 80 100],100*mean_F500_results,100*std_F500_results./sqrt(numRuns),'LineWidth',1.5);
                plot([0 20 40 60 80 100],sig_plot,'k*')

                xlim([-10 110])
                ylim([0 105])


                if subplot_count1 == 1
                    legend('MetaboPAC','Random','500','Location','Northwest');
                end

                title({sprintf('%s',model_full_name),sprintf('< log_2(%.1f) error',log2_error)})
                xlabel(string_for_xlabel);
                if isequal(model_idx,start_idx)
                    ylabel({'% of predicted response','factors within error range'});
                end
                set(gca,'FontSize',12)
                

            a = annotation('textbox', [0.0 0.85 0 0],'String',{'nT=50','CoV=0.05'},'EdgeColor','none');
            b = annotation('textbox', [0.0 0.63 0 0],'String',{'nT=50','CoV=0.15'},'EdgeColor','none');
            c = annotation('textbox', [0.0 0.41 0 0],'String',{'nT=15','CoV=0.05'},'EdgeColor','none');
            d = annotation('textbox', [0.0 0.19 0 0],'String',{'nT=15','CoV=0.15'},'EdgeColor','none');
            a.FontSize = 12;
            a.FontWeight = 'bold';
            b.FontSize = 12;
            b.FontWeight = 'bold';
            c.FontSize = 12;
            c.FontWeight = 'bold';
            d.FontSize = 12;
            d.FontWeight = 'bold';

        subplot_count1 = subplot_count1 + 2;
        end
        
    end
    subplot_count2 = subplot_count2 + 1;
end
