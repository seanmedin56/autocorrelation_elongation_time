% makes SI plot showing how the combination of noise and interpolation
% affects the autocorrelation and its derivatives for our simplest model of
% the constituitively active Poisson promoter

% makes folder for figures if it doesn't exist already
if ~exist('./../fig/interp_figs', 'dir')
    mkdir('./../fig/interp_figs');
end

% key variables
elong_time = 120;
time_res = 10;
rise_time = 0;
points_per_trace = 300;
rna_per_sec = 0.2; 
noises = [0,750];
num_traces = 500;
fluo_per_rna = 350;
prob_start_on = 1; 
cut = 15;

max_delay = 20;
% loops through noise scenarios
for noise = noises
    % generates simple poisson promoter plots
    traces = cell(1,num_traces);
    traces2 = cell(1,num_traces);
    traces3 = cell(1,num_traces);
    for i = 1:num_traces
        traces{i} = gillespie_gen(elong_time, time_res, points_per_trace, ...
                                     1, 0, rna_per_sec, ...
                                     fluo_per_rna, 0,1,noise);
        traces{i} = traces{i}(1+cut:end);
        
        % deletes and interpolates
        traces2{i} = traces{i};
        traces3{i} = traces{i};
        for j = 1:length(traces{i})
            if rand() < 0.2
                traces2{i}(j) = NaN;
            end
            if rand() < 0.4
                traces3{i}(j) = NaN;
            end
        end
        idxes2 = find(~isnan(traces2{i}));
        idxes3 = find(~isnan(traces3{i}));
        traces2{i} = interp1(idxes2, traces2{i}(idxes2), 1:length(traces{i}));
        traces3{i} = interp1(idxes3, traces3{i}(idxes3), 1:length(traces{i}));
        
    end
    

    % simulated autocorrelation
    auto = calc_auto_corr(traces, max_delay);
    auto2 = calc_auto_corr(traces2, max_delay);
    auto3 = calc_auto_corr(traces3, max_delay);
    
    % simulates autocorrelation bootstraps
    auto_stds = corr_bootstraps(traces, max_delay, 100, 0);
    auto_stds2 = corr_bootstraps(traces2, max_delay, 100, 0);
    auto_stds3 = corr_bootstraps(traces3, max_delay, 100, 0);

    figure('DefaultAxesFontSize',10)
    errorbar(0:length(auto)-1, auto, auto_stds{1}, 'o-');
    hold on
    errorbar(0:length(auto)-1, auto2, auto_stds2{1}, 'o-');
    hold on
    errorbar(0:length(auto)-1, auto3, auto_stds3{1}, 'o-');

    xline(12, '-', 'Elongation Time');
    xlabel('time delay');
    ylabel(['correlation']);
    grid on
    legend({'Normal', '20% Interpolated', '40% Interpolated'});
    saveas(gcf, ['./../fig/interp_figs/auto_noise' int2str(noise) '.svg']);
    %close;

    cur_corr = auto;
    cur_corr2 = auto2;
    cur_corr3 = auto3;
    for i = 1:3
        cur_corr = diff(cur_corr); % calculates latest deriv
        cur_corr2 = diff(cur_corr2);
        cur_corr3 = diff(cur_corr3);
        
        corr_stds = corr_bootstraps(traces, max_delay, 100, i);
        corr_stds2 = corr_bootstraps(traces2, max_delay, 100, i);
        corr_stds3 = corr_bootstraps(traces3, max_delay, 100, i);

        figure('DefaultAxesFontSize',10)
        errorbar(1:length(cur_corr), cur_corr, corr_stds{i+1}, 'o-');
        hold on
        errorbar(1:length(cur_corr), cur_corr2, corr_stds2{i+1}, 'o-');
        hold on
        errorbar(1:length(cur_corr), cur_corr3, corr_stds3{i+1}, 'o-');
        
        xline(12, '-', 'Elongation Time');
        xlabel('time delay');
        ylabel(['\Delta^' int2str(i) ' correlation']);
        grid on
        legend({'Normal', '20% Interpolated', '40% Interpolated'}, 'Location', 'northwest');
        saveas(gcf, ['./../fig/interp_figs/deriv' int2str(i) '_noise' int2str(noise) '.svg']);
        %close;
    end
end