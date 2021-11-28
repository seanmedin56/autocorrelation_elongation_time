% Calculates autocorrelation and its derivatives for each of our models for
% the paper.

% makes folder for figures if it doesn't exist already
if ~exist('./../fig/big_plot_figs', 'dir')
    mkdir('./../fig/big_plot_figs');
end

% key variables
elong_time = 80;
time_res = 10;
rise_time = 40; % when applicable
points_per_trace = 40 * 60 / time_res;
rna_per_sec = 0.2; % ON state
noise = 0;
num_traces = 400;%
fluo_per_rna = 350;
prob_start_on = 1; 
k_on = 0.01; % when applicable
k_off = 0.05; % when applicable
cut = 10;

max_delay = 20;

% generates simple poisson promoter plots
traces = cell(1,num_traces);
for i = 1:num_traces
    traces{i} = gillespie_gen(elong_time, time_res, points_per_trace, ...
                                 1, 0, rna_per_sec, ...
                                 fluo_per_rna, 0,1,noise);
    traces{i} = traces{i}(1+cut:end);
end

% simulated autocorrelation
auto = calc_auto_corr(traces, max_delay);
auto_stds = corr_bootstraps(traces, max_delay, 100, 0);

% expected autocorrelation
auto_eq = zeros(1, max_delay);
norm_auto = full_func_cor(elong_time / time_res, 0, 0, [], []);
for i = 0:(max_delay-1)
    auto_eq(i + 1) = full_func_cor(elong_time / time_res, ...
        0,i,[],[]) / norm_auto;
end

figure('DefaultAxesFontSize',10)
errorbar(0:length(auto)-1, auto, auto_stds{1}, 'o-');
hold on
plot(0:length(auto)-1, auto_eq);
hold on
%plot(8, auto(9), 'x', 'MarkerEdgeColor', 'Red', ...
%    'MarkerSize', 12, 'LineWidth', 3);
xline(8, '-', 'Elongation Time');
xlabel('time delay');
ylabel(['correlation']);
grid on
legend({'Simulated', 'Expected Value'});
saveas(gcf, './../fig/big_plot_figs/poisson_auto.svg');
%close;

cur_corr = auto;
cur_corr_eq = auto_eq;
for i = 1:3
    cur_corr = diff(cur_corr); % calculates latest deriv
    corr_stds = corr_bootstraps(traces, max_delay, 100, i);
    cur_corr_eq = diff(cur_corr_eq);   
    
    figure('DefaultAxesFontSize',10)
    errorbar(1:length(cur_corr), cur_corr, corr_stds{i+1}, 'o-');
    hold on
    plot(1:length(cur_corr), cur_corr_eq);
    hold on
    %plot(8, cur_corr(8), 'x', 'MarkerEdgeColor', 'Red', ...
    %'MarkerSize', 12, 'LineWidth', 3);
    xline(8, '-', 'Elongation Time');
    xlabel('time delay');
    ylabel(['\Delta^' int2str(i) ' correlation']);
    grid on
    legend({'Simulated', 'Expected Value'});
    saveas(gcf, ['./../fig/big_plot_figs/poisson_deriv' int2str(i) '.svg']);
    %close;
end

% generates rise time plots
traces = cell(1,num_traces);
for i = 1:num_traces
    traces{i} = gillespie_gen(elong_time, time_res, points_per_trace, ...
                                 1, 0, rna_per_sec, ...
                                 fluo_per_rna, rise_time,1,noise);
    traces{i} = traces{i}(1+cut:end);
end

% simulated autocorrelation
auto = calc_auto_corr(traces, max_delay);
auto_stds = corr_bootstraps(traces, max_delay, 100, 0);

% expected autocorrelation
auto_eq = zeros(1, max_delay);
norm_auto = full_func_cor(elong_time / time_res, rise_time / elong_time, 0, [], []);
for i = 0:(max_delay-1)
    auto_eq(i + 1) = full_func_cor(elong_time / time_res, ...
        rise_time / elong_time,i,[],[]) / norm_auto;
end

figure('DefaultAxesFontSize',10)
errorbar(0:length(auto)-1, auto, auto_stds{1}, 'o-');
hold on
plot(0:length(auto)-1, auto_eq);
hold on
%plot([2,6,8], auto([3,7,9]), 'x', 'MarkerEdgeColor', 'Red', ...
%    'MarkerSize', 12, 'LineWidth', 3);
xline(4, '-', 'Rise Time');
xline(8, '-', 'Elongation Time');
xlabel('time delay');
ylabel(['correlation']);
grid on
legend({'Simulated', 'Expected Value'});
saveas(gcf, './../fig/big_plot_figs/rise_auto.svg');
%close;

cur_corr = auto;
cur_corr_eq = auto_eq;
for i = 1:3
    cur_corr = diff(cur_corr); % calculates latest deriv
    corr_stds = corr_bootstraps(traces, max_delay, 100, i);
    cur_corr_eq = diff(cur_corr_eq);   
    
    figure('DefaultAxesFontSize',10)
    errorbar(1:length(cur_corr), cur_corr, corr_stds{i+1}, 'o-');
    hold on
    plot(1:length(cur_corr), cur_corr_eq);
    hold on
    %plot([2,6,8], cur_corr([2,6,8]), 'x', 'MarkerEdgeColor', 'Red', ...
    %    'MarkerSize', 12, 'LineWidth', 3);
    xline(4, '-', 'Rise Time');
    xline(8, '-', 'Elongation Time');
    xlabel('time delay');
    ylabel(['\Delta^' int2str(i) ' correlation']);
    grid on
    legend({'Simulated', 'Expected Value'});
    saveas(gcf, ['./../fig/big_plot_figs/rise_deriv' int2str(i) '.svg']);
    %close;
end

% generates dynamics plots
traces = cell(1,num_traces);
for i = 1:num_traces
    traces{i} = gillespie_gen(elong_time, time_res, points_per_trace, ...
                                 2, [-k_on, k_off; k_on, -k_off], [.001, rna_per_sec], ...
                                 fluo_per_rna, rise_time,[0,1],noise);
    traces{i} = traces{i}(1+cut:end);
end

% simulated autocorrelation
auto = calc_auto_corr(traces, max_delay);
auto_stds = corr_bootstraps(traces, max_delay, 100, 0);

% expected autocorrelation
auto_eq = zeros(1, max_delay);
[aes, bes] = decompose_matrix([-k_on, k_off; k_on, -k_off]*10, [.001, rna_per_sec]*10);
[M,I] = sort(abs(bes));
aes = aes(I(2:end));
bes = abs(bes(I(2:end)));
norm_auto = full_func_cor(elong_time / time_res, rise_time / elong_time, 0, aes,bes);
for i = 0:(max_delay-1)
    auto_eq(i + 1) = full_func_cor(elong_time / time_res, ...
        rise_time / elong_time,i,aes,bes) / norm_auto;
end

figure('DefaultAxesFontSize',10)
errorbar(0:length(auto)-1, auto, auto_stds{1}, 'o-');
hold on
plot(0:length(auto)-1, auto_eq);
hold on
%plot([2,6,8], auto([3,7,9]), 'x', 'MarkerEdgeColor', 'Red', ...
%    'MarkerSize', 12, 'LineWidth', 3);
xline(4, '-', 'Rise Time');
xline(8, '-', 'Elongation Time');
xlabel('time delay');
ylabel(['correlation']);
grid on
legend({'Simulated', 'Expected Value'});
saveas(gcf, './../fig/big_plot_figs/dynam_auto.svg');
%close;

cur_corr = auto;
cur_corr_eq = auto_eq;
for i = 1:3
    cur_corr = diff(cur_corr); % calculates latest deriv
    corr_stds = corr_bootstraps(traces, max_delay, 100, i);
    cur_corr_eq = diff(cur_corr_eq);   
    
    figure('DefaultAxesFontSize',10)
    errorbar(1:length(cur_corr), cur_corr, corr_stds{i+1}, 'o-');
    hold on
    plot(1:length(cur_corr), cur_corr_eq);
    hold on
    %plot([2,6,8], cur_corr([2,6,8]), 'x', 'MarkerEdgeColor', 'Red', ...
    %    'MarkerSize', 12, 'LineWidth', 3);
    xline(4, '-', 'Rise Time');
    xline(8, '-', 'Elongation Time');
    xlabel('time delay');
    ylabel(['\Delta^' int2str(i) ' correlation']);
    grid on
    legend({'Simulated', 'Expected Value'});
    saveas(gcf, ['./../fig/big_plot_figs/dynam_deriv' int2str(i) '.svg']);
    %close;
end

