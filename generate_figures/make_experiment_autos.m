% plots autocorrelation and derivatives of traces from 1.9kb, 2.65kb, and
% 3.4kb experimental data along with a boostrapped distribution of 3rd
% derivative global minimum

% makes folder for figures if it doesn't exist already
if ~exist('./../fig/experimental_system', 'dir')
    mkdir('./../fig/experimental_system');
end

% data sets
names = {'3.4kb', '2.65kb', '1.9kb'};
file_names = {'3_4kb', '2_65kb', '1_9kb'};
addresses = {'./../dat/3_4kb_21_04_24/nucleus_struct.mat',...
    './../dat/2_65kb_21_04_24/nucleus_struct.mat', ...
    './../dat/1_9kb_21_04_24/nucleus_struct.mat'};

% parameters
min_valid = 40;
max_del_xy = 8;
f = 1.5;
thresh = 99.5;
nc = 14;
data_type = 'fluo3DN';
max_delay = 20;

% loops through stuff
for n = 1:length(names)
    raw_data = load(addresses{n});
    nucleus_struct = raw_data.nucleus_struct;
    subset = nucleus_struct([nucleus_struct.ncStart] == nc);
    
    jump_thresh = prctile(abs(diff(diff([subset.(data_type)]))), thresh);
    time_res = nanmedian(diff([subset.time]));
    traces = {};
    idx = 1;
    for i = 1:length(subset)
        fluo = subset(i).(data_type);
        time = subset(i).time;
        xPos = subset(i).xPosParticle;
        yPos = subset(i).yPosParticle;

        % figures out timing for interpolation
        t_start = time(1);
        t_end = round((time(end) - t_start) ...
            / time_res) * time_res + t_start;
        time_interp = t_start:time_res:t_end;

        % nans too big movements
        del_xy = sqrt(diff(xPos).^2 + diff(yPos).^2);
        tr_del_xy = [0 max(del_xy(1:end-1), del_xy(2:end)) 0];
        fluo(tr_del_xy >= max_del_xy) = NaN;  

        tr_dd1 = abs([0 diff(diff(fluo)) 0]);
        fluo(tr_dd1>jump_thresh) = NaN;

        % interpolates
        non_nans = find(~isnan(fluo));
        if length(non_nans) < min_valid
            continue
        end
        trace = interp1(time(non_nans), ...
            fluo(non_nans), time_interp);

        % nans below thresh
        trace(trace < thresh) = NaN;

        % nans anything not close to known point
        new_times = time(non_nans);
        for j = 1:length(time_interp)
            higher = find(new_times > time_interp(j));
            lower = find(new_times < time_interp(j));
            if isempty(higher) || isempty(lower)
                continue
            end
            dist_high = new_times(higher(1)) - time_interp(j);
            dist_low = time_interp(j) - new_times(lower(end));
            if dist_low > f*time_res && dist_high > f*time_res
                trace(j) = NaN;
            end
        end

        % divide into segments
        bad_pts = [1 find(isnan(trace)) length(trace)];
        for j = 2:length(bad_pts)
            if bad_pts(j) - bad_pts(j-1) > min_valid
                new_tr = trace(bad_pts(j-1)+1:bad_pts(j)-1);
                traces{idx} = new_tr;
                idx = idx + 1;
            end
        end
    end
    
    %----calculates/plots 2nd, 3rd derivatives of the autocorrelation------
    auto = calc_auto_corr(traces, max_delay);
    deriv1 = diff(auto);
    deriv2 = diff(deriv1);
    deriv3 = diff(deriv2);
    auto_stds = corr_bootstraps(traces, max_delay, 100, 3);
    
    % smoothing
    deriv3_smooth = deriv3;
    for i = 4:(length(deriv3) - 1)
        deriv3_smooth(i) = 0.25*deriv3(i-1) + ...
            0.5*deriv3(i) + 0.25*deriv3(i+1);
    end
    deriv3_smooth(end) = (deriv3_smooth(end-1) + ...
        deriv3_smooth(end) * 2) / 3;
    
    [M,I] = min(deriv3_smooth(4:end));
    glob_min = I(1) + 3;
    
    % plots autocorrelation
    figure('DefaultAxesFontSize',10)
    errorbar([0:(length(auto)-1)] * time_res, auto, auto_stds{1}, 'o-');
    xline(glob_min*time_res, '-', 'Estimated Elongation Time', ...
        'LabelVerticalAlignment', 'bottom');
    xlabel('time delay (seconds)');
    ylabel('correlation');
    grid on
    saveas(gcf, ['./../fig/experimental_system/auto_' file_names{n} '.svg']);
    
    % plots 1st derivative
    figure('DefaultAxesFontSize',10)
    errorbar([1:length(deriv1)] * time_res, deriv1, auto_stds{2}, 'o-');
    xline(glob_min*time_res, '-', 'Estimated Elongation Time', ...
        'LabelVerticalAlignment', 'bottom');
    xlabel('time delay (seconds)');
    ylabel('\Delta correlation');
    grid on
    saveas(gcf, ['./../fig/experimental_system/deriv1_' file_names{n} '.svg']);
    
    % plots 2nd derivative
    figure('DefaultAxesFontSize',10)
    errorbar([1:length(deriv2)] * time_res, deriv2, auto_stds{3}, 'o-');
    xline(glob_min*time_res, '-', 'Estimated Elongation Time', ...
        'LabelVerticalAlignment', 'bottom');
    xlabel('time delay (seconds)');
    ylabel('\Delta^2 correlation');
    grid on
    saveas(gcf, ['./../fig/experimental_system/deriv2_' file_names{n} '.svg']);
    
    % plots 3rd derivative
    figure('DefaultAxesFontSize',10)
    errorbar([3:length(deriv3)] * time_res, deriv3(3:end), auto_stds{4}(3:end), 'o-');
    hold on
    plot([3:length(deriv3)] * time_res, deriv3_smooth(3:end), 'o-');
    xline(glob_min*time_res, '-', 'Estimated Elongation Time');
    xlabel('time delay (seconds)');
    ylabel('\Delta^3 correlation');
    legend({'Original', 'Smoothed'});
    grid on
    saveas(gcf, ['./../fig/experimental_system/deriv3_' file_names{n} '.svg']);
    
    %------------- bootstraps 3rd derivative global minimum---------------
    mins_lst = zeros(1, 500);
    for j = 1:length(mins_lst)
        sample_idx = randi([1 length(traces)], 1, length(traces));
        sample = cell([1 length(traces)]);
        for k = 1:length(sample)
            sample{k} = traces{sample_idx(k)};
        end
        auto = calc_auto_corr(sample, max_delay);
        deriv3 = diff(diff(diff(auto)));
        
        deriv3_smooth = deriv3;
        for i = 4:(length(deriv3) - 1)
            deriv3_smooth(i) = 0.25*deriv3(i-1) + ...
                0.5*deriv3(i) + 0.25*deriv3(i+1);
        end
        deriv3_smooth(end) = (deriv3_smooth(end-1) + ...
            deriv3_smooth(end) * 2) / 3;

        [M,I] = min(deriv3_smooth(4:end));
        temp_min = I(1) + 3;
        mins_lst(j) = temp_min;
    end
    options = unique(mins_lst);
    vals = zeros(1, length(options));
    for j = 1:length(options)
        vals(j) = length(find(mins_lst == options(j))) / length(mins_lst);
    end
    figure();
    plot(options * time_res, vals, 'o-');
    xline(glob_min * time_res, '-', 'Estimated Elongation Time', ...
        'LabelVerticalAlignment', 'bottom');
    xlabel('Global Minimum of 3rd Derivative');
    ylabel('Probability');
    grid on
    saveas(gcf, ['./../fig/experimental_system/min_dist_' file_names{n} '.svg']);
end
