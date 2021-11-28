% script simulates fluorescence traces for a range of parameters
% (elongation times, rise times, dynamics) and tests if the autocorrelation
% still has expected features when we use the data cleaning methods we use
% for experimental data. Data is analyzed in live script "exp_clean_analysis.mlx"

% makes folder if it doesn't exist already
if ~exist('./../dat/sim_clean_res', 'dir')
    mkdir('./../dat/sim_clean_res');
end

% range of intrinstic parameters
Ts = [6:15];
max_delay = 20;
aes = [0:0.25:1];
num_states = 2;
num_traces = 400;
k_ons = [0.001, 0.01, 0.02, 0.05, 0.1];
rna_per_sec = [.0001, 0.34];
p0 = [0,1];
time_res = 10;
fluo_per_rna = 350;
total_time = 60 * 40;
points_per_trace = total_time / time_res;
time_interp = [1:points_per_trace] * time_res;
f = 1.5;

% adjusted spread of time points collected
dt_spread1 = [time_res, 0.4];
dt_spread2 = [15, 5];
p_spread2 = .03;

% range of cleaning/dirtying parameters
noises = [0, 500, 1000, 1500];
threshes = 1000;
jump_percs = [97, 98, 99, 99.5, 99.9];
p_nans = 0.02;
min_valids = [25:5:50];

all_data = struct;

% loop through combinations
data_idx = 1;
for T = Ts
    for a0 = aes
        for k_on = k_ons
            % create transition rate matrix
            R = [-k_on, k_on;
                 k_on, -k_on];

            % generates underlying polymerase data
            base_traces = cell(1, num_traces);
            base_times = cell(1, num_traces);
            for i = 1:num_traces
                [trace, arrival_times] = gillespie_gen(T * time_res, time_res, ...
                    points_per_trace, num_states, R, rna_per_sec, ...
                    fluo_per_rna, a0 * T * time_res,p0,0);
                % uses the arrival times to generate the traces
                t_cur = time_res;
                trace2 = [];
                time2 = [];
                while t_cur < total_time
                    t_end = t_cur;
                    t_start = max([0, t_end - T * time_res]);
                    new_val = 0;
                    ind_start_before = find(arrival_times >= t_start);
                    if ~isempty(ind_start_before)
                        i_start = ind_start_before(1);

                        ind_end_before = find(arrival_times <= t_end);

                        if ~isempty(ind_end_before)

                            i_end = ind_end_before(end);

                            times_window = arrival_times(i_start:i_end);

                            for ii = 1:(length(times_window))

                                t2 = t_end - times_window(ii);
                                if t2 > a0 * T * time_res
                                    new_val = new_val + fluo_per_rna;
                                else
                                    new_val = new_val + fluo_per_rna * ...
                                              t2 / (a0 * T * time_res);
                                end
                            end
                        end
                    end
                    trace2 = [trace2 new_val];
                    time2 = [time2, t_cur];
                    time_pass = normrnd(dt_spread1(1), dt_spread1(2));
                    if rand < p_spread2
                        time_pass = time_pass + normrnd(dt_spread2(1),dt_spread2(2));
                    end
                    t_cur = time_pass + t_cur;
                end
                base_traces{i} = trace2;
                base_times{i} = time2;
            end

            % loops through combinations of cleaning/dirtying parameters
            for noise = noises
                for thresh = threshes
                    for jump_perc = jump_percs
                        for p_nan = p_nans
                            for min_valid = min_valids

                                % adds noise and nans out low values
                                traces_used = cell(1, length(base_traces));
                                for i = 1:length(base_traces)
                                    trace = base_traces{i};
                                    gauss_noise = normrnd(0,noise,1,length(trace));
                                    trace = trace + gauss_noise;
                                    trace(trace < thresh) = NaN;
                                    for j = 1:length(trace)
                                        if rand < p_nan
                                            trace(j) = NaN;
                                        end
                                    end
                                    traces_used{i} = trace;
                                end

                                % figures out threshhold for fluo change
                                all_del_fluo = [];
                                for i = 1:length(traces_used)
                                    all_del_fluo = [all_del_fluo abs(diff(diff(traces_used{i})))];
                                end
                                num_big_jumps = 0;
                                big_jump = prctile(all_del_fluo, jump_perc);

                                % interpolates and divides traces
                                fin_traces = {};
                                idx = 1;
                                for i = 1:length(traces_used)
                                    if isempty(~isnan(traces_used{i}))
                                        continue
                                    end
                                    non_nans = find(~isnan(traces_used{i}));
                                    if length(non_nans) < min_valid
                                        continue
                                    end
                                    trace = interp1(base_times{i}(non_nans), ...
                                        traces_used{i}(non_nans), time_interp);
                                    new_times = base_times{i}(non_nans);
                                    for j = 1:length(time_interp)
                                        higher = find(new_times > time_interp(j));
                                        lower = find(new_times < time_interp(j));
                                        if isempty(higher) || isempty(lower)
                                            continue
                                        end
                                        dist_high = new_times(higher(1)) - time_interp(j);
                                        dist_low = time_interp(j) - new_times(lower(end));
                                        if dist_low > time_res*f && dist_high > time_res*f
                                            trace(j) = NaN;
                                        end
                                    end
                                    trace(trace < thresh) = NaN;

                                    % axes any changes that are too big
                                    tr_dd1 = abs([0 diff(diff(trace)) 0]);
                                    trace(tr_dd1>big_jump) = NaN;
                                    num_big_jumps = num_big_jumps + length(find(tr_dd1>big_jump));

                                    bad_pts = [1 find(isnan(trace)) length(trace)];
                                    for j = 2:length(bad_pts)
                                        if bad_pts(j) - bad_pts(j-1) > min_valid
                                            new_tr = trace(bad_pts(j-1):bad_pts(j));
                                            fin_traces{idx} = new_tr;
                                            idx = idx + 1;
                                        end
                                    end
                                end

                                % runs autocorrelation
                                auto = calc_auto_corr(fin_traces, max_delay);
                                deriv1 = diff(auto);
                                deriv2 = diff(deriv1);
                                deriv3 = diff(deriv2);

                                % smoothing
                                deriv3_smooth = deriv3;
                                deriv3_smooth(3) = (deriv3(3) * 2 + deriv3(3)) / 3;
                                for i = 4:(length(deriv3) - 1)
                                    deriv3_smooth(i) = 0.25*deriv3(i-1) + ...
                                        0.5*deriv3(i) + 0.25*deriv3(i+1);
                                end
                                deriv3_smooth(end) = (deriv3_smooth(end-1) + ...
                                    deriv3_smooth(end) * 2) / 3;

                                % needs offset due to noise correlations
                                % created by interpolation
                                [~,min_idx] = min(deriv3(3:end));
                                min_idx = min_idx + 2;

                                [~,min_idx_s] = min(deriv3_smooth(3:end));
                                min_idx_s = min_idx_s + 2;

                                all_data(data_idx).T = T;
                                all_data(data_idx).a0 = a0;
                                all_data(data_idx).R = R;
                                all_data(data_idx).noise = noise;
                                all_data(data_idx).thresh = thresh;
                                all_data(data_idx).jump_perc = jump_perc;
                                all_data(data_idx).p_nan = p_nan;
                                all_data(data_idx).min_valid = min_valid;
                                all_data(data_idx).min3 = min_idx;
                                all_data(data_idx).min3_s = min_idx_s;
                                all_data(data_idx).num_traces = length(fin_traces);
                                data_idx = data_idx + 1;
                            end
                        end
                    end
                end
            end        
        end
    end
end

save('./../dat/sim_clean_res/simulated_cleaning_06_13_21.mat', 'all_data');

