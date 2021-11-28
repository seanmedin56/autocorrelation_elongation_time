% script finds the minimum of the third derivative for a range of cleaning
% parameters for a given experimental data set. Results are analyzed in the
% live script "exp_clean_analysis.mlx"

% makes folder for figures if it doesn't exist already
if ~exist('./../dat/exp_clean_res', 'dir')
    mkdir('./../dat/exp_clean_res');
end

% autocorrelaiton parameters
max_delay = 20;
cut = 0;

% range of cleaning/dirtying parameters
threshes = 0;
jump_percs = [97, 98, 99, 99.5, 99.9];
min_valids = 25:5:80;
data_type = 'fluo3DN'; % options: 'fluo3D2', 'fluo3DN', 'fluo3DNRaw'

max_del_z = 3;
max_del_xy = 8;
f = 1.5;

% imports experimental data to analyze
raw_data = load('./../dat/3_4kb_21_04_24/nucleus_struct.mat');
%raw_data = load('./../dat/2_65kb_21_04_24/nucleus_struct.mat');
%raw_data = load('./../dat/1_9kb_21_04_24/nucleus_struct.mat');

nc = 14;
nucleus_struct = raw_data.nucleus_struct;
nucleus_struct = nucleus_struct([nucleus_struct.setID] >= 0);
subset = nucleus_struct([nucleus_struct.ncStart] == nc);

% figures out time resolution
time_res = nanmedian(diff([subset.time]));

all_data = struct;

% loop through combinations
data_idx = 1;
% loops through combinations of cleaning/dirtying parameters
for thresh = threshes
    for jump_perc = jump_percs
        for min_valid = min_valids
            num_big = 0;
            % figures out threshhold for fluo change
            jump_thresh = prctile(abs(diff(diff([subset.(data_type)]))), jump_perc);

            % interpolates and divides traces
            tracesp = {};
            idx = 1;
            idx2 = 1;
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
                del_z = abs(diff(zPos));
                tr_del_xy = [0 max(del_xy(1:end-1), del_xy(2:end)) 0];
                fluo(tr_del_xy >= max_del_xy) = NaN; 

                tr_dd1 = abs([0 diff(diff(fluo)) 0]);
                fluo(tr_dd1>jump_thresh) = NaN;
                num_big = num_big + length(find(tr_dd1 > jump_thresh));

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
                        tracesp{idx} = new_tr;
                        idx = idx + 1;
                    end
                end
                
                % processes traces for the nan autocorrelation methods
                trace = interp1(time, fluo, time_interp);

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
                    if bad_pts(j) - bad_pts(j-1) <= min_valid
                        trace(bad_pts(j-1)+1:bad_pts(j)-1) = nan;

                    end
                end
            end
            
            %--------- runs autocorrelation-------------------%
            auto = calc_auto_corr(tracesp, max_delay);
            deriv3 = diff(diff(diff(auto)));

            % smoothing
            deriv3_smooth = deriv3;
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

            all_data(data_idx).time_res = time_res;
            all_data(data_idx).thresh = thresh;
            all_data(data_idx).jump_perc = jump_perc;
            all_data(data_idx).min_valid = min_valid;
            all_data(data_idx).min3 = min_idx;
            all_data(data_idx).min3_s = min_idx_s;
            all_data(data_idx).num_traces = length(tracesp);
            data_idx = data_idx + 1;
        end
    end
end     

save('./../dat/exp_clean_res/3_4kb_fluo3DN_06_13_21.mat', 'all_data');

