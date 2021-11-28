function auto = calc_auto_corr(traces, max_delay)

% -calculates central moment of the autocorrelation of *traces*
% -only calculates time delays up until *max_delay* 
% -normalizes individual correlations and then performs weighted average
% according to number of points per trace

% ----------------------calculates the global means-------------------------

    global_mean = 0;
    count = 0;
    for i = 1:length(traces)
        global_mean = global_mean + nansum(traces{i});
        count = count + sum(~isnan(traces{i}));
    end
    global_mean = global_mean / count;
    
%--------------subtracts means----------

    new_traces = cell([1 length(traces)]);
    for i = 1:length(traces)
        new_traces{i} = traces{i} - global_mean;
    end
    
    traces = new_traces;

% ----------calculates individual correaltions (with weights)----------
   
    corrs = cell([1 length(traces)]);
    counts = zeros([1 max_delay]);
    for i = 1:length(traces)
        limit = min([max_delay, length(traces{i})]);
        len = length(traces{i});
        corr = zeros([1 limit]);
        num_points = zeros([1 limit]);
        for j = 1:limit
            temp_mult = traces{i}(1:len - j + 1) .* traces{i}(j:len);
            corr(j) = nansum(temp_mult);
            num_points(j) = sum(~isnan(temp_mult));
        end
        counts = counts + num_points;
        corr = corr ./ num_points;
        corrs{i} = corr / corr(1);
    end
    
% -----------------combines correaltions together------------------------    

    auto = zeros([1 max_delay]);
    for i = 1:length(corrs)
        auto(1:length(corrs{i})) = auto(1:length(corrs{i})) ...
            + corrs{i} * length(traces{i});
    end
    auto = auto / auto(1); 

end