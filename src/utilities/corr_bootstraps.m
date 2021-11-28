function std_derivs = corr_bootstraps(traces, max_delay, num_times, ...
    num_derivs)

% takes random selections of traces and calculates the standard deviation
% of the bootstrapped autocorrelation and its derivatives

%   max_delay: number of time delay points to take in the auto correlation
%   num_times: number of times to sample the traces
%   num_derivs: number of derivatives to return

    vals = zeros(num_times, max_delay);
    std_derivs = cell(1, num_derivs + 1);
    
    %iterates through random subsamples of traces
    for i = 1:num_times
        sample_idx = randi([1 length(traces)], 1, length(traces));
        samples = cell([1 length(traces)]);
        for j = 1:length(samples)
            samples{j} = traces{sample_idx(j)};
        end

        corr = calc_auto_corr(samples, max_delay);

        vals(i,:) = corr;
    end
    for deriv = 0:num_derivs
        new_vals = vals;
        for idx = 1:deriv
            new_vals = new_vals(:,2:end) - new_vals(:,1:end-1);
        end
        std_derivs{deriv + 1} = std(new_vals, 0, 1);
    end
end