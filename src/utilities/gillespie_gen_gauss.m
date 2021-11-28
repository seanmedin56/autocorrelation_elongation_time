function [trace, arrival_times] = gillespie_gen_gauss(elong_time, time_res, points_per_trace, ...
                                 num_states, trans_mat, rna_per_sec, ...
                                 fluo_per_rna, MS2_rise_time,init_dist,noise,elong_std)
%Generates individual traces with the gillespie algorithm with elongation
%time being drawn from a Gaussian distribution
% elong_time: how long it takes an rna polymerase to traverse the gene
% time_res: How often the fluorescence is measured
% points_per_trace: How many measurements are taken
% num_state: Number of distinct rates of transcription
% trans_mat: The transition flow matrix between distrinct rates of
%            transcription
% rna_per_sec: Vector of how many rnas are loaded per second for each state
% fluo_per_rna: How much fluorescence each rna produces when its fully
% loaded
% MS2_rise_time: How long it takes for a polymerase to fully fluoresce
% init_dist: The initial probability distribution for being in a given
% state
% noise: Standard deviation of the Gaussian noise
% elong_std: Standard deviation of the Gaussian distribution where the
% elongation time for each polymerase is drawn from
%

% the trace that will be returned
trace = zeros(1,points_per_trace);

% uniformly distributed time points with resolution deltaT 
% and length seq_length
times_unif = (1:points_per_trace) * time_res;

%array for keeping track of polymerase arrivals
arrival_times = [];

%keep track of the state when each polymerase came
naive_states = [];

% duration of the simulated process
t_max = points_per_trace * time_res;

% first state obtained by sampling the initial state pmf
naive_states(1) = randsample(1:num_states, 1, true, init_dist);

% variable to keep track of the current reaction time
t = 0;

while (t < t_max)
    
    % determine when transition to another state occurs
    lambda = -trans_mat(naive_states(end),naive_states(end));
    if lambda == 0
        dt = t_max - t;
    else
        dt = min(exprnd(1/lambda), t_max - t);
    end
    t = t + dt;
    
    %determine probability of transitioning to each of the other states
    rates = trans_mat(:,naive_states(end));
    rates(naive_states(end)) = 0;
    if lambda ~= 0
        probs = rates / lambda;

        %determine next state
        naive_states = [naive_states randsample(1:num_states, 1, true, probs)];
    else
        naive_states = [naive_states naive_states(end)];
    end
    %determine how many polymerases arrived before the transition occurred
    avg_time = 1 / rna_per_sec(naive_states(end-1));
    ddt = exprnd(avg_time);
    while ddt < dt
        arrival_times = [arrival_times t - dt + ddt];
        next = exprnd(avg_time);
        ddt = ddt + next;
    end   
end

% array for keeping track for how far along each polymerase is
pol_perc_covered = zeros(1, length(arrival_times));

% standard deviation of elongation at each time step
step_std = sqrt(elong_std^2 * time_res / elong_time);

% uses the arrival times to generate the traces
for k = 1:points_per_trace
    t_end = times_unif(k);
    t_start = max([0, t_end - time_res]);

    ind_start_before = find(arrival_times >= t_start);
    if isempty(ind_start_before)
        continue;
    end
    i_start = ind_start_before(1);

    ind_end_before = find(arrival_times <= t_end);
    if isempty(ind_end_before)
        continue;
    end
    i_end = ind_end_before(end);

    partial_times = t_end - arrival_times(i_start:i_end);
    
    pol_perc_covered(1:i_start-1) = pol_perc_covered(1:i_start-1) + ... 
        max(0,normrnd(time_res, step_std, 1, i_start-1)) / elong_time;
    
    for newie = 1:length(partial_times)
        adj_step_std = sqrt(partial_times(newie) / time_res * step_std^2);
        idx = newie - 1 + i_start;
        pol_perc_covered(idx) = max(0,normrnd(partial_times(newie), adj_step_std)) / elong_time;
    end
    
    pol_perc_covered = min(1, pol_perc_covered);
    
    belowMS2 = pol_perc_covered(pol_perc_covered < MS2_rise_time / elong_time);
    num_above_MS2 = length(pol_perc_covered(pol_perc_covered >= ...
        MS2_rise_time / elong_time & pol_perc_covered < 1));
    trace(k) = trace(k) + (sum(belowMS2 / (MS2_rise_time / elong_time))...
        +num_above_MS2)*fluo_per_rna;

end

% adds Gaussian noise
gauss_noise = normrnd(0,noise,1,points_per_trace);
trace = trace + gauss_noise;
trace(trace < 0) = 0;

end

