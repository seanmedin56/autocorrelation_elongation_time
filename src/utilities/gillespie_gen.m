function [trace, arrival_times] = gillespie_gen(elong_time, time_res, points_per_trace, ...
                                 num_states, trans_mat, rna_per_sec, ...
                                 fluo_per_rna, MS2_rise_time,init_dist,noise)
%Generates individual traces with the gillespie algorithm
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
% noise: Guassian noise standard deviation

%temp

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

% uses the arrival times to generate the traces
for k = 1:points_per_trace
    t_end = times_unif(k);
    t_start = max([0, t_end - elong_time]);

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

    times_window = arrival_times(i_start:i_end);

    for i = 1:(length(times_window))

        t2 = t_end - times_window(i);
        if t2 > MS2_rise_time
            trace(k) = trace(k) + fluo_per_rna;
        else
            trace(k) = trace(k) + fluo_per_rna * ...
                      t2 / MS2_rise_time;
        end
    end
end

% adds Gaussian noise
gauss_noise = normrnd(0,noise,1,points_per_trace);
trace = trace + gauss_noise;
trace(trace < 0) = 0;

end

