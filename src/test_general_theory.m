% Loops through combination of elongation times, rise times, and dynamics
% and calculates the expected autocorrelation for those parameters. Checks
% if the third derivative global minimum is truly at the elongation time.
% Generates .mat file that can be used to help check results and generate
% figures in the live script "make_paper_plots.mlx"

% makes folder if it doesn't exist already
if ~exist('./../dat/comp_explore_res', 'dir')
    mkdir('./../dat/comp_explore_res');
end

% establish cases to check
Ts = [4:0.2:60]; % elongation times
kons = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, ...
    10, 20, 50, 100];
koffs = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, ...
    10, 20, 50, 100];
aes = 0:0.04:1; % rise times (as percentage of elongation time)

dists = struct;
idx = 1;

% check all cases for feature
for T = Ts
    for a0 = aes
        for kon = kons
            for koff = koffs
                auto = zeros(1, ceil(T) + 5);
                c = koff / sqrt(2 * (kon^2 + koff^2));
                b = kon + koff;
                for delay =  0:length(auto)-1
                    [cor_tot, p_term, d_term] = full_func_cor(T,a0,delay,1,b);
                    auto(delay + 1) = p_term + c * d_term;
                end

                deriv1 = diff(auto);
                deriv2 = diff(deriv1);
                deriv3 = diff(deriv2);

                dists(idx).T = T;
                dists(idx).a0 = a0;
                dists(idx).b = b;
                dists(idx).kon = kon;
                dists(idx).koff = koff;
                dists(idx).mins = find(islocalmin(deriv3));
                [~, real_min] = min(deriv3);
                dists(idx).abs_min = real_min;
                if length(real_min) > 1
                    stop = true;
                end
                if isempty(dists(idx).mins)
                    stop = true;
                end
                %dists(idx).maxs = find(islocalmax(deriv3));
                idx = idx + 1;
            end
        end
    end
end

save('./../dat/com_explore_res/comp_explorationg_07_10_21.mat', 'dists');
