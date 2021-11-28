function hs = gen_auto_cor(traces, derivs, ...
                                    bootstraps, max_delay,cut, varargin)
% Generates an autocorrelation and/or the derivatives of an
% autocorrelation for the traces
% traces: The traces we are taking the autocorrelation of
% auto: Boolean saying whether or not to include the autocorrelation
% first: Boolean determining if to include 1st derivative
% second: Boolean determining if to include 2nd derivative
% bootstraps: Boolean determining if to include bootstraps
% max_delay: How many points to take of the autocorrelation
% returns handles for plots

hs = [];
%cuts the traces by amount cut
for i=1:length(traces)
    traces{i} = traces{i}(1 + cut:end);
end

corr = calc_auto_corr(traces, max_delay);

if bootstraps
    std_derivs = corr_bootstraps(traces, max_delay,100, max(derivs));
end

for deriv = derivs
    h = figure;
    corr_deriv = corr;
    for i = 1:deriv
        corr_deriv = diff(corr_deriv);
    end
    if deriv == 0
       times = 0:length(corr_deriv) - 1;
    else
        times = 1:length(corr_deriv);
    end
    if ~isempty(varargin)
        times = times * varargin{1};
    end
    if bootstraps
        errorbar(times, corr_deriv,std_derivs{deriv + 1}, '-o');
    else
        plot(times, corr_deriv, '-o');
    end
    title(['Central Moment Derivative ' num2str(deriv)], 'FontSize', 14);
    xlabel('time delay', 'FontSize', 14);
    if deriv == 0
        ylabel('Correlation')
    else
        ylabel(['\Delta^' int2str(deriv) 'Correlation']); 
    end
    grid on
    hs(deriv + 1) = h;
end
hs(~ismember(1:length(hs), derivs + 1)) = -1;
