function [cor_tot, p_term, d_term] = full_func_cor(elong,alph_perc,tau,aes,bes)
% calculates autocorrelation at time delay tau
%   elong: elongation time
%   alph: rise time percentage of elongation time
%   tau: time delay
%   aes: coefficiencts for exponential terms
%   bes: decay values for exponential terms
    alph_perc = max(min(1, alph_perc),0.000001);
    alph = alph_perc * elong;
    
    %calculates poisson term
    
    p_term = 0;
    if tau <= min(alph,elong - alph)
        p_term = -2 * alph / 3 + tau^3 / 6 / alph^2 - tau / 2 ...
            - tau^2 / 2 / alph + elong;
    elseif tau > alph && tau < elong - alph
        p_term = elong - alph / 2 - tau;
    elseif tau > elong - alph && tau < alph
        p_term = -alph /6 + elong^2 / 2 / alph + tau^3 / 6 / alph^2 ...
            + tau / 2 - tau * elong / alph;
    elseif tau < elong
        p_term = elong^2 / 2 / alph - tau * elong / alph + tau^2 / 2 / alph;
    end
    
    %calculates dynamics terms
    d_term = 0;
    for i = 1:length(aes)
        a = aes(i);
        b = bes(i);
        cor = 0;
        if elong > alph

            %t > alph and t' > alph

            % t' + tau - t > 0
            % part 1
            if elong > alph + tau
                cor = cor + (elong - alph - tau) / b;
                cor = cor + 1 / b^2;
            else
                cor = cor + exp(-b * (alph + tau - elong)) / b^2;
            end
            if tau > alph
                if tau > elong
                    cor = cor + (exp(-b * tau) - exp(-b * (alph + tau - elong))) / b^2;
                else
                    cor = cor + (exp(-b * elong) - exp(-b * alph)) / b^2;
                end
            else
                cor = cor + (exp(-b * (elong + tau - alph)) - exp(-b * tau)) / b^2;
            end
            cor = cor - exp(-b * tau) / b^2;

            %part 2

            if tau > alph
                if elong > tau
                    cor = cor + (exp(-b* alph) - exp(-b * tau) + ...
                        exp(-b * (elong + tau - alph)) - exp(-b * elong)) / b^2;
                else
                    cor = cor + (exp(-b* (alph + tau - elong)) - exp(-b * tau) + ...
                        exp(-b * (elong + tau - alph)) - exp(-b * tau)) / b^2;
                end
            end

            %t' + tau - t < 0

            if alph + tau < elong
                cor = cor + (elong - alph - tau) / b;
                cor = cor + (exp(b * (alph + tau - elong)) - 1) / b^2;
            end

            %t > alph and t' < alph
            if alph > 0.1
                %t' + tau - t > 0


                piece1 = min(max(alph,tau),elong);
                piece2 = min(alph + tau,elong);

                if alph < tau
                    cor = cor + (exp(-b * (tau - piece1)) - exp(-b * (tau - alph))) ...
                    / (alph * b^3);
                end
                cor = cor + (exp(-b * tau) - exp(-b * (alph + tau - piece2))) ...
                    / b^2 + (exp(-b * tau) - ...
                    exp(-b * (alph + tau - piece2))) / b^3 / alph + (piece2 - piece1) / b^2 / alph + ...
                    (piece2^2 - piece1^2 + 2*tau * piece1 - 2*tau * piece2) / (2 * alph * b);

                % t' + tau - t < 0
                if elong > tau
                    cor = cor + (exp(b * (tau - piece1)) - exp(b * (tau - elong))) ...
                        / (alph * b^3) - (piece2 - piece1) / (alph * b^2) + ...
                        (piece2^2 - piece1^2) / (2 * alph * b) + tau * ...
                        (piece1 - piece2) / (alph * b);
                    if elong > tau + alph
                        cor = cor + (1/b^2 - 1 / alph / b^3) * ...
                            (exp(b * (alph + tau - piece2)) - exp(b * (alph + tau - elong)));
                    end
                end

                % t < alph and t' > alph

                % t' + tau - t > 0

                cor = cor + (exp(-b * tau) - exp(-b * (elong + tau - alph))) / b^2 ...
                    + (exp(-b * (elong + tau - alph)) + exp(-b * (alph + tau)) - ...
                    exp(-b * (elong + tau)) - exp(-b * tau)) / (alph * b^3);
            end
        end

        % t < alph and t' < alph

        if alph > 0.1 || alph >= elong
            alph = min(alph,elong);
            %t' + tau - t > 0

            piece3 = min(tau,alph);
            piece4 = max(tau - alph, 0);

            cor = cor + ((-1 / b^2) + (-1 / (b^3 * alph))) * exp(-b * tau) + ...
                ((1 / (b^3 * alph)) + (1 / (b^4 * alph^2))) * ...
                (exp(-b * tau) - exp(-b * (alph + tau))) + piece3 * ...
                exp(-b * piece4) / (b^3 * alph^2) + (exp(-b * tau) - ...
                exp(-b * piece4)) / (b^4 * alph^2) + 1 / (2 * b^2) - ...
                piece3^2 / (2 * b^2 * alph^2) + alph / (3 * b) - tau / (2 * b) ...
                + piece3^2 * tau / (2 * b * alph^2) - piece3^3 / (3 * b * alph^2);

            %t' + tau - t < 0
            if tau < alph
                cor = cor - exp(b * (tau - alph)) / b^3 / alph + piece3 * ...
                    exp(b * (tau - piece3)) / b^3 / alph^2 + (exp(b * (tau - piece3)) ...
                    - exp(b * (tau - alph))) / b^4 / alph^2 - 1 / 2 / b^2 + ...
                    piece3^2 / 2 / b^2 / alph^2 + alph / 3 / b - tau / 2 / b + ...
                    piece3^2 * tau / 2 / b / alph^2 - piece3^3 / 3 / b / alph^2;
            end

        end

        if (isnan(cor) || cor == inf || cor == -inf)
            debug = true;
        end
        d_term = d_term + a * cor;
    end
    cor_tot = p_term + d_term;
end

