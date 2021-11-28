function [a_array, b_array] = decompose_matrix(R, lambdas)
% takes a transitation rate matrix *R* and a list of mean rates of intiation
% *lambdas* and finds what the values of a_i and b_i should be in the
    
    lambda_mat = diag(lambdas);
    [V,D] = eig(R);
    V_inv = inv(V);
    a_array = zeros(1, length(lambdas));
    b_array = zeros(1, length(lambdas));
    [~, idx0] = min(abs(diag(D)));
    p = V(:,idx0(1)) / sum(V(:,idx0(1)));
    for i = 1:length(lambdas)
        a_array(i) = lambdas * V(:,i) * V_inv(i,:) * lambda_mat * p;
        a_array(i) = abs(a_array(i) / sum(p .* lambdas'));
        b_array(i) = D(i,i);
    end
end

