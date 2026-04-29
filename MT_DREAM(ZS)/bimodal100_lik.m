function loglik = bimodal100_lik(x)
%BIMODAL100_LIK Log target density for a d-D bimodal Gaussian mixture
%
% INPUT
%   x : N x d matrix, each row is one point in R^d
%
% OUTPUT
%   loglik : 1 x N vector, log density for each row of x
%
% Target:
%   pi(x) = (1/3) N(-5*1_d, I_d) + (2/3) N(5*1_d, I_d)

    if isvector(x)
        x = reshape(x, 1, []);
    end

    [N, d] = size(x);

    % mixture weights and means
    w1  = 1/3;
    w2  = 2/3;
    mu1 = -5 * ones(1, d);
    mu2 =  5 * ones(1, d);

    % log-density constant for N(mu, I_d)
    c0 = -0.5 * d * log(2*pi);

    % log component densities
    diff1 = x - mu1;
    l1 = log(w1) + c0 - 0.5 * sum(diff1.^2, 2);

    diff2 = x - mu2;
    l2 = log(w2) + c0 - 0.5 * sum(diff2.^2, 2);

    % log-sum-exp for numerical stability
    m = max([l1, l2], [], 2);
    loglik_col = m + log(exp(l1 - m) + exp(l2 - m));

    % DREAM-Suite examples usually return 1 x N
    loglik = loglik_col.';
end