function loglik = mix20_lik(x)
%MIX20_LIK Log target density for a 2D, 20-mode Gaussian mixture
%
% INPUT
%   x : N x 2 matrix, each row is one point in R^2
%
% OUTPUT
%   loglik : 1 x N vector, log density for each row of x

% ============================================================
% Count likelihood/model calls
% MIX20_NFUN  : number of times this function is entered
% MIX20_NCALL : number of evaluated points, i.e., model/likelihood calls
% ============================================================
global MIX20_NFUN MIX20_NCALL

if isempty(MIX20_NFUN)
    MIX20_NFUN = 0;
end
if isempty(MIX20_NCALL)
    MIX20_NCALL = 0;
end

if isvector(x)
    x = reshape(x, 1, []);
end

N = size(x,1);

MIX20_NFUN  = MIX20_NFUN + 1;
MIX20_NCALL = MIX20_NCALL + N;

% ============================================================
% Original likelihood calculation
% ============================================================
Nm    = 20;
sigma = ones(1, Nm) * (0.1^2);

mu = [2.18 5.76;
      8.67 9.59;
      4.24 8.48;
      8.41 1.68;
      3.93 7.82;
      3.25 3.47;
      1.70 0.50;
      4.59 5.60;
      6.91 5.81;
      6.87 5.40;
      5.41 2.65;
      2.70 7.88;
      4.98 3.70;
      1.14 2.39;
      8.33 9.50;
      4.93 1.50;
      1.83 0.09;
      2.26 0.31;
      5.54 6.86;
      1.69 8.11];

w  = ones(1, Nm) * 0.05;

loglik = zeros(1,N);

for n = 1:N
    f = 0;
    xn = x(n,:);
    for i = 1:Nm
        diff = xn - mu(i,:);
        f = f + w(i) / (2*pi*sigma(i)) * exp( - (diff*diff') / (2*sigma(i)) );
    end
    loglik(n) = log(f);
end

end