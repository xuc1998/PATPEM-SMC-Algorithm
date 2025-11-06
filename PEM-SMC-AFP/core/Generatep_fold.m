function theta_new = Generatep_fold(particles, k, LB, UB, DEMH)
%GENERATEP_FOLD  DE–MH proposal with folding back into the box.
%
% Input
%   particles : N x d matrix of current particles
%   k         : index of the particle to mutate (1..N)
%   LB, UB    : 1 x d lower/upper bounds
%   DEMH      : struct with fields:
%                 .Gamma    : DE jump rate (default 2.38/sqrt(2*d))
%                 .NoiseSD  : Gaussian jitter std (default 1e-4)
%                 .FoldType : 'fold' or 'reflect'
%
% Output
%   theta_new : 1 x d DE proposal after folding
%
% Notes
% - Parents (r1,r2) are drawn as an ordered pair uniformly (r1~=r2, r1/r2~=k).
% - With symmetric jitter and ordered-pair symmetry, the folded proposal is
%   symmetric; MH acceptance uses standard ratio with pi_s.

[N, d] = size(particles);

if ~isfield(DEMH,'Gamma'),    DEMH.Gamma   = 2.38 / sqrt(2*d); end
if ~isfield(DEMH,'NoiseSD'),  DEMH.NoiseSD = 1e-4; end
if ~isfield(DEMH,'FoldType'), DEMH.FoldType = 'fold'; end

% choose two distinct parents excluding k
B = setdiff(1:N, k);
R = B(randperm(N-1, 2));

prop = particles(k,:) ...
     + DEMH.Gamma .* (particles(R(1),:) - particles(R(2),:)) ...
     + DEMH.NoiseSD .* randn(1,d);

switch lower(DEMH.FoldType)
    case 'fold'
        theta_new = wrap_into_box(prop, LB, UB);
    case 'reflect'
        theta_new = reflect_into_box(prop, LB, UB);
    otherwise
        error('DEMH.FoldType must be ''fold'' or ''reflect''.');
end
end

% --------- local folding utilities ---------
function y = wrap_into_box(z, LB, UB)
L = UB - LB;
y = LB + mod(z - LB, L);
end

function y = reflect_into_box(z, LB, UB)
L = UB - LB;
y = (z - LB) ./ L;
y = abs(mod(y,2) - 1);
y = LB + y .* L;
end
