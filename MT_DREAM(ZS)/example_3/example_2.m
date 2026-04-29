clc; clear; close all;

%% ============================================================
%  Example 2 for DREAM-Suite:
%  100D bimodal Gaussian mixture target on [-10,10]^100
%  pi(x) = (1/3) N(-5*1_d, I_d) + (2/3) N(5*1_d, I_d)
% =============================================================

%% User settings
DREAMPar.d        = 30;        % problem dimension
DREAMPar.lik      = 2;          % target returns log-likelihood / log-density
DREAMPar.thinning = 1;          % keep every sample

% Initial sampling: use Latin hypercube in same box as PATPEMS
Par_info.initial       = 'latin';
Par_info.min           = -10 * ones(1, DREAMPar.d);
Par_info.max           =  10 * ones(1, DREAMPar.d);
Par_info.boundhandling = 'none';   % same style as Example 1

% Target function name
Func_name = 'bimodal100_lik';

% DREAM method
method = 'mtdream_zs';
% method = 'dream_zs';
% method = 'dream';

switch lower(method)
    case {'dream','dream_d'}
        DREAMPar.N = 10;        % number of chains
        DREAMPar.T = 50000;     % number of generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 10;        % same style as Example 1
        DREAMPar.T = 100000;     % can increase later for stricter benchmark
end
DREAMPar.mt=10;       % Number of multi-try proposals
% Optional save / print
options.save  = 'no';
options.print = 'no';

%% Random seed
rng(1,'twister');

%% Run DREAM-Suite
tic
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,Par_info,[],options);
toc

%% ------------------------------------------------------------
% Postprocess samples
%% ------------------------------------------------------------
P = extract_converged_samples(chain,DREAMPar,output);
fprintf('\nTotal retained posterior samples inside box: %d\n', size(P,1));

%% ------------------------------------------------------------
% True mixture settings
%% ------------------------------------------------------------
truth = [-5 * ones(1, DREAMPar.d);
          5 * ones(1, DREAMPar.d)];

w1 = 1/3;
w2 = 2/3;

comp_sd  = 1;         % component std
comp_var = comp_sd^2; % = 1

%% ------------------------------------------------------------
% Accuracy metrics: Ds and marginal W2
%% ------------------------------------------------------------

% ---- True marginal mean and std of the 100D mixture ----
mu_ref = w1 * truth(1,:) + w2 * truth(2,:);
var_between = w1 * (truth(1,:) - mu_ref).^2 + ...
              w2 * (truth(2,:) - mu_ref).^2;
var_ref   = comp_var + var_between;
sigma_ref = sqrt(var_ref);

% ---- Estimated marginal mean and std from DREAM posterior samples ----
mu_est    = mean(P, 1);
sigma_est = std(P, 0, 1);

% ---- Ds metric ----
Ds = Ds_metric(mu_ref, sigma_ref, mu_est, sigma_est);

% ---- Reference sample from the true mixture for marginal W2 ----
Mref = 100000;
rng(2026, 'twister');

u   = rand(Mref,1);
idx = 2 * ones(Mref,1);
idx(u <= w1) = 1;

Yref = truth(idx,:) + comp_sd * randn(Mref, DREAMPar.d);

% ---- Marginal W2 ----
[W2_per_dim, W2_mean] = marginal_w2(P, Yref, 512);

fprintf('\n===== Accuracy metrics after convergence =====\n');
fprintf('Ds = %.6f\n', Ds);
fprintf('W2_mean = %.6f\n', W2_mean);
fprintf('Mean of W2_per_dim = %.6f\n', mean(W2_per_dim));

%% ------------------------------------------------------------
% Save posterior samples
%% ------------------------------------------------------------
samples_rep = P;
save Posterior_Samples_Run8.mat samples_rep

%% ------------------------------------------------------------
% Posterior comparison plots for Example 2
%% ------------------------------------------------------------

% Common x-grid and true 1D marginal
xgrid = linspace(-10, 10, 600);
true_pdf = (1/3) * normpdf(xgrid, -5, 1) + (2/3) * normpdf(xgrid,  5, 1);

% -------- Figure 1: randomly choose 20 dimensions (fixed seed) --------
rng(2026, 'twister');                 % fixed seed for reproducibility
dims20 = sort(randperm(DREAMPar.d,20));

fprintf('\nRandom 20 dimensions used in Figure 1:\n');
disp(dims20);

fig20 = plot_bimodal_random20_marginals(samples_rep, dims20, xgrid, true_pdf, ...
    'NBins', 35, ...
    'FigureTitle', sprintf('Random 20 marginals — MT-DREAM(ZS) posterior vs truth'));

% -------- Figure 2: overall 100D marginal band --------
figBand = plot_bimodal_100d_band(samples_rep, xgrid, true_pdf, ...
    'Bandwidth', 0.30, ...
    'BandPct', [2.5 97.5], ...
    'CenterStat', 'median', ...
    'FigureTitle', sprintf('100D marginal density band — DREAM posterior vs truth (%s)', upper(method)));

% Optional save
% exportgraphics(fig20,  'Posterior_20dims_DREAM_zs.tif');
% exportgraphics(figBand,'Posterior_100d_band_DREAM_zs.tif');