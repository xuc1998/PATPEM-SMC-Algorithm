clc; clear; close all;

%% ============================================================
%  Example 1 for DREAM-Suite:
%  2D 20-mode Gaussian mixture target on [-1,10]^2
% =============================================================

%% User settings
DREAMPar.d        = 2;          % problem dimension
DREAMPar.lik      = 2;          % target returns log-likelihood / log-density
DREAMPar.thinning = 1;          % keep every sample (can change later)

% Initial sampling: use Latin hypercube in same box as PATPEMS
Par_info.initial       = 'latin';
Par_info.min           = [-1 -1];
Par_info.max           = [10 10];
Par_info.boundhandling = 'none';   % target itself handles outside-box points

% Target function name
Func_name = 'mix20_lik';

% DREAM method
method = 'mtdream_zs';
% method = 'dream_zs';
% method = 'dream';

switch lower(method)
    case {'dream','dream_d'}
        DREAMPar.N = 10;         % number of chains
        DREAMPar.T = 50000;      % number of generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 10;          % MT-DREAM(ZS) often uses fewer chains
        DREAMPar.T = 1000;      % can increase for stricter benchmark
end
DREAMPar.mt=5;
% Optional save / print
options.save  = 'no';
options.print = 'no';

%% Random seed
rng(1,'twister');

%% Reset likelihood/model-call counter
global MIX20_NFUN MIX20_NCALL
MIX20_NFUN  = 0;
MIX20_NCALL = 0;

%% Run DREAM-Suite
tic
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,Par_info,[],options);
toc

fprintf('\n===== Likelihood/model-call counter =====\n');
fprintf('Number of mix20_lik function entries = %d\n', MIX20_NFUN);
fprintf('Number of evaluated target points    = %d\n', MIX20_NCALL);

%% ------------------------------------------------------------
% Postprocess samples
% genparset is included in DREAM-Suite examples
% P columns are typically: parameters + diagnostics columns
%% ------------------------------------------------------------
P = extract_converged_samples(chain,DREAMPar,output);

fprintf('\nTotal retained posterior samples inside box: %d\n', size(P,1));

%% True mode locations
truth = [2.18 5.76;
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

%% ------------------------------------------------------------
% Simple scatter plot
%% ------------------------------------------------------------
figure('Color','w');
scatter(P(:,1), P(:,2), 8, 'b', 'filled'); hold on;
plot(truth(:,1), truth(:,2), 'rp', 'MarkerSize', 12, 'LineWidth', 1.5);
xlim([-1 10]); ylim([-1 10]);
xlabel('\theta_1'); ylabel('\theta_2');
legend('DREAM samples','True mode centers','Location','best');
title(sprintf('%s on 20-mode target', upper(method)));
axis square; box on;

%% ------------------------------------------------------------
% 2D histogram / density view
%% ------------------------------------------------------------
figure('Color','w');
histogram2(P(:,1), P(:,2), [80 80], 'DisplayStyle','tile', 'ShowEmptyBins','off');
hold on;
plot(truth(:,1), truth(:,2), 'rp', 'MarkerSize', 10, 'LineWidth', 1.2);
xlim([-1 10]); ylim([-1 10]);
xlabel('\theta_1'); ylabel('\theta_2');
title('Posterior sample density from DREAM');
axis square; box on; colorbar;

%% ------------------------------------------------------------
% Accuracy metrics: Ds and marginal W2
%% ------------------------------------------------------------

% ---- True mixture settings ----
Nm = size(truth,1);
w  = ones(Nm,1) / Nm;     % equal weights
comp_sd = 0.1;            % component std
comp_var = comp_sd^2;

% ---- True marginal mean and std of the 2D mixture ----
mu_ref = sum(truth .* w, 1);                       % 1 x d
var_between = sum(((truth - mu_ref).^2) .* w, 1); % 1 x d
var_ref = comp_var + var_between;                  % within + between
sigma_ref = sqrt(var_ref);                         % 1 x d

% ---- Estimated marginal mean and std from DREAM posterior samples ----
mu_est = mean(P, 1);
sigma_est = std(P, 0, 1);

% ---- Ds metric ----
Ds = Ds_metric(mu_ref, sigma_ref, mu_est, sigma_est);

% ---- Reference sample from the true mixture for marginal W2 ----
Mref = 100000;        % reference sample size for empirical W2
rng(2026, 'twister'); % fixed seed for reproducibility
idx = randi(Nm, Mref, 1);          % equal-weight mixture
Yref = truth(idx,:) + comp_sd * randn(Mref, DREAMPar.d);

% ---- Marginal W2 ----
[W2_per_dim, W2_mean] = marginal_w2(P, Yref, 512);

fprintf('\n===== Accuracy metrics after convergence =====\n');
fprintf('Ds = %.6f\n', Ds);
fprintf('W2_mean = %.6f\n', W2_mean);
fprintf('W2_per_dim = [%.6f, %.6f]\n', W2_per_dim(1), W2_per_dim(2));

%% ------------------------------------------------------------
% Draw 1D marginals – DREAM posterior vs true mixture
%% ------------------------------------------------------------
w = ones(20,1) / 20;    % mixture weights
samples_rep = P;        % convergence-retained DREAM posterior samples
% save Posterior_Samples_Run9.mat  samples_rep

fig_marg = plot_20mode_marginals(samples_rep, truth, 0.1, w, ...
    'KDEBandwidth', 0.08, ...
    'ShowSeedBand', false, ...
    'FigureTitle', sprintf('1D marginals — true vs DREAM posterior (%s)', upper(method)));

%% ------------------------------------------------------------
% Accuracy vs. number of model calls
%% ------------------------------------------------------------
lb = [-1 -1];
ub = [10 10];
DREAM_W2_result = plot_dream_w2_vs_calls(chain, output, DREAMPar, Yref, lb, ub, method, ...
    'RhatThreshold', 1.2, ...
    'W2Bins', 512, ...
    'MinSamples', 10, ...
    'SavePrefix', sprintf('%s_W2_vs_model_calls', upper(method)), ...
    'ExportFigure', true);