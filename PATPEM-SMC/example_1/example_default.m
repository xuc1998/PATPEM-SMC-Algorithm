clc; clear; close all;

%% ============================================================
% Example 1: 2D 20-mode Gaussian mixture target
%
% Metrics:
%   - Ds: based on overall marginal mean and standard deviation
%   - W2: marginal W2 distance against reference samples drawn from
%         the true 20-mode Gaussian mixture
%
% Additional analysis:
%   - Accuracy vs. number of model calls for PATPEMS
% =============================================================

%% ============================================================
% 1. User-defined random seeds
% =============================================================
seeds = 1;
% seeds = [1, 7, 13, 23, 37, 41, 59, 83, 101, 211];

%% ============================================================
% 2. Common PATPEMS settings
% =============================================================
Np    = 800;                 % number of particles
S     = 500;                  % maximum number of beta stages
bound = [-1 -1;               % lower bounds
          10 10];             % upper bounds

d = size(bound, 2);

% Target log-density
logpdf = @(th) target(th);

% PATPEMS options
opts = struct();

% Example optional settings:
% opts.Parallel.Enabled = true;
% opts.Parallel.NumWorkers = 20;
% opts.Thinning = 1;

%% ============================================================
% 3. True reference settings for the 20-mode mixture
% =============================================================
Nm = 20;

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

w = ones(Nm, 1) / Nm;         % equal mixture weights
comp_sd  = 0.1;               % component standard deviation
comp_var = comp_sd^2;

% True marginal mean and std of the 2D mixture
mu_ref = sum(truth .* w, 1);                        % 1 x d
var_between = sum(((truth - mu_ref).^2) .* w, 1);  % 1 x d
var_ref = comp_var + var_between;                  % within + between
sigma_ref = sqrt(var_ref);                         % 1 x d

%% ============================================================
% 4. Add postprocessing path
% =============================================================
postproc_rel = fullfile('..', 'postprocessing');

if exist(postproc_rel, 'dir')
    addpath(postproc_rel);
end

%% ============================================================
% 5. Reference samples for W2
% =============================================================
Mref = 1e5;

rng(2026, 'twister');
idx_ref = randi(Nm, Mref, 1);
Yref = truth(idx_ref, :) + comp_sd * randn(Mref, d);

%% ============================================================
% 6. Storage
% =============================================================
nS = numel(seeds);

seedsCell     = cell(nS, 1);      % final posterior samples for each seed
paramIterCell = cell(nS, 1);      % stored particle trajectories
outCell       = cell(nS, 1);      % PATPEMS diagnostic outputs

Ds_vals      = nan(nS, 1);
W2_vals      = nan(nS, 1);
runtime_vals = nan(nS, 1);

%% ============================================================
% 7. Batch loop over seeds
% =============================================================
for i = 1:nS

    seed_i = seeds(i);
    fprintf('\n=== Running seed %d of %d: %d ===\n', i, nS, seed_i);

    rng(seed_i, 'twister');

    % Run PATPEMS
    tStart = tic;
    [parameter_iteration, out] = PATPEMS(Np, S, bound, logpdf, opts);
    runtime_vals(i) = toc(tStart);

    % Store full outputs for later plotting
    paramIterCell{i} = parameter_iteration;
    outCell{i}       = out;

    % Final-stage particles
    S_K = out.S_K;
    theta_final = squeeze(parameter_iteration(:, :, S_K));   % [Np x d]

    % Store final posterior samples
    seedsCell{i} = theta_final;

    % Ds metric
    mu_est = mean(theta_final, 1);
    sigma_est = std(theta_final, 0, 1);
    Ds_vals(i) = Ds_metric(mu_ref, sigma_ref, mu_est, sigma_est);

    % Marginal W2 metric
    [~, W2_vals(i)] = marginal_w2(theta_final, Yref, 512);

    fprintf('Seed %d finished: Runtime = %.3f s, Ds = %.6f, W2 = %.6f\n', ...
        seed_i, runtime_vals(i), Ds_vals(i), W2_vals(i));
end

%% ============================================================
% 8. Aggregate statistics across seeds
% =============================================================
Ds_mean_over_seeds = mean(Ds_vals, 'omitnan');
Ds_std_over_seeds  = std(Ds_vals,  'omitnan');

W2_mean_over_seeds = mean(W2_vals, 'omitnan');
W2_std_over_seeds  = std(W2_vals,  'omitnan');

runtime_mean_over_seeds = mean(runtime_vals, 'omitnan');
runtime_std_over_seeds  = std(runtime_vals,  'omitnan');

SummaryTable = table(seeds(:), runtime_vals, Ds_vals, W2_vals, ...
    'VariableNames', {'Seed', 'Runtime_s', 'Ds', 'W2_marginal_mean'});

fprintf('\n===== Summary over %d seed(s) =====\n', nS);
fprintf('Runtime:          mean = %.6g s, std = %.6g s\n', ...
    runtime_mean_over_seeds, runtime_std_over_seeds);
fprintf('Ds:               mean = %.6g, std = %.6g\n', ...
    Ds_mean_over_seeds, Ds_std_over_seeds);
fprintf('W2_marginal_mean: mean = %.6g, std = %.6g\n', ...
    W2_mean_over_seeds, W2_std_over_seeds);

disp(SummaryTable);

%% ============================================================
% 9. Representative run for plotting and saving
% =============================================================
idx_plot  = 1;
seed_plot = seeds(idx_plot);

samples_rep             = seedsCell{idx_plot};
parameter_iteration_rep = paramIterCell{idx_plot};
out_rep                 = outCell{idx_plot};

fprintf('\n>>> Representative run for plotting: seed = %d <<<\n', seed_plot);

% Save representative final posterior samples
% save('Posterior_Samples_PATPEMS_Run5.mat', 'samples_rep');

%% ============================================================
% 10. Accuracy vs. number of model calls
% =============================================================
plot_call_accuracy = true;

if plot_call_accuracy

    PATPEMS_W2_result = plot_patpems_w2_vs_calls( ...
        parameter_iteration_rep, out_rep, Yref, ...
        'Bound', bound, ...
        'W2Bins', 512, ...
        'MinSamples', 10, ...
        'BetaMin', [], ...
        'SavePrefix', 'PATPEMS_W2_vs_model_calls', ...
        'ExportFigure', true, ...
        'FigureTitle', 'Accuracy vs. model calls — PATPEMS');
%     exportgraphics(gcf, 'PATPEMS_ZS_W2_vs_model_calls.tif', 'Resolution', 600);
else
    PATPEMS_W2_result = [];
end

%% ============================================================
% 11. Optional PATPEMS diagnostic figures
% =============================================================
plot_diagnostics = true;

if plot_diagnostics

    [figs, metrics_plot] = plot_pem_smc_all( ...
        parameter_iteration_rep, out_rep, bound, opts, ...
        'NumTraj', 30, ...
        'Seed', seed_plot, ...
        'Truth', truth, ...
        'SigmaRef', ones(Nm, d) * comp_var, ...
        'W2Method', 'bimodal-ref', ...
        'W2RefM', 1e5, ...
        'Plots', {'trajectories', 'diagnostics', 'posteriors', ...
                  'Distance_eval', 'scatter2d'});

    % exportgraphics(gcf, 'Diagnostics_PATPEMS.tif', 'Resolution', 600);

else
    figs = [];
    metrics_plot = [];
end

%% ============================================================
% 12. 1D marginal posterior plots
% =============================================================
plot_marginals = true;

if plot_marginals

    % Single representative seed
    fig1 = plot_20mode_marginals(samples_rep, truth, comp_sd, w, ...
        'KDEBandwidth', 0.08, ...
        'ShowSeedBand', false, ...
        'FigureTitle', sprintf( ...
            '1D marginals — truth vs PATPEMS posterior (seed = %d)', seed_plot));

    % Multiple-seed posterior mean ± SD, only when repeated seeds are used
    if nS > 1
        fig2 = plot_20mode_marginals(samples_rep, truth, comp_sd, w, ...
            'SeedSamples', seedsCell, ...
            'KDEBandwidth', 0.07, ...
            'FigureTitle', ...
            '1D marginals — truth vs PATPEMS posterior (repeated seeds: mean ± SD)');

        % exportgraphics(gcf, 'Multi_Seeds_Results.tif', 'Resolution', 600);
    else
        fig2 = [];
    end

else
    fig1 = [];
    fig2 = [];
end

%% ============================================================
% 13. Save main results
% =============================================================
% save('case1_patpems_results.mat', ...
%      'seeds', 'seed_plot', ...
%      'Np', 'S', 'bound', 'opts', ...
%      'truth', 'w', 'comp_sd', 'comp_var', ...
%      'mu_ref', 'sigma_ref', 'Yref', ...
%      'samples_rep', 'seedsCell', ...
%      'parameter_iteration_rep', 'out_rep', ...
%      'Ds_vals', 'W2_vals', 'runtime_vals', 'SummaryTable', ...
%      'PATPEMS_W2_result');
% 
% fprintf('\nAll PATPEMS case-1 results have been saved.\n');