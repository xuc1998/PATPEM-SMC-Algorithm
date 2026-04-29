function result = plot_patpems_w2_vs_calls(parameter_iteration, out, Yref, varargin)
%PLOT_PATPEMS_W2_VS_CALLS
% Plot posterior approximation accuracy vs. number of model calls for PATPEMS.
%
% For PATPEMS, each stored stage contains a fixed-size particle population.
% The W2 distance is therefore computed from the particle population at each
% stored stage, while the x-axis uses the cumulative number of logpdf/model
% evaluations recorded in out.calls.
%
% INPUTS
%   parameter_iteration : [Np x d x K] stored particle populations
%   out                 : PATPEMS output struct with out.calls and out.beta
%   Yref                : reference samples from the true/reference posterior
%
% OPTIONAL NAME-VALUE PAIRS
%   'Bound'        : 2 x d array, [LB; UB], used to remove outside-box samples
%   'W2Bins'       : number of bins/quantiles used in marginal_w2, default = 512
%   'MinSamples'   : minimum number of particles required for W2, default = 10
%   'BetaMin'      : only plot stages with beta >= BetaMin; default = []
%   'SavePrefix'   : prefix for saved files, default = 'PATPEMS_W2_vs_model_calls'
%   'ExportFigure' : true/false, default = true
%   'FigureTitle'  : figure title
%
% OUTPUT
%   result : struct with model calls, beta, W2 values, and sample counts

%% ------------------------------------------------------------
% Parse inputs
%% ------------------------------------------------------------
p = inputParser;
p.addParameter('Bound', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('W2Bins', 512, @(x) isnumeric(x) && isscalar(x));
p.addParameter('MinSamples', 10, @(x) isnumeric(x) && isscalar(x));
p.addParameter('BetaMin', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('SavePrefix', 'PATPEMS_W2_vs_model_calls', @(x) ischar(x) || isstring(x));
p.addParameter('ExportFigure', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('FigureTitle', 'Accuracy vs. model calls — PATPEMS', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

bound      = p.Results.Bound;
W2Bins     = p.Results.W2Bins;
MinSamples = p.Results.MinSamples;
BetaMin    = p.Results.BetaMin;
SavePrefix = char(p.Results.SavePrefix);
ExportFig  = logical(p.Results.ExportFigure);
FigureTitle = char(p.Results.FigureTitle);

%% ------------------------------------------------------------
% Basic checks
%% ------------------------------------------------------------
if ~isfield(out, 'calls')
    error('out.calls is missing. Please run the updated PATPEMS with call counting enabled.');
end

if ~isfield(out, 'beta')
    error('out.beta is missing.');
end

K_particles = size(parameter_iteration, 3);
K_out = min(numel(out.calls), numel(out.beta));
K = min(K_particles, K_out);

calls = out.calls(1:K);
beta  = out.beta(1:K);

% Remove invalid stages
valid_stage = isfinite(calls) & isfinite(beta);

% Optional: plot only stages after a beta threshold
if ~isempty(BetaMin)
    valid_stage = valid_stage & (beta >= BetaMin);
end

idx_stage = find(valid_stage);

if isempty(idx_stage)
    error('No valid stored stages found for W2-vs-calls analysis.');
end

%% ------------------------------------------------------------
% Compute W2 at each stored stage
%% ------------------------------------------------------------
d = size(parameter_iteration, 2);
W2_mean = nan(numel(idx_stage), 1);
W2_dim  = nan(numel(idx_stage), d);
nUsed   = nan(numel(idx_stage), 1);

for kk = 1:numel(idx_stage)

    k = idx_stage(kk);

    P_current = squeeze(parameter_iteration(:,:,k));

    % Remove NaN rows
    P_current = P_current(all(isfinite(P_current), 2), :);

    % Optional: remove outside-box particles
    if ~isempty(bound)
        lb = bound(1,:);
        ub = bound(2,:);
        inside = all(P_current >= lb, 2) & all(P_current <= ub, 2);
        P_current = P_current(inside, :);
    end

    nUsed(kk) = size(P_current, 1);

    if nUsed(kk) >= MinSamples
        [W2_per_dim_k, W2_mean_k] = marginal_w2(P_current, Yref, W2Bins);
        W2_dim(kk,:) = W2_per_dim_k;
        W2_mean(kk)  = W2_mean_k;
    end
end

ok = isfinite(W2_mean);

%% ------------------------------------------------------------
% Plot
%% ------------------------------------------------------------
fig = figure('Color', 'w');

plot(calls(idx_stage(ok)), W2_mean(ok), '-o', ...
    'LineWidth', 1.4, ...
    'MarkerSize', 5);

xlabel('Number of model calls');
ylabel('Marginal $W_2$ distance', 'Interpreter', 'latex');
title(FigureTitle);
box on;
grid on;

% Optional annotation for final beta = 1 point
hold on;
idx_final = find(beta(idx_stage(ok)) >= 1 - 1e-12, 1, 'first');
if ~isempty(idx_final)
    x_final = calls(idx_stage(ok));
    y_final = W2_mean(ok);
    plot(x_final(idx_final), y_final(idx_final), 'o', ...
        'MarkerSize', 7, ...
        'LineWidth', 1.5);
end

%% ------------------------------------------------------------
% Save result
%% ------------------------------------------------------------
result = struct();
result.success      = true;
result.model_calls  = calls(idx_stage(ok));
result.beta         = beta(idx_stage(ok));
result.W2_mean      = W2_mean(ok);
result.W2_per_dim   = W2_dim(ok,:);
result.nParticles   = nUsed(ok);
result.W2Bins       = W2Bins;
result.BetaMin      = BetaMin;

save([SavePrefix, '.mat'], 'result');

if ExportFig
    exportgraphics(fig, [SavePrefix, '.png'], 'Resolution', 600);
end

fprintf('\n===== PATPEMS W2 vs. model calls analysis finished =====\n');
fprintf('First plotted call: %d, beta = %.6f, W2 = %.6f\n', ...
    result.model_calls(1), result.beta(1), result.W2_mean(1));
fprintf('Final plotted call: %d, beta = %.6f, W2 = %.6f\n', ...
    result.model_calls(end), result.beta(end), result.W2_mean(end));

end