function result = plot_dream_w2_vs_calls(chain, output, DREAMPar, Yref, lb, ub, method, varargin)
%PLOT_DREAM_W2_VS_CALLS
% Plot posterior approximation accuracy vs. number of model calls for DREAM.
%
% The convergence start is determined by the sustained multivariate R-hat
% criterion: the earliest point after which all subsequent R-hat values remain
% below the specified threshold.
%
% INPUTS
%   chain    : DREAM output chain
%   output   : DREAM output structure containing MR_stat or MR
%   DREAMPar : DREAM parameter structure
%   Yref     : reference samples from true/reference posterior
%   method   : DREAM method name, e.g., 'dream_zs'
%
% OPTIONAL NAME-VALUE PAIRS
%   'RhatThreshold' : threshold for multivariate R-hat, default = 1.2
%   'W2Bins'        : number of bins/quantiles used in marginal_w2, default = 512
%   'MinSamples'    : minimum accumulated posterior samples for W2, default = 10
%   'SavePrefix'    : prefix for saved .mat and .png files
%   'ExportFigure'  : true/false, whether to export figure, default = true
%
% OUTPUT
%   result : structure containing model calls, W2 values, convergence point, etc.

%% ------------------------------------------------------------
% Parse optional inputs
%% ------------------------------------------------------------
p = inputParser;
p.addParameter('RhatThreshold', 1.2, @(x) isnumeric(x) && isscalar(x));
p.addParameter('W2Bins', 512, @(x) isnumeric(x) && isscalar(x));
p.addParameter('MinSamples', 10, @(x) isnumeric(x) && isscalar(x));
p.addParameter('SavePrefix', 'DREAM_W2_vs_model_calls', @(x) ischar(x) || isstring(x));
p.addParameter('ExportFigure', true, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});

Rhat_thr    = p.Results.RhatThreshold;
W2Bins      = p.Results.W2Bins;
MinSamples  = p.Results.MinSamples;
SavePrefix  = char(p.Results.SavePrefix);
ExportFig   = logical(p.Results.ExportFigure);

%% ------------------------------------------------------------
% Reconstruct full DREAM sample set before burn-in truncation
%% ------------------------------------------------------------
P_all = genparset(chain);
P_all = P_all(:, 1:DREAMPar.d);

nTotal = size(P_all, 1);

%% ------------------------------------------------------------
% Read multivariate R-hat diagnostic
%% ------------------------------------------------------------
if isfield(output, 'MR_stat')
    MR = output.MR_stat;
elseif isfield(output, 'MR')
    MR = output.MR;
else
    error('Cannot find output.MR_stat or output.MR in DREAM output.');
end

model_calls_MR = MR(:, 1);
Rhat_MR        = MR(:, 2);

valid = isfinite(model_calls_MR) & isfinite(Rhat_MR);
model_calls_MR = model_calls_MR(valid);
Rhat_MR        = Rhat_MR(valid);

% Keep only evaluation points within available sample range
keep = model_calls_MR <= nTotal;
model_calls_MR = model_calls_MR(keep);
Rhat_MR        = Rhat_MR(keep);

%% ------------------------------------------------------------
% Find sustained convergence point
% Not the first Rhat < threshold, but the first point after which
% all subsequent Rhat values remain below the threshold.
%% ------------------------------------------------------------
Rhat_tmp = Rhat_MR;
Rhat_tmp(~isfinite(Rhat_tmp)) = Inf;

futureMaxR = flipud(cummax(flipud(Rhat_tmp)));
idx_conv = find(futureMaxR < Rhat_thr, 1, 'first');

if isempty(idx_conv)
    warning('No sustained convergence point with multivariate R-hat < %.2f was found.', Rhat_thr);

    result = struct();
    result.success = false;
    result.message = sprintf('No sustained convergence point with Rhat < %.2f.', Rhat_thr);
    return
end

conv_call = round(model_calls_MR(idx_conv));
conv_Rhat = Rhat_MR(idx_conv);

fprintf('\n===== DREAM sustained convergence diagnostic =====\n');
fprintf('Sustained convergence starts at %d model calls.\n', conv_call);
fprintf('Multivariate R-hat at this point = %.4f\n', conv_Rhat);

%% ------------------------------------------------------------
% Evaluate W2 using progressively accumulated post-convergence samples
%% ------------------------------------------------------------
eval_calls = round(model_calls_MR(idx_conv+1:end));
eval_calls = eval_calls(eval_calls > conv_call & eval_calls <= nTotal);
eval_calls = unique(eval_calls, 'stable');

W2_curve     = nan(length(eval_calls), 1);
W2_dim_curve = nan(length(eval_calls), DREAMPar.d);
nPost_curve  = nan(length(eval_calls), 1);

for k = 1:length(eval_calls)

    c = eval_calls(k);

    % Accumulated post-convergence samples
    P_current = P_all(conv_call+1:c, :);

    % Remove outside-box samples, consistent with the main analysis
    inside_current = all(P_current >= lb, 2) & all(P_current <= ub, 2);
    P_current = P_current(inside_current, :);

    nPost_curve(k) = size(P_current, 1);

    if size(P_current, 1) >= MinSamples
        [W2_per_dim_k, W2_mean_k] = marginal_w2(P_current, Yref, W2Bins);

        W2_curve(k) = W2_mean_k;
        W2_dim_curve(k, :) = W2_per_dim_k;
    end
end

ok = isfinite(W2_curve);

%% ------------------------------------------------------------
% Plot
%% ------------------------------------------------------------
fig = figure('Color', 'w');

plot(eval_calls(ok), W2_curve(ok), '-o', ...
    'LineWidth', 1.4, ...
    'MarkerSize', 5);
hold on;

% Mark the sustained convergence point
xline(conv_call, '--', 'LineWidth', 1.0);

% Add a clearer text annotation
xl = xlim;
yl = ylim;

text(conv_call + 0.015*diff(xl), yl(1) + 0.04*diff(yl), ...
    'Start of post-convergence samples', ...
    'FontSize', 11, ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'left');

xlabel('Number of model calls');
ylabel('Marginal $W_2$ distance', 'Interpreter', 'latex');
title(sprintf('Accuracy vs. model calls — %s', upper(method)));
box on;
grid on;

%% ------------------------------------------------------------
% Save results
%% ------------------------------------------------------------
result = struct();
result.success        = true;
result.method         = method;
result.model_calls    = eval_calls(ok);
result.W2_mean        = W2_curve(ok);
result.W2_per_dim     = W2_dim_curve(ok, :);
result.nPostSamples   = nPost_curve(ok);
result.conv_call      = conv_call;
result.conv_Rhat      = conv_Rhat;
result.Rhat_threshold = Rhat_thr;
result.MR_stat        = [model_calls_MR, Rhat_MR];

save([SavePrefix, '.mat'], 'result');

if ExportFig
    exportgraphics(fig, [SavePrefix, '.png'], 'Resolution', 600);
end

fprintf('\n===== W2 vs. model calls analysis finished =====\n');
fprintf('Convergence start: %d model calls\n', result.conv_call);
fprintf('Final evaluation point: %d model calls\n', result.model_calls(end));
fprintf('Final W2_mean: %.6f\n', result.W2_mean(end));

end
