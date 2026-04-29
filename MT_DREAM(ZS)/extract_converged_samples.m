function [P] = extract_converged_samples(chain,DREAMPar,output)
%% ------------------------------------------------------------
% Postprocess samples
% Automatically keep post-convergence samples based on multivariate R-hat
%% ------------------------------------------------------------

P_all = genparset(chain);
P_all = P_all(:, 1:DREAMPar.d);

% Read multivariate R-hat diagnostic
if isfield(output, 'MR_stat')
    MR = output.MR_stat;
elseif isfield(output, 'MR')
    MR = output.MR;
else
    error('Cannot find output.MR_stat or output.MR in DREAM output.');
end

model_calls_MR = MR(:,1);   % cumulative model calls / sample numbers
Rhat_MR        = MR(:,2);   % multivariate R-hat

% Remove invalid rows
valid = isfinite(model_calls_MR) & isfinite(Rhat_MR);
model_calls_MR = model_calls_MR(valid);
Rhat_MR        = Rhat_MR(valid);

% Convergence threshold
Rhat_thr = 1.2;

% Find the earliest point after which all subsequent Rhat values remain < 1.2
% This avoids treating a temporary crossing of the threshold as convergence.
futureMaxR = flipud(cummax(flipud(Rhat_MR)));
idx_conv = find(futureMaxR < Rhat_thr, 1, 'first');

if isempty(idx_conv)
    error('No sustained convergence point with multivariate R-hat < %.2f was found.', Rhat_thr);
end

conv_call = round(model_calls_MR(idx_conv));
conv_Rhat = Rhat_MR(idx_conv);

fprintf('\n===== DREAM sustained convergence diagnostic =====\n');
fprintf('Sustained convergence starts after %d model calls.\n', conv_call);
fprintf('Multivariate R-hat at this point = %.4f\n', conv_Rhat);

% Keep only post-convergence samples
% Because conv_call corresponds to the cumulative sample/model-call count,
% samples after this point are retained.
P = P_all(conv_call+1:end, :);
end