function fig = plot_marginal_pdf_band_100d(theta, varargin)
% Plot 100D marginal posterior PDFs as an envelope band + median curve,
% overlay the 1D analytic reference mixture density, and annotate W2/Ds.
%
% theta: [N x d] posterior samples at beta=1 (e.g., resampled + moved particles)
%
% Name-Value:
%   'XLim'      [-10 10]
%   'Nx'        number of grid points (default 600)
%   'BandPct'   [5 95] percentile band for PDF envelope
%   'ShowDims'  true/false, if true plots faint lines for each dim (default false)
%   'Seed'      RNG seed for W2 reference sampling (default 1)

p = inputParser;
p.addParameter('XLim', [-10 10]);
p.addParameter('Nx', 600);
p.addParameter('BandPct', [5 95]);
p.addParameter('ShowDims', false);
p.addParameter('Seed', 1);
p.parse(varargin{:});
opt = p.Results;

[N, d] = size(theta);
x = linspace(opt.XLim(1), opt.XLim(2), opt.Nx);

% -------- Reference 1D density: (1/3)N(-5,1) + (2/3)N(5,1) --------
ref_pdf = (1/3)*normpdf(x, -5, 1) + (2/3)*normpdf(x, +5, 1);

% -------- Reference CDF for Ds (KS distance) --------
ref_cdf = @(t) (1/3)*normcdf(t, -5, 1) + (2/3)*normcdf(t, +5, 1);

% Storage
pdf_mat = zeros(d, opt.Nx);
Ds = zeros(d, 1);
W2 = zeros(d, 1);

% For W2: draw a reference sample (same N) per dim (cheap and stable here)
rng(opt.Seed);
% Generate reference samples via mixture
% (vectorized: generate N*d samples in one shot then reshape)
uMix = rand(N*d,1);
z = randn(N*d,1);
ref_samp_all = (-5 + z).*(uMix <= 1/3) + ( 5 + z).*(uMix >  1/3);
ref_samp_all = reshape(ref_samp_all, [N, d]);

for j = 1:d
    sj = theta(:,j);

    % --- KDE on common grid (PDF) ---
    % Use ksdensity for smooth curve; it returns a PDF estimate
    pdf_mat(j,:) = ksdensity(sj, x, 'Function','pdf');

    % --- Ds: KS distance between empirical CDF and reference CDF ---
    [Femp, Xemp] = ecdf(sj);
    Fref = ref_cdf(Xemp);
    Ds(j) = max(abs(Femp - Fref));

    % --- W2 (approx): 1D Wasserstein-2 between samples and ref-samples ---
    % In 1D, W2^2 = mean((sorted_s - sorted_r).^2)
    sr = ref_samp_all(:,j);
    sj_sort = sort(sj);
    sr_sort = sort(sr);
    W2(j) = sqrt(mean((sj_sort - sr_sort).^2));
end

% -------- Build band + median across dims at each x --------
lo = prctile(pdf_mat, opt.BandPct(1), 1);   % 1 x Nx
hi = prctile(pdf_mat, opt.BandPct(2), 1);
med = prctile(pdf_mat, 50, 1);

% -------- Plot --------
fig = figure('Color','w','Name','100D marginal PDF band');
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');

% Optional: show each-dim pdf faintly (can clutter; default off)
if opt.ShowDims
    for j = 1:d
        plot(ax, x, pdf_mat(j,:), 'LineWidth', 0.5);
    end
end

% Band as filled area
fill(ax, [x, fliplr(x)], [lo, fliplr(hi)], ...
    [0.7 0.7 0.7], 'FaceAlpha', 0.85, 'EdgeColor', 'none');

% Median curve
plot(ax, x, med, 'k-', 'LineWidth', 1.8);

% Reference curve (red)
plot(ax, x, ref_pdf, 'r-', 'LineWidth', 2.0);

% Labels
xlabel(ax, '参数值 \theta');
ylabel(ax, '概率密度（PDF）');
xlim(ax, opt.XLim);

% Annotate metrics (mean ± std)
txt = sprintf(['W_2 = %.4f \\pm %.4f\n' ...
               'D_s = %.4f \\pm %.4f\n'], ...
    mean(W2), std(W2), mean(Ds), std(Ds));

text(ax, 0.02, 0.7, txt, 'Units','normalized', ...
    'VerticalAlignment','top', 'FontSize', 11, 'BackgroundColor','w', 'Margin', 6);

legend(ax, {'PDF95%区间带','PDF中位数','参考密度'}, 'Location','northeast');
end
