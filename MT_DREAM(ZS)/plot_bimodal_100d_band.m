function fig = plot_bimodal_100d_band(samples_rep, xgrid, true_pdf, varargin)
% Plot the overall band of 100 marginal posterior density estimates
% and compare with the true 1D marginal.
%
% INPUT
%   samples_rep : Ns x d posterior sample matrix
%   xgrid       : 1 x Nx grid
%   true_pdf    : 1 x Nx true 1D marginal density
%
% OPTIONAL
%   'Bandwidth'   : KDE bandwidth (default 0.30)
%   'BandPct'     : percentile band, e.g. [2.5 97.5]
%   'CenterStat'  : 'median' or 'mean'
%   'FigureTitle' : figure title

    p = inputParser;
    addParameter(p, 'Bandwidth', 0.30);
    addParameter(p, 'BandPct', [2.5 97.5]);
    addParameter(p, 'CenterStat', 'median');
    addParameter(p, 'FigureTitle', '100D marginal density band');
    parse(p, varargin{:});

    bw         = p.Results.Bandwidth;
    bandPct    = p.Results.BandPct;
    centerStat = lower(p.Results.CenterStat);
    figTitle   = p.Results.FigureTitle;

    [~, d] = size(samples_rep);
    nx = numel(xgrid);

    % KDE for each dimension on the same grid
    fhat_all = zeros(d, nx);
    for j = 1:d
        fhat_all(j,:) = ksdensity(samples_rep(:,j), xgrid, ...
            'Bandwidth', bw);
    end

    % Pointwise band across 100 dimensions
    f_lo = prctile(fhat_all, bandPct(1), 1);
    f_hi = prctile(fhat_all, bandPct(2), 1);

    switch centerStat
        case 'mean'
            f_center = mean(fhat_all, 1);
            centerName = 'Mean estimated density';
        otherwise
            f_center = median(fhat_all, 1);
            centerName = 'Median estimated density';
    end

    fig = figure('Color','w', 'Position',[180 120 900 520]);
    hold on;

    % band
    fill([xgrid, fliplr(xgrid)], [f_lo, fliplr(f_hi)], ...
        [0.55 0.75 0.95], ...
        'FaceAlpha', 0.45, ...
        'EdgeColor', 'none');

    % center line
    plot(xgrid, f_center, ...
        'Color', [0.00 0.45 0.74], ...
        'LineWidth', 1.8);

    % true line
    plot(xgrid, true_pdf, ...
        'Color', [0.85 0.33 0.10], ...
        'LineWidth', 2.0);

    xlim([-10 10]);
    xticks([-10 -5 0 5 10]);
    grid on;
    box on;

    xlabel('\theta');
    ylabel('pdf');
    title(figTitle, 'FontWeight', 'bold');

    legend({sprintf('Estimated %g–%g%% band', bandPct(1), bandPct(2)), ...
            centerName, ...
            'True marginal density'}, ...
            'Location', 'northwest');
end