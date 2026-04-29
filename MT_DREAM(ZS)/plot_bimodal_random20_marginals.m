function fig = plot_bimodal_random20_marginals(samples_rep, dims20, xgrid, true_pdf, varargin)
% Plot 20 randomly selected marginal posteriors against the true 1D marginal.
%
% INPUT
%   samples_rep : Ns x d posterior sample matrix
%   dims20      : 1 x 20 selected dimensions
%   xgrid       : 1 x Nx grid for true curve
%   true_pdf    : 1 x Nx true 1D marginal density
%
% OPTIONAL
%   'NBins'       : histogram bins (default 35)
%   'FigureTitle' : figure title

    p = inputParser;
    addParameter(p, 'NBins', 35);
    addParameter(p, 'FigureTitle', 'Random 20 marginals');
    parse(p, varargin{:});

    NBins = p.Results.NBins;
    figTitle = p.Results.FigureTitle;

    if numel(dims20) ~= 20
        error('dims20 must contain exactly 20 dimensions.');
    end

    fig = figure('Color','w', 'Position',[60 60 1500 760]);
    tl = tiledlayout(4,5,'TileSpacing','compact','Padding','compact');
    title(tl, figTitle, 'FontWeight','bold');

    for k = 1:20
        ax = nexttile;
        xd = samples_rep(:, dims20(k));

        histogram(ax, xd, NBins, ...
            'Normalization', 'pdf', ...
            'FaceColor', [0.35 0.60 0.85], ...
            'EdgeColor', [0.20 0.20 0.20], ...
            'LineWidth', 0.3);
        hold(ax, 'on');

        plot(ax, xgrid, true_pdf, ...
            'Color', [0.85 0.33 0.10], ...
            'LineWidth', 1.6);

        xlim(ax, [-10 10]);
        xticks(ax, [-10 -5 0 5 10]);
        grid(ax, 'on');
        box(ax, 'on');

        xlabel(ax, sprintf('\\theta_{%d}', dims20(k)), 'Interpreter', 'tex');
        ylabel(ax, 'pdf');

        set(ax, 'FontSize', 9);
    end
end