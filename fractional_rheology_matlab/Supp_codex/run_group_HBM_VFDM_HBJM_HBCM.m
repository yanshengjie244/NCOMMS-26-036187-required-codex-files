function run_group_HBM_VFDM_HBJM_HBCM(groupId, dataFile)
% Fit one two-column rheology group with HBM, VFDM, HBJM, and HBCM.
% Figures are displayed only; no files are saved.

clearvars -except groupId dataFile; clc; close all;

if nargin < 2 || isempty(dataFile)
    DATA_FILE = char([27969 21464 25968 25454 46 120 108 115 120]); % Encoded legacy rheology data filename
else
    DATA_FILE = dataFile;
end
NUM_STARTS = 8;
RANDOM_SEED = 20260513;
USE_LOG_X = false;
FONT_SIZE = 35;

rng(RANDOM_SEED);

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

if exist(DATA_FILE, 'file') == 2
    dataPath = DATA_FILE;
elseif exist(fullfile(scriptDir, DATA_FILE), 'file') == 2
    dataPath = fullfile(scriptDir, DATA_FILE);
else
    error('Cannot find data file: %s', DATA_FILE);
end

numData = load_numeric_matrix(dataPath);
if size(numData, 2) < 2 * groupId
    error('Group %d requires columns %d-%d, but the data file has only %d columns.', ...
        groupId, 2 * groupId - 1, 2 * groupId, size(numData, 2));
end

raw = numData(:, (2 * groupId - 1):(2 * groupId));
raw = raw(all(isfinite(raw), 2), :);
raw = sortrows(raw, 1);

gd = raw(:, 1);
tau = raw(:, 2);
valid = gd > 0 & tau > 0;
gd = gd(valid);
tau = tau(valid);

if numel(gd) < 6
    error('Group %d has fewer than six valid points.', groupId);
end

if USE_LOG_X
    gdFit = logspace(log10(min(gd)), log10(max(gd)), 700).';
else
    gdFit = linspace(0, max(gd), 700).';
end

fitVFDM = fit_model('VFDM', @(p, x) model_vfdm(p, x), ...
    initial_vfdm(gd, tau), lower_vfdm(tau, gd), upper_vfdm(gd, tau), gd, tau, NUM_STARTS);

fitHBJM = fit_model('HBJM', @(p, x) model_hbjm(p, x), ...
    initial_hbjm(gd, tau), lower_hbjm(tau), upper_hbjm(gd, tau), gd, tau, NUM_STARTS);

fitHBCM = fit_model('HBCM', @(p, x) model_hbcm(p, x), ...
    initial_hbcm(gd, tau), lower_hbcm(tau), upper_hbcm(gd, tau), gd, tau, NUM_STARTS);

fits = [fitVFDM, fitHBJM, fitHBCM];
statsTable = build_stats_table(fits, tau);
paramTable = build_param_table(fits, groupId);

fprintf('\n%s, Group %d completed.\n', DATA_FILE, groupId);
disp(statsTable);
disp(paramTable);
diagnose_parameters(fits, gd);

plot_fit_figure(groupId, gd, tau, gdFit, fits, USE_LOG_X, FONT_SIZE);
plot_residual_figure(groupId, gd, fits, USE_LOG_X, FONT_SIZE);
plot_order_figure(groupId, gd, gdFit, fits, USE_LOG_X, FONT_SIZE);

fprintf('\nFigures are displayed only. No image or figure files were saved.\n');
end

%% ===================== Data import =====================
function numData = load_numeric_matrix(fileName)
    [~, ~, ext] = fileparts(fileName);
    if strcmpi(ext, '.csv')
        imported = importdata(fileName);
        if isstruct(imported)
            numData = imported.data;
        else
            numData = imported;
        end
    else
        [numData, ~, ~] = xlsread(fileName);
    end
    if isempty(numData)
        error('No numeric data were found in %s.', fileName);
    end
end

%% ===================== Models =====================
function y = model_vfdm(p, gd)
    % p = [tau_y0, K_alpha, alpha_0, alpha_inf, k, gamma_dot_c]
    % gamma_dot is used directly; no dimensionless shear-rate normalization.
    alpha = alpha_sigmoid(p(3:6), gd);
    y = p(1) + p(2) .* gd .^ alpha ./ gamma(2 - alpha);
end

function y = model_hbjm(p, gd)
    % p = [tau_0, eta_inf, eta_0, a, b, m, n]
    eta = p(2) + (p(3) - p(2)) .* p(5) ./ (p(4) + p(5) .* gd .^ p(6));
    y = p(1) + eta .* gd .^ p(7);
end

function y = model_hbcm(p, gd)
    % p = [tau_0, eta_inf, delta_eta, theta, beta, delta, alpha_3]
    etaInf = p(2);
    eta0 = p(2) + p(3);
    eta = etaInf + (eta0 - etaInf) .* (1 + (p(4) .* gd) .^ p(5)) .^ ((p(6) - 1) ./ p(5));
    y = p(1) + eta .* gd .^ p(7);
end

function alpha = alpha_sigmoid(q, gd)
    z = -q(3) .* (gd - q(4));
    z = max(min(z, 60), -60);
    alpha = q(1) + (q(2) - q(1)) ./ (1 + exp(z));
    alpha = min(max(alpha, 1e-4), 0.9999);
end

%% ===================== Fitting =====================
function fit = fit_model(name, modelFun, p0, lb, ub, x, y, nStarts)
    starts = make_starts(p0, lb, ub, nStarts);
    bestSSE = inf;
    bestP = p0;
    exitFlag = NaN;

    hasLsqcurvefit = exist('lsqcurvefit', 'file') == 2;
    if hasLsqcurvefit
        opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
            'MaxFunctionEvaluations', 2e4, 'MaxIterations', 2e3, ...
            'FunctionTolerance', 1e-11, 'StepTolerance', 1e-11);
    else
        opts = optimset('Display', 'off', 'MaxFunEvals', 2e4, 'MaxIter', 2e3);
    end

    for i = 1:size(starts, 1)
        try
            if hasLsqcurvefit
                [p, ~, residual, flag] = lsqcurvefit(modelFun, starts(i, :), x, y, lb, ub, opts);
                sse = sum(residual .^ 2);
            else
                obj = @(pp) sum((modelFun(bound_params(pp, lb, ub), x) - y) .^ 2);
                [pRaw, sse, flag] = fminsearch(obj, starts(i, :), opts);
                p = bound_params(pRaw, lb, ub);
            end
            if isfinite(sse) && sse < bestSSE
                bestSSE = sse;
                bestP = p;
                exitFlag = flag;
            end
        catch
        end
    end

    yhat = modelFun(bestP, x);
    fit.name = name;
    fit.params = bestP;
    fit.yhat = yhat;
    fit.residual = y - yhat;
    fit.sse = sum((y - yhat) .^ 2);
    fit.exitFlag = exitFlag;
    fit.modelFun = modelFun;
end

function p = bound_params(p, lb, ub)
    p = min(max(p, lb), ub);
end

function starts = make_starts(p0, lb, ub, nStarts)
    n = numel(p0);
    starts = zeros(nStarts, n);
    starts(1, :) = p0;
    finiteSpan = isfinite(lb) & isfinite(ub);
    for i = 2:nStarts
        p = p0;
        for j = 1:n
            if finiteSpan(j)
                if lb(j) >= 0 && ub(j) > 0 && ub(j) / max(lb(j), eps) > 100
                    lo = log10(max(lb(j), eps));
                    hi = log10(max(ub(j), eps));
                    p(j) = 10 ^ (lo + rand * (hi - lo));
                else
                    p(j) = lb(j) + rand * (ub(j) - lb(j));
                end
            end
        end
        starts(i, :) = p;
    end
end

%% ===================== Initial values and bounds =====================
function p0 = initial_vfdm(gd, tau)
    tauY = max(0, min(tau) * 0.75);
    k = max(range(tau), max(tau) * 0.05);
    p0 = [tauY, k, 0.90, 0.20, 1 / max(range(gd), eps), median(gd)];
end

function lb = lower_vfdm(tau, gd)
    lb = [0, 0, 0.01, 0.01, 1e-5, min(gd)];
end

function ub = upper_vfdm(gd, tau)
    ub = [max(tau) * 1.2, max(tau) * 500, 0.99, 0.99, 100 / max(range(gd), eps), max(gd)];
end

function p0 = initial_hbjm(gd, tau)
    tauY = max(0, min(tau) * 0.75);
    etaInf = max(range(tau) / max(max(gd), eps), eps);
    eta0 = max(max(tau) / max(min(gd), eps), etaInf * 2);
    p0 = [tauY, etaInf, eta0, 1, 1, 1, 0.5];
end

function lb = lower_hbjm(tau)
    lb = [0, 0, 0, 1e-6, 1e-6, 0.05, 0.01];
end

function ub = upper_hbjm(gd, tau)
    ub = [max(tau) * 1.2, max(tau) * 200, max(tau) * 500, 1e6, 1e6, 3.0, 2.0];
end

function p0 = initial_hbcm(gd, tau)
    tauY = max(0, min(tau) * 0.75);
    etaInf = max(range(tau) / max(max(gd), eps), eps);
    deltaEta = max(max(tau) / max(min(gd), eps) - etaInf, etaInf);
    p0 = [tauY, etaInf, deltaEta, 1 / max(median(gd), eps), 1, 0.3, 0.5];
end

function lb = lower_hbcm(tau)
    lb = [0, 0, 0, 1e-4, 0.05, 0.01, 0.01];
end

function ub = upper_hbcm(gd, tau)
    ub = [max(tau) * 1.2, max(tau) * 200, max(tau) * 500, 1e4, 5.0, 2.0, 2.0];
end

%% ===================== Tables and diagnostics =====================
function stats = build_stats_table(fits, y)
    names = cell(numel(fits), 1);
    sse = zeros(numel(fits), 1);
    mse = zeros(numel(fits), 1);
    mae = zeros(numel(fits), 1);
    r2 = zeros(numel(fits), 1);
    n = numel(y);
    sst = sum((y - mean(y)) .^ 2);
    for i = 1:numel(fits)
        names{i} = fits(i).name;
        sse(i) = fits(i).sse;
        mse(i) = sse(i) / n;
        mae(i) = mean(abs(fits(i).residual));
        r2(i) = 1 - sse(i) / max(sst, eps);
    end
    stats = table(names, sse, mse, mae, r2, ...
        'VariableNames', {'Model', 'SSE', 'MSE', 'MAE', 'R2'});
end

function params = build_param_table(fits, groupId)
    rows = {};
    for i = 1:numel(fits)
        p = fits(i).params;
        switch fits(i).name
            case 'VFDM'
                labels = {'tau_0', 'K_alpha', 'alpha_0', 'alpha_inf', 'k', 'gamma_dot_c'};
            case 'HBJM'
                labels = {'tau_0', 'eta_inf', 'eta_0', 'a', 'b', 'm', 'n'};
            case 'HBCM'
                labels = {'tau_0', 'eta_inf', 'delta_eta', 'theta', 'beta', 'delta', 'alpha_3'};
        end
        for j = 1:numel(p)
            rows(end + 1, :) = {groupId, fits(i).name, labels{j}, p(j)}; %#ok<AGROW>
        end
    end
    params = cell2table(rows, 'VariableNames', {'Group', 'Model', 'Parameter', 'Value'});
end

function diagnose_parameters(fits, gd)
    fprintf('Parameter check:\n');
    for i = 1:numel(fits)
        p = fits(i).params;
        switch fits(i).name
            case 'VFDM'
                alpha = alpha_sigmoid(p(3:6), gd);
                if abs(p(3) - 0.99) < 1e-3 || abs(p(4) - 0.01) < 1e-3
                    fprintf('  - VFDM alpha limit is near a bound: alpha_0 = %.4g, alpha_inf = %.4g.\n', p(3), p(4));
                end
                if range(alpha) < 0.03
                    fprintf('  - VFDM alpha changes weakly over the tested shear-rate range.\n');
                end
            case 'HBJM'
                if p(4) > 1e5 || p(5) > 1e5
                    fprintf('  - HBJM a or b is very large; the transition term may be weakly identifiable.\n');
                end
                if p(7) > 1.95 || p(7) < 0.03
                    fprintf('  - HBJM n = %.4g is close to a fitting bound.\n', p(7));
                end
            case 'HBCM'
                if p(4) > 1e3 || p(5) > 4.9
                    fprintf('  - HBCM theta or beta is near a high bound; check parameter identifiability.\n');
                end
                if p(6) < 0.02 || p(6) > 1.95
                    fprintf('  - HBCM delta = %.4g is close to a fitting bound.\n', p(6));
                end
        end
    end
end

%% ===================== Plotting =====================
function plot_fit_figure(groupId, gd, tau, gdFit, fits, useLogX, fontSize)
    colors = model_colors();
    [lineStyles, lineWidth] = model_line_styles();
    fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 28, 21], 'Renderer', 'painters');
    ax = axes('Position', [0.16, 0.16, 0.78, 0.74]);
    hold(ax, 'on');
    scatter(ax, gd, tau, 160, 'o', 'MarkerFaceColor', colors.exp, ...
        'MarkerEdgeColor', 'w', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    plot(ax, NaN, NaN, 'o', ...
        'MarkerSize', 18, ...
        'MarkerFaceColor', colors.exp, ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 1.2, ...
        'DisplayName', 'Experimental data');
    for i = 1:numel(fits)
        style = lineStyles.(lower(fits(i).name));
        plot(ax, gdFit, fits(i).modelFun(fits(i).params, gdFit), ...
            'Color', colors.(lower(fits(i).name)), 'LineStyle', style, ...
            'LineWidth', lineWidth, 'DisplayName', fits(i).name);
    end
    format_axes(ax, useLogX, fontSize);
    xlabel(ax, 'Shear rate (s^{-1})', 'FontSize', fontSize);
    ylabel(ax, '\tau (Pa)', 'FontSize', fontSize);
    xlim(ax, get_x_limits(gd, useLogX));
    ylim_with_margin(ax, [tau; all_predictions(fits, gdFit)]);
    lgd = legend(ax, 'Location', 'best');
    set(lgd, 'Box', 'off', 'FontSize', fontSize);
    enlarge_legend_markers(lgd, 18);
    set(fig, 'Name', sprintf('Group %d fit', groupId));
    drawnow;
end

function plot_residual_figure(groupId, gd, fits, useLogX, fontSize)
    colors = model_colors();
    [lineStyles, lineWidth] = model_line_styles();
    fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [3, 3, 28, 21], 'Renderer', 'painters');
    ax = axes('Position', [0.16, 0.16, 0.78, 0.74]);
    hold(ax, 'on');
    for i = 1:numel(fits)
        plot(ax, gd, fits(i).residual, 'o', ...
            'Color', colors.(lower(fits(i).name)), ...
            'MarkerFaceColor', colors.(lower(fits(i).name)), ...
            'MarkerEdgeColor', 'w', 'MarkerSize', 9, ...
            'LineStyle', lineStyles.(lower(fits(i).name)), ...
            'LineWidth', lineWidth * 0.65, 'DisplayName', fits(i).name);
    end
    plot(ax, [min(gd), max(gd)], [0, 0], ':', 'Color', [0.35, 0.35, 0.35], ...
        'LineWidth', 2.0, 'HandleVisibility', 'off');
    format_axes(ax, useLogX, fontSize);
    xlabel(ax, 'Shear rate (s^{-1})', 'FontSize', fontSize);
    ylabel(ax, '\Delta\tau (Pa)', 'FontSize', fontSize);
    xlim(ax, get_x_limits(gd, useLogX));
    lgd = legend(ax, 'Location', 'best');
    set(lgd, 'Box', 'off', 'FontSize', fontSize);
    set(fig, 'Name', sprintf('Group %d residual', groupId));
    drawnow;
end

function plot_order_figure(groupId, gd, gdFit, fits, useLogX, fontSize)
    colors = model_colors();
    [lineStyles, lineWidth] = model_line_styles();
    fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [4, 4, 28, 21], 'Renderer', 'painters');
    ax = axes('Position', [0.16, 0.16, 0.78, 0.74]);
    hold(ax, 'on');
    orderData = [];
    vfdmParams = [];
    for i = 1:numel(fits)
        p = fits(i).params;
        switch fits(i).name
            case 'VFDM'
                order = alpha_sigmoid(p(3:6), gdFit);
                vfdmParams = p;
            otherwise
                continue;
        end
        orderData = [orderData; order(:)]; %#ok<AGROW>
        plot(ax, gdFit, order, 'Color', colors.(lower(fits(i).name)), ...
            'LineStyle', lineStyles.(lower(fits(i).name)), ...
            'LineWidth', lineWidth, 'DisplayName', fits(i).name);
    end
    format_axes(ax, useLogX, fontSize);
    xlabel(ax, 'Shear rate (s^{-1})', 'FontSize', fontSize);
    ylabel(ax, 'Order', 'FontSize', fontSize);
    xlim(ax, get_x_limits(gd, useLogX));
    ylim(ax, order_limits(orderData));
    add_gamma_c_marker(ax, vfdmParams);
    add_vfdm_parameter_text(ax, vfdmParams, fontSize);
    lgd = legend(ax, 'Location', 'best');
    set(lgd, 'Box', 'off', 'FontSize', fontSize);
    set(fig, 'Name', sprintf('Group %d order', groupId));
    drawnow;
end

function colors = model_colors()
    colors.exp = [0.08, 0.08, 0.08];
    colors.vfdm = [0.85, 0.33, 0.10];
    colors.hbjm = [0.00, 0.45, 0.74];
    colors.hbcm = [0.49, 0.18, 0.56];
end

function [lineStyles, lineWidth] = model_line_styles()
    lineStyles.vfdm = '-';
    lineStyles.hbjm = '-.';
    lineStyles.hbcm = ':';
    lineWidth = 4.0;
end

function format_axes(ax, useLogX, fontSize)
    set(ax, 'FontName', 'Arial', 'FontSize', fontSize, 'LineWidth', 2.0, ...
        'Box', 'on', 'TickDir', 'out', 'TickLength', [0.018, 0.018], ...
        'Layer', 'top');
    if useLogX
        set(ax, 'XScale', 'log');
    else
        set(ax, 'XScale', 'linear');
    end
    grid(ax, 'off');
end

function xl = get_x_limits(gd, useLogX)
    if useLogX
        xl = [min(gd), max(gd)];
    else
        xl = [0, max(gd)];
    end
end

function add_vfdm_parameter_text(ax, p, fontSize)
    if isempty(p)
        return;
    end
    txt = { ...
        sprintf('$\\alpha_0 = %.3f$', p(3)), ...
        sprintf('$\\alpha_{\\infty} = %.3f$', p(4)), ...
        sprintf('$k = %.3g\\ \\mathrm{s}$', p(5)), ...
        sprintf('$\\dot{\\gamma}_c = %.3g\\ \\mathrm{s}^{-1}$', p(6))};
    text(ax, 0.62, 0.92, txt, ...
        'Units', 'normalized', ...
        'Interpreter', 'latex', ...
        'VerticalAlignment', 'top', ...
        'FontName', 'Arial', ...
        'FontSize', max(18, round(fontSize * 0.68)), ...
        'BackgroundColor', 'w', ...
        'Margin', 6, ...
        'EdgeColor', [0.82, 0.82, 0.82]);
end

function enlarge_legend_markers(lgd, markerSize)
    icons = findobj(lgd, 'Type', 'line');
    for i = 1:numel(icons)
        marker = get(icons(i), 'Marker');
        if ~strcmp(marker, 'none')
            set(icons(i), 'MarkerSize', markerSize);
        end
    end
end

function add_gamma_c_marker(ax, p)
    if isempty(p)
        return;
    end
    gammaC = p(6);
    alphaC = alpha_sigmoid(p(3:6), gammaC);
    plot(ax, gammaC, alphaC, 'o', ...
        'MarkerSize', 12, ...
        'MarkerFaceColor', [0.85, 0.33, 0.10], ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 1.4, ...
        'HandleVisibility', 'off');
end

function y = all_predictions(fits, x)
    y = [];
    for i = 1:numel(fits)
        y = [y; fits(i).modelFun(fits(i).params, x)]; %#ok<AGROW>
    end
end

function ylim_with_margin(ax, y)
    y = y(isfinite(y));
    dy = max(y) - min(y);
    if dy <= eps
        dy = max(abs(mean(y)), 1) * 0.1;
    end
    ylim(ax, [min(y) - 0.08 * dy, max(y) + 0.12 * dy]);
end

function yl = order_limits(orderData)
    orderData = orderData(isfinite(orderData));
    if isempty(orderData)
        yl = [0, 1];
        return;
    end
    ymin = min(orderData);
    ymax = max(orderData);
    dy = max(ymax - ymin, 0.08);
    yl = [max(0, ymin - 0.15 * dy), min(2.05, ymax + 0.18 * dy)];
end
