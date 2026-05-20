%% Clear workspace
clear all;
clc;

%% === 1. Constant settings: Nature/Science-style publication aesthetics ===

% --- Core parameter: time split point ---
T_SPLIT = 60; % Time point separating fitting and prediction regions
T_MAX   = 120; % Total time

% --- Fonts and sizes ---
FONT_NAME = 'Arial'; 
FONT_SIZE_AX = 32;
FONT_SIZE_LABEL = 38;
FONT_SIZE_LEGEND = 35;

% --- Core palette: optimized Morandi broad-spectrum colors ---
COL_IFDM  = [0.88, 0.40, 0.50] ; % Dark coral pink (IFDM) - slightly deepened for contrast
COL_VMLM  = [0.45, 0.55, 0.80] ; % Hazy gray-blue (VMLM) - slightly deepened for contrast

% --- Experimental data point settings (larger markers, finer edges) ---
MARKER_SIZE = 520;     % Significantly enlarge experimental data points
MARKER_ALPHA = 0.55;   % Adjust transparency to reveal the underlying grid or background
COL_EXP_FACE = [0.60 0.60 0.60]; % Textured light gray
COL_EXP_EDGE = [0.35 0.35 0.35]; % Dark gray edge for clearer outlines

% --- Line width settings ---
LINE_WIDTH_MAIN = 10.0; % Main line width (kept visually strong for this figure size)

%% -------------------- Import data --------------------
r = xlsread("Fig.4c-d_data.csv", "A3:A80");  % Time (t), column vector
ss = xlsread("Fig.4c-d_data.csv", "B3:B80"); % Shear Stress (tau), column vector

% Keep the full data range from 0 to T_MAX.
idx_total = r <= T_MAX;
r_total = r(idx_total);
ss_total = ss(idx_total);

% Split the data into fitting and prediction regions.
idx_fit  = r_total <= T_SPLIT;
idx_pred = r_total >  T_SPLIT;

r_fit = r_total(idx_fit);      ss_fit = ss_total(idx_fit);
r_pred = r_total(idx_pred);    ss_pred = ss_total(idx_pred);

% Generate smooth plotting time vectors starting from 0.1 to avoid singularities at 0.
t_fit  = linspace(0.1, T_SPLIT, 500)';      
t_pred = linspace(T_SPLIT, T_MAX, 500)';  

%% -------------------- Model fitting and prediction --------------------
options = optimoptions('lsqcurvefit','Display','off', 'MaxFunctionEvaluations', 2000);

% === 1. VMLM (nonlinear variable order) ===
a0_vmlm = [0.1, 2000, 3, 1, 1, 0.1];
p_vmlm = lsqcurvefit(@aefun, a0_vmlm, r_fit, ss_fit, [0 0 0 0 0 0], [Inf Inf Inf Inf Inf Inf], options);

z_fit_vmlm  = aefun(p_vmlm, t_fit);  
z_pred_vmlm = aefun(p_vmlm, t_pred); 

y_true_pred = ss_pred(:);
y_pred_vmlm = aefun(p_vmlm, r_pred); y_pred_vmlm = y_pred_vmlm(:);

MAPE_vmlm_pred = mean(abs((y_true_pred - y_pred_vmlm) ./ y_true_pred)) * 100;

% === 2. IFDM (nonlinear variable order) ===
c0_ifdm = [1, 2000, 1, 0.2, 1, 0.1, 1];
p_ifdm = lsqcurvefit(@cfun, c0_ifdm, r_fit, ss_fit, [0 0 0 0 0 0 0], [Inf Inf Inf Inf 2 Inf Inf], options);

z_fit_ifdm  = cfun(p_ifdm, t_fit);   
z_pred_ifdm = cfun(p_ifdm, t_pred);  

y_pred_ifdm = cfun(p_ifdm, r_pred); y_pred_ifdm = y_pred_ifdm(:);

MAPE_ifdm_pred = mean(abs((y_true_pred - y_pred_ifdm) ./ y_true_pred)) * 100;

%% -------------------- Main figure plotting and styling --------------------
figure('Color','w','Position',[100, 100, 1150, 800]); % Slightly wider figure to leave more whitespace on the right.
hold on;

% 1. Draw refined background shading.
ylim_range = [600, 2600]; 
% Fitting-region background: very light warm gray-white to stabilize the data points visually.
patch([0 T_SPLIT T_SPLIT 0], [ylim_range(1) ylim_range(1) ylim_range(2) ylim_range(2)], ...
      [0.97 0.97 0.96], 'EdgeColor', 'none', 'HandleVisibility', 'off'); 
% Prediction-region background: very soft Morandi ice blue to suggest extrapolation/uncertainty.
patch([T_SPLIT T_MAX T_MAX T_SPLIT], [ylim_range(1) ylim_range(1) ylim_range(2) ylim_range(2)], ...
      [0.92 0.95 0.98], 'EdgeColor', 'none', 'HandleVisibility', 'off'); 

% Vertical divider: dark gray instead of pure black for a softer appearance.
xline(T_SPLIT, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2.5, 'Alpha', 0.8, 'HandleVisibility', 'off');

% Top labels: muted color so they do not dominate the figure.
text(0.25, 0.93, 'Fitting Region', 'Units', 'normalized', 'FontSize', 35, ...
    'FontWeight', 'bold', 'FontName', FONT_NAME, 'Color', 'k', 'HorizontalAlignment', 'center');
text(0.75, 0.93, 'Prediction Region', 'Units', 'normalized', 'FontSize', 35, ...
    'FontWeight', 'bold', 'FontName', FONT_NAME, 'Color', 'k', 'HorizontalAlignment', 'center');

% 2. Plot experimental data points as large circles with thin edges.
h0 = scatter(r_total, ss_total, MARKER_SIZE, ...
    'MarkerEdgeColor', COL_EXP_EDGE, 'MarkerFaceColor', COL_EXP_FACE, ...
    'LineWidth', 1.5, 'MarkerFaceAlpha', MARKER_ALPHA);

% 3. Plot VMLM with solid fitting and dashed prediction curves.
h1_fit  = plot(t_fit, z_fit_vmlm, '-', 'Color', COL_VMLM, 'LineWidth', LINE_WIDTH_MAIN);
h1_pred = plot(t_pred, z_pred_vmlm, '--', 'Color', COL_VMLM, 'LineWidth', LINE_WIDTH_MAIN, 'HandleVisibility', 'off');

% 4. Plot IFDM with solid fitting and dashed prediction curves.
h2_fit  = plot(t_fit, z_fit_ifdm, '-', 'Color', COL_IFDM, 'LineWidth', LINE_WIDTH_MAIN);
h2_pred = plot(t_pred, z_pred_ifdm, '--', 'Color', COL_IFDM, 'LineWidth', LINE_WIDTH_MAIN, 'HandleVisibility', 'off');

% 5. Accuracy text box: borderless white background with subtle transparency.
acc_str = sprintf('Prediction Accuracy:\nIFDM: MAPE = %.1f%%\nVMLM: MAPE = %.1f%%', MAPE_ifdm_pred, MAPE_vmlm_pred);
text(0.55, 0.15, acc_str, 'Units', 'normalized', 'FontSize', 24, ...
    'FontName', FONT_NAME, 'BackgroundColor', 'w', 'EdgeColor', 'none', ...
    'Margin', 12, 'Color', [0.2 0.2 0.2]);

% 6. Axis settings and styling.
ax = gca;
set(ax, 'FontName', FONT_NAME, 'FontSize', FONT_SIZE_AX, 'LineWidth', 2.5); 
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.015, 0.015], 'XMinorTick', 'on', 'YMinorTick', 'on'); 

% Axis labels (LaTeX rendering)
xlabel('$t\ (\mathrm{s})$', 'Interpreter', 'latex', 'FontSize', FONT_SIZE_LABEL);
ylabel('$\tau\ (\mathrm{Pa})$', 'Interpreter', 'latex', 'FontSize', FONT_SIZE_LABEL);
xlim([0 T_MAX]); 
ylim([600 2600]); 

% 7. Legend: remove border and place it in a suitable position.
lgd = legend([ h1_fit, h2_fit], { 'HMLM', 'IFDM'}, ...
    'Location', 'northwest', 'FontSize', FONT_SIZE_LEGEND, 'Box', 'off');
lgd.ItemTokenSize = [40, 18]; % Make legend line samples longer for readability.

%% -------------------- Model function definitions --------------------
function ae_out = aefun(a, t)
    ae_out = a(1) + (a(6) + (a(2)-a(6)) .* ((a(4).*t.^a(3)).^(a(5)./t))) .* 0.1;
end

function c_out = cfun(c, t)
    a = c(4) + (c(5)-c(4))./(1 + exp(-c(6).*(t - c(7))));
    c_out = c(1) + (c(2).*c(3).^a .* 0.1 .* t.^(1-a) ./ gamma(2-a));
end
