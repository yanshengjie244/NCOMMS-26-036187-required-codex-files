%% Clear workspace
clear all;
clc;

%% === 1. Constant settings: Vibrant Soft Palette ===

% --- 1. Font settings (large Times New Roman) ---
FONT_NAME= 'Times New Roman';
FONT_SIZE= 24;
FONT_SIZE_INSET = 18;

% --- 2. Data point settings (soft marker style) ---
MARKER_SIZE= 120;
MARKER_ALPHA = 0.55;
COL_EXP_FACE= [0.60 0.60 0.60]; % Fill: neutral gray
COL_EXP_EDGE = [0.25 0.25 0.25]; % Edge: dark charcoal gray
LINE_WIDTH_EXP = 1.2;

% --- 3. Line colors (vivid but soft candy palette) ---
COL_IFDM= [0.95 0.35 0.40];% Soft Strawberry (IFDM) - solid line
COL_VFDM= [0.20 0.75 0.70];% Bright Teal (VFDM) - dashed line
COL_FDM = [0.55 0.75 0.25];% Fresh Lime (FDM alpha <= 1) - dash-dot line
COL_FDM2= [0.85 0.60 0.15];% Warm Gold (FDM alpha > 1) - dotted line
COL_SKM = [1.00 0.55 0.20];% Golden Ochre (SKM) - dotted line

% --- 4. Line widths (thicker lines for the soft colors) ---
LINE_WIDTH_MAIN = 3.5;
LINE_WIDTH_SUB = 2.5;

%% Import data
try
    % Try to read data (make sure the Excel/CSV file exists).
    r_raw = xlsread('20RS100.csv','A2:A104');
    ss_raw = xlsread('20RS100.csv','B2:B104');
    
    % Ensure that r and ss are numeric column vectors.
    r = double(r_raw(:)); % Force conversion to numeric arrays for better stability.
    ss = double(ss_raw(:));
    
    % Check for NaN values.
    if any(isnan(r)) || any(isnan(ss))
        warning('The data contain NaN values. Please check the Excel/CSV file for text or blank cells.');
        r = r(~isnan(r));
        ss = ss(~isnan(ss));
    end

catch
    disp('!! WARNING: Unable to read the data file. Random simulated data will be used instead. !!');
    % Use simulated data so the code can still run.
    r = linspace(0.1, 100, 102)';
    ss = 1000 + 500 * exp(-0.05 * r) + 500 * (r.^0.5) + 100 * randn(size(r));
end

%% Build time vector
t_data_min = min(r(r>0)); 
t = linspace(0, 100, 1000);

%% Figure window settings
figure('Color','w','Position',[100, 100, 900, 750]);
hold on;
box on;

% ---------- 1. Plot large textured experimental data points ----------
h0 = scatter(r, ss, MARKER_SIZE, ...
'MarkerEdgeColor', COL_EXP_EDGE, ...
'MarkerFaceColor', COL_EXP_FACE, ...
'LineWidth', 1.5, 'MarkerFaceAlpha', MARKER_ALPHA);

h_main = gca; % Get the main axes handle

%% -------------------- Model fitting and plotting --------------------

% --- 1. IFDM (nonlinear variable-order fractional derivative model) ---
c0 = [0.1, 1000, 1, 0.2, 0.1, 0.1, 0.1];
lb_ifdm = [0 0 0 0 0 0 0];
ub_ifdm = [Inf Inf Inf 1 2 Inf Inf];
% Ensure that input data r and ss are double during fitting.
c = lsqcurvefit(@cfun, c0, double(r), double(ss), lb_ifdm, ub_ifdm); 
alpha_nonlin_t = c(4) + (c(5)-c(4))./(1 + exp(-c(6).*(t - c(7))));
z_nonlin= c(1) + (c(2).*c(3).^alpha_nonlin_t .* 0.2 .* t.^(1-alpha_nonlin_t) ./ gamma(2-alpha_nonlin_t));
alpha_nonlin_r = c(4) + (c(5) - c(4)) ./ (1 + exp(-c(6) .* (r - c(7))));
z1=c(1)+(c(2).*c(3).^(alpha_nonlin_r).*0.2.*r.^(1-alpha_nonlin_r)./(gamma(2-alpha_nonlin_r)));
X2=z1;
Y2=ss;
A1=(sum((X2-mean(Y2)).^2))./(sum((Y2-mean(Y2)).^2));
h1 = plot(t, z_nonlin, '-', 'Color', COL_IFDM, 'LineWidth', LINE_WIDTH_MAIN);

% --- 2. VFDM (linear variable-order fractional derivative model) ---
a0 = [0.1, 5000, 0.01, 0.2, 0.1];
lb_vfdm = [0 0 0 0 0];
ub_vfdm = [Inf Inf Inf 0.01 1.5];
a = lsqcurvefit(@afun, a0, double(r), double(ss), lb_vfdm, ub_vfdm);
alpha_linear_t = a(4).*t + a(5);
z_linear= a(1) + (a(2).*a(3).^alpha_linear_t .* 0.2 .* t.^(1-alpha_linear_t) ./ gamma(2-alpha_linear_t));
h2 = plot(t, z_linear, '--', 'Color', COL_VFDM, 'LineWidth', LINE_WIDTH_MAIN);

% --- 3. FDM (fixed order, order constraint $\alpha \in [0, 1]$) ---
b0 = [0.1,1000,0.001, 0.1];
lb_fdm_0_1 = [0 0 0 0];
ub_fdm_0_1 = [Inf Inf Inf 2];
b = lsqcurvefit(@bfun, b0, double(r), double(ss), lb_fdm_0_1, ub_fdm_0_1);
alpha_fixed_0_1 = b(4);
z_fixed_0_1 = b(1) + (b(2).*b(3).^alpha_fixed_0_1 .* 0.2 .* t.^(1-alpha_fixed_0_1) ./ gamma(2-alpha_fixed_0_1));
h3 = plot(t, z_fixed_0_1, '-.', 'Color', COL_FDM, 'LineWidth', LINE_WIDTH_MAIN);

% --- 4. FDM (fixed order, order constraint $\alpha \in [1, 2]$) ---

% --- 5. SKM (structural kinetic model) ---
d0 = [0.1, 0.01,6000, 0.2, 1];
d = lsqcurvefit(@dfun, d0, double(r), double(ss), [], [Inf Inf Inf Inf Inf]);
z_skm = ((d(2)./(d(1).*0.2) + d(2)) + (1 - (d(2)./(d(1).*0.2) + d(2))).*exp(-(d(1).*0.2 + d(2)).*t)) .* (d(3) + d(4).*0.2.^d(5));
h5 = plot(t, z_skm, ':', 'Color', COL_SKM, 'LineWidth', LINE_WIDTH_MAIN);


%% -------------------- Main figure styling --------------------
axes(h_main);

% Axis labels (Times New Roman + LaTeX)
xlabel('$$t \ \rm{(s)}$$', 'Interpreter','latex', 'FontSize', FONT_SIZE+4, 'FontName', FONT_NAME);
ylabel('$$\tau \ \rm{(Pa)}$$', 'Interpreter','latex', 'FontSize', FONT_SIZE+4, 'FontName', FONT_NAME);
xlim([0 100]);
ylim([0 5200]);
set(h_main, ...
'FontName', FONT_NAME, ...
'FontSize', FONT_SIZE, ...
'LineWidth', 2.0, ...
'TickDir','in', ...
'TickLength',[0.015 0.015], ...
'XMinorTick','on', ...
'YMinorTick','on', ...
'Box','on');

% Legend
h_all = [h0 h1 h2 h3 h5];
legend_text = {'Test data', 'IFDM', 'VFDM', 'FDM', 'SKM'}; 
[lgd, icons] = legend(h_all, legend_text, ...
'Location','northwest', ...
'Interpreter','latex', ... 
'FontSize', FONT_SIZE-2, ...
'FontName', FONT_NAME, ...
'Box','off');

% Adjust data point size in the legend to fix undersized scatter markers.
lgd.ItemTokenSize = [30, 18];
h_icon_patch = findobj(icons, 'Type', 'patch');
set(h_icon_patch, 'MarkerSize', 14);
set(h_icon_patch, 'LineWidth', 1.5);

% Panel label (a)
text(0.02, 0.95, '\bf{a}', 'Units','normalized', ...
'FontName', FONT_NAME, 'FontSize', FONT_SIZE+6);

%% -------------------- Inset: order trend alpha(t) --------------------
ax_inset = axes('Position',[0.52 0.45 0.40 0.40]);
box on; hold on;

% Plot inset curves
p1 = plot(t, alpha_nonlin_t, '-', 'Color', COL_IFDM, 'LineWidth', LINE_WIDTH_SUB);
p2 = plot(t, alpha_linear_t, '--', 'Color', COL_VFDM, 'LineWidth', LINE_WIDTH_SUB);
p3 = plot(t, alpha_fixed_0_1 * ones(size(t)), '-.', 'Color', COL_FDM, 'LineWidth', LINE_WIDTH_SUB);


% --- Build alpha_all using a stable step-by-step assignment method ---
N = length(t);
alpha_all_size = 4 * N; 
alpha_all = zeros(alpha_all_size, 1); 

alpha_all(1:N) = alpha_nonlin_t(:);
alpha_all((N+1):(2*N)) = alpha_linear_t(:);
alpha_all((2*N+1):(3*N)) = alpha_fixed_0_1 * ones(N, 1);
% --- End alpha_all construction ---

% Precompute the Y-axis range and margin for text positioning.
y_margin = 0.1 * (max(alpha_all) - min(alpha_all) + eps);

% Automatically adjust the Y-axis range and fix the X-axis so guide lines reach the edges correctly.
ylim([min(alpha_all) - y_margin, max(alpha_all) + y_margin]);
xlim([0 100]);

% Get inset axes edge values.
x_min_inset = min(xlim(ax_inset));
y_min_inset = min(ylim(ax_inset));


%% === Annotate key times of the IFDM model ===

% --- 1. Inflection/center time t_c (50% transition) ---
t_c = c(7);
alpha_at_tc = c(4) + (c(5)-c(4))./(1 + exp(-c(6).*(t_c - c(7))));

% 1.1. Plot the inflection marker (pentagram, slightly larger MarkerSize).
plot(t_c, alpha_at_tc, 'p', 'MarkerSize', 14, 'MarkerFaceColor', COL_IFDM, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% 1.2. Plot the vertical guide line to the X-axis.
plot([t_c t_c], [y_min_inset alpha_at_tc], '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);

% 1.3. Plot the horizontal guide line to the Y-axis.
plot([x_min_inset t_c], [alpha_at_tc alpha_at_tc], '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);

% 1.4. Annotate center time $t_c$ with raised text to avoid overlap with t_stab.
text(t_c, y_min_inset + y_margin/1.2, ['$$t_c = ' num2str(t_c, '%.1f') '$$ s'], ...
    'Interpreter', 'latex', 'FontSize', FONT_SIZE_INSET-4, 'FontName', FONT_NAME, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% 1.5. Annotate center order $\alpha(t_c)$ near the Y-axis.
text(x_min_inset + (max(xlim(ax_inset))-x_min_inset)*0.05, alpha_at_tc, ['$$\alpha(t_c) = ' num2str(alpha_at_tc, '%.3f') '$$'], ...
    'Interpreter', 'latex', 'FontSize', FONT_SIZE_INSET-6, 'FontName', FONT_NAME, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');


% --- 2. Stabilization time t_stab (99% convergence) ---
P_stab = 0.99; % 99% convergence toward stability
if c(6) ~= 0
    R_frac = P_stab; 
    % Calculate stabilization time t_stab at 99% of the range between c4 and c5.
    t_stab = c(7) - (1/c(6)) * log((1 - R_frac) / R_frac);
    
    % Alpha value at the stabilization time for plotting.
    alpha_stab = c(4) + (c(5)-c(4))./(1 + exp(-c(6).*(t_stab - c(7))));
else
    t_stab = t(end);
    alpha_stab = c(5);
end

% If t_stab exceeds the plotting range, draw only up to t = 105 s.
if t_stab > 105
    t_stab_plot = 105; % Draw to the boundary
    alpha_stab_plot = c(4) + (c(5)-c(4))./(1 + exp(-c(6).*(t_stab_plot - c(7))));
    t_stab_text = ['$$t_{stab} > ' num2str(105, '%.1f') '$$ s']; % Text shows t > 105
else
    t_stab_plot = t_stab;
    alpha_stab_plot = alpha_stab;
    t_stab_text = ['$$t_{stab} = ' num2str(t_stab, '%.1f') '$$ s'];
end

% 2.1. Plot the stabilization marker as a filled circle.
plot(t_stab_plot, alpha_stab_plot, 'o', 'MarkerSize', 10, 'MarkerFaceColor', COL_IFDM, 'MarkerEdgeColor', 'k');

% 2.2. Plot the vertical guide line to the X-axis.
plot([t_stab_plot t_stab_plot], [y_min_inset alpha_stab_plot], ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0);

% 2.3. Plot the horizontal guide line to the Y-axis.
plot([x_min_inset t_stab_plot], [alpha_stab_plot alpha_stab_plot], ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0);

% 2.4. Annotate stabilization time $t_{stab}$ with text near the X-axis.
text(t_stab_plot, y_min_inset + y_margin/8, t_stab_text, ...
    'Interpreter', 'latex', 'FontSize', FONT_SIZE_INSET-6, 'FontName', FONT_NAME, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% 2.5. Annotate stabilization order $\alpha_{stab}$ near the Y-axis.
text(x_min_inset + (max(xlim(ax_inset))-x_min_inset)*0.05, alpha_stab_plot, ['$$\alpha_{stab} \approx ' num2str(alpha_stab, '%.3f') '$$'], ...
    'Interpreter', 'latex', 'FontSize', FONT_SIZE_INSET-8, 'FontName', FONT_NAME, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');


%% --- Inset styling ---
set(ax_inset, ...
'FontName', FONT_NAME, ...
'FontSize', FONT_SIZE_INSET, ...
'LineWidth', 1.5, ...
'TickDir','in', ...
'XMinorTick','on', ...
'YMinorTick','on', ...
'Color', 'none');

xlabel('$$t \ \rm{(s)}$$', 'Interpreter','latex','FontName', FONT_NAME, 'FontSize', FONT_SIZE_INSET);
ylabel('$$\alpha(t)$$', 'Interpreter','latex','FontName', FONT_NAME, 'FontSize', FONT_SIZE_INSET);

% Inset legend
legend([p1 p2 p3], {'IFDM','VFDM','FDM'}, ...
'Location','southwest', ...
'Interpreter','latex', ...
'FontSize', FONT_SIZE_INSET-4, ...
'FontName', FONT_NAME, ...
'Box','off');

%% -------------------- Fitting function definitions (unchanged) --------------------
function c_out = cfun(c, t)
alpha = c(4) + (c(5)-c(4))./(1 + exp(-c(6).*(t - c(7))));
c_out = c(1) + (c(2).*c(3).^alpha .* 0.2 .* t.^(1-alpha) ./ gamma(2-alpha));
end

function b_out = bfun(b, t)
alpha_1 = b(4);
b_out= b(1) + (b(2).*b(3).^alpha_1 .* 0.2 .* t.^(1-alpha_1) ./ gamma(2-alpha_1));
end
function a_out = afun(a, t)
alpha_2 = a(4).*t + a(5);
a_out= a(1) + (a(2).*a(3).^alpha_2 .* 0.2 .* t.^(1-alpha_2) ./ gamma(2-alpha_2));
end

function d_out = dfun(d, t)
d_out = ((d(2)./(d(1).*0.2) + d(2)) + ...
(1 - (d(2)./(d(1).*0.2) + d(2))).*exp(-(d(1).*0.2 + d(2)).*t)) .* ...
(d(3) + d(4).*0.2.^d(5));
end
