%% Clear workspace
clear all;
clc;

%% === 1. Constant settings: matched to the second Nature/Science style section ===

% --- Fonts and sizes ---
FONT_NAME = 'Arial'; 
FONT_SIZE_AX = 40;
FONT_SIZE_LABEL = 45;
FONT_SIZE_LEGEND = 40;

% --- Core palette: Morandi broad-spectrum colors matched to the second section ---
COL_IFDM  = [1, 0.6, 0.78] ; % Light blue
COL_VFDM  = [1, 0.84, 0] ; % Light green
COL_FDM   = [0.68, 0.92, 1] ; % Light pink
COL_SKM   = [0.85, 0.7, 1] ; % Coral
COL_OTHER = [166, 166, 210] / 255; % Light purple (backup)

% --- Experimental data point settings (muted gray to emphasize curves) ---
MARKER_SIZE = 120;
MARKER_ALPHA =0.01;
COL_EXP_FACE = [0.70 0.70 0.70]; % Neutral light gray
COL_EXP_EDGE = [0.70 0.70 0.70]; % Dark gray edge

% --- Line width settings ---
LINE_WIDTH_MAIN = 10.0; % Extra-thick line matched to the second section

%% Import data
r = xlsread("csb30.csv", "A3:A80"); % Time (t)
ss = xlsread("csb30.csv", "B3:B80"); % Shear Stress (ta
 
t = linspace(1, 120, 1000);

%% -------------------- Model fitting and plotting --------------------
figure('Color','w','Position',[100, 100, 1000, 850]);
hold on;

% Plot experimental data points
h0 = scatter(r, ss, MARKER_SIZE, ...
    'MarkerEdgeColor', COL_EXP_EDGE, ...
    'MarkerFaceColor', COL_EXP_FACE, ...
    'LineWidth', 18.0, 'MarkerFaceAlpha', MARKER_ALPHA);

% --- 1. IFDM (nonlinear variable order) ---
c = lsqcurvefit(@cfun, [1,200, 1, 0.2, 1, 0.1, 5], r, ss, [0 0 0 0 0 0 0], [Inf Inf Inf Inf 2 Inf Inf]); 
alpha_nonlin_t = c(4) + (c(5)-c(4))./(1 + exp(-c(6).*(t - c(7))));
z_nonlin = c(1) + (c(2).*c(3).^alpha_nonlin_t .* 0.1 .* t.^(1-alpha_nonlin_t) ./ gamma(2-alpha_nonlin_t));
alpha_nonlin_r = c(4) + (c(5)-c(4))./(1 + exp(-c(6).*(r - c(7))));
z3=c(1) + (c(2).*c(3).^alpha_nonlin_r .* 0.1 .* r.^(1-alpha_nonlin_r) ./ gamma(2-alpha_nonlin_r));
X1=z3
Y1=ss
A0=(sum((X1-mean(Y1)).^2))./(sum((Y1-mean(Y1)).^2))
h1 = plot(t, z_nonlin, '-', 'Color', COL_IFDM, 'LineWidth', LINE_WIDTH_MAIN);

% --- 2. VMLM (nonlinear variable order) ---
a0=[0.1,1,4,1,3,0.1];
a=lsqcurvefit(@aefun,a0,r,ss,[0 0 0 0 0 0],[Inf Inf Inf Inf Inf]);
z=a(1)+(a(6)+(a(2)-a(6)).*((a(4).*t.^a(3)).^(a(5)./t))).*0.1
z1=a(1)+(a(6)+(a(2)-a(6)).*((a(4).*r.^a(3)).^(a(5)./r))).*0.1
X1=z1
Y1=ss
A1=(sum((X1-mean(Y1)).^2))./(sum((Y1-mean(Y1)).^2))
%h2 = plot(t, z, '-', 'Color', COL_SKM, 'LineWidth', LINE_WIDTH_MAIN);
% --- 2. VFDM (linear variable order) ---
a = lsqcurvefit(@afun, [0.1, 100, 0.01, 0.2, 0.1], r, ss, [0 0 0 -2 0], [Inf Inf Inf Inf Inf]);
alpha_linear_t = a(4).*t + a(5);
z_linear = a(1) + (a(2).*a(3).^alpha_linear_t .* 0.1 .* t.^(1-alpha_linear_t) ./ gamma(2-alpha_linear_t));
alpha_linear_r = a(4).*r + a(5);
z2=a(1) + (a(2).*a(3).^alpha_linear_r .* 0.1 .* r.^(1-alpha_linear_r) ./ gamma(2-alpha_linear_r));
X1=z2
Y1=ss
A2=(sum((X1-mean(Y1)).^2))./(sum((Y1-mean(Y1)).^2))
h3 = plot(t, z_linear, '-', 'Color', COL_VFDM, 'LineWidth', LINE_WIDTH_MAIN);

% --- 3. FDM (constant fractional order) ---
b = lsqcurvefit(@bfun, [0.1,5000,0.1, 0], r, ss, [0 0 0 1.1], [Inf Inf Inf Inf]);
alpha_fixed_0_1 = b(4);
z_fixed_0_1 = b(1) + (b(2).*b(3).^alpha_fixed_0_1 .* 0.1 .* t.^(1-alpha_fixed_0_1) ./ gamma(2-alpha_fixed_0_1));
z4=b(1) + (b(2).*b(3).^alpha_fixed_0_1 .* 0.1 .* r.^(1-alpha_fixed_0_1) ./ gamma(2-alpha_fixed_0_1));
X1=z4
Y1=ss
A3=(sum((X1-mean(Y1)).^2))./(sum((Y1-mean(Y1)).^2))
h4 = plot(t, z_fixed_0_1, '-', 'Color', COL_FDM, 'LineWidth', LINE_WIDTH_MAIN);

% --- 5. SKM (structural kinetic model) ---
d = lsqcurvefit(@dfun, [0.001, 0.01, 1800, 0.2, 1], r, ss, [], []);
z_skm = ((d(2)./(d(1).*100) + d(2)) + (1 - (d(2)./(d(1).*100) + d(2))).*exp(-(d(1).*100 + d(2)).*t)) .* (d(3) + d(4).*100.^d(5));
z5=((d(2)./(d(1).*100) + d(2)) + (1 - (d(2)./(d(1).*100) + d(2))).*exp(-(d(1).*100 + d(2)).*r)) .* (d(3) + d(4).*100.^d(5));
X1=z5
Y1=ss
A4=(sum((X1-mean(Y1)).^2))./(sum((Y1-mean(Y1)).^2))
%h5 = plot(t, z_skm, '-', 'Color', COL_SKM, 'LineWidth', LINE_WIDTH_MAIN);

%% -------------------- Main figure styling (box off) --------------------
ax = gca;
set(ax, 'FontName', FONT_NAME, 'FontSize', FONT_SIZE_AX, 'LineWidth', 2.0); 
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.02, 0.02]); % Semi-open design

xlabel(' {\it t} (s)', 'FontSize', FONT_SIZE_LABEL);
ylabel('\tau (Pa)', 'FontSize', FONT_SIZE_LABEL);
xlim([0 60]); ylim([40 150]);

% Legend settings
lgd = legend([h0 h1 h3 h4], {'Test data', 'IFDM','VFDM','FDM'}, ...
    'Location','northwest', 'FontSize', FONT_SIZE_LEGEND, 'Box','off');

%% -------------------- Figure 2: independent plot of order trends --------------------
figure('Color','w','Position',[1050, 100, 900, 750]);
hold on; 

p1 = plot(t, alpha_nonlin_t, '-', 'Color', COL_IFDM, 'LineWidth', LINE_WIDTH_MAIN);
p2 = plot(t, alpha_linear_t, '--', 'Color', COL_VFDM, 'LineWidth', LINE_WIDTH_MAIN);
p3 = plot(t, alpha_fixed_0_1 * ones(size(t)), '-.', 'Color', COL_FDM, 'LineWidth', LINE_WIDTH_MAIN);

set(gca, 'FontName', FONT_NAME, 'FontSize', FONT_SIZE_AX, 'LineWidth', 2.0, ...
    'TickDir','out', 'Box','off');

xlabel('{\it t} (s)', 'FontSize', FONT_SIZE_LABEL);
ylabel(' \alpha({\it t})', 'FontSize', FONT_SIZE_LABEL);
xlim([0 60]); ylim([0, 1.7]);
legend([p1 p2 p3], {'IFDM','VFDM','FDM'}, 'Location','best', 'Box','off');

%% -------------------- Function definitions --------------------
function c_out = cfun(c, t)
    a = c(4) + (c(5)-c(4))./(1 + exp(-c(6).*(t - c(7))));
    c_out = c(1) + (c(2).*c(3).^a .* 0.1 .* t.^(1-a) ./ gamma(2-a));
end
function b_out = bfun(b, t)
    b_out = b(1) + (b(2).*b(3).^b(4) .* 0.1 .* t.^(1-b(4)) ./ gamma(2-b(4)));
end
function a_out = afun(a, t)
    al = a(4).*t + a(5);
    a_out = a(1) + (a(2).*a(3).^al .* 0.1 .* t.^(1-al) ./ gamma(2-al));
end
function d_out = dfun(d, t)
    d_out = ((d(2)./(d(1).*100) + d(2)) + (1 - (d(2)./(d(1).*100) + d(2))).*exp(-(d(1).*100 + d(2)).*t)) .* (d(3) + d(4).*100.^d(5));
end
function ae_out = aefun(a, t)
    ae_out =a(1)+(a(6)+(a(2)-a(6)).*((a(4).*t.^a(3)).^(a(5)./t))).*0.1
end
