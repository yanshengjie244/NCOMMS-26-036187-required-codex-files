clear; clc; close all;

%% 1. Data definition
alpha_list = [0.2,0.4,0.6, 0.8,1.0]; 
t = 0.1:0.01:28;
sigma_ref = 1.0; E = 1.0; gamma_dot = 0.1;
a(4)=0.2; a(5)=2; a(6)=0.1
%% 2. Color assignment (based on the provided palette image)
% Colors are extracted from the uploaded palette and assigned to alpha_list in order.
colors = [
    133, 200, 221; % Light blue (alpha=2)
    212, 228, 173; % Light green (alpha=4)
    250, 224, 221; % Light pink (alpha=6)
    243, 178, 169; % Coral (alpha=8)
    166, 166, 210; % Light purple (alpha=10)
] / 255; % Normalize to [0, 1]

%% 3. Canvas settings (no background)
fig_width = 18.0; fig_height = 14.0;
figure('Units', 'centimeters', 'Position', [5, 5, fig_width, fig_height], 'Color', 'w');
hold on;

%% 4. Main plotting (without shadows or background layers)
for i = 1:length(alpha_list)
    % Calculate data
    z = sigma_ref + (E .* (gamma_dot .^alpha_list(i)).* t.^(1-alpha_list(i))) ./ gamma(2 - alpha_list(i) + eps);
    %z = sigma_ref+(t.^0.8.*0.8).^(alpha_list(i)./t)
    %z= sigma_ref+0.1.^(exp(-t.*0.2).*log(2./t)).*t.^alpha_list(i)
    % Keep only the main curves and increase line width for clarity.
    plot(t, z, 'Color', colors(i,:), 'LineWidth', 10.0); 
end

%% 5. Publication-style layout (Nature/Science style)
ax = gca;
set(ax, 'Layer', 'top', ...
        'FontName', 'Arial', ...
        'FontSize', 45, ... % Slightly reduced to avoid oversized text
        'LineWidth', 2.0, ...
        'Box', 'off', ...      % Remove the top and right borders
        'TickDir', 'out', ...
        'Color', 'none');      % Keep the axes background transparent/white

% Axis labels
xlabel('Time, {\it t} (s)', 'FontSize', 50);
ylabel('Shear stress, \tau (Pa)', 'FontSize', 50);

% Axis ranges and ticks
xlim([0, 15]); ylim([1, 5]);
set(ax, 'XTick', 0:5:15, 'YTick', 0:1:7);

%% 6. Colorbar (optional)
cb = colorbar;
colormap(ax, colors); 
set(cb, 'Units', 'normalized', 'Position', [0.88, 0.25, 0.02, 0.5], 'Box', 'off');
set(cb, 'Ticks', linspace(0.1, 0.9, length(alpha_list)), ...
        'TickLabels', {'0.2','0.4','0.6','0.8','1.0'}, ...
        'FontSize', 50, 'TickDirection', 'out');

% Final cleanup: remove the background grid and any remaining patches.
grid off;
