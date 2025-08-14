%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Run_FigureS2.m
% Purpose: Generate Figure S2 a,b,c from the Supplementary Information of
%          "Decoherence in Majorana Qubits by 1/f Noise"
%
% Author: Marcus C. Goffage
% Date: 12-Aug-2025
% Affiliation: University of New South Wales
%
% Paper: "Decoherence in Majorana Qubits by 1/f Noise"
% Paper Authors: A. Alase^1, M. C. Goffage^2, M. C. Cassidy^2, 
%                S. N. Coppersmith^{2*}
% Affiliations:  ^1 University of Sydney
%                ^2 University of New South Wales
%                *  Corresponding Author
%
% -------------------------------------------------------------------------
% ABOUT THIS SCRIPT
% -------------------------------------------------------------------------
% Description:
%   - Generates Figure S2a,b,c from the supplementary information of our 
%     paper.
%
% Dependencies:
%   - Results from Run_Figure2c.m must be run and saved in '../results'
%
% Requirements:
%   - MATLAB R2024 or newer
%   - Dependencies: functions in /QPP_Library directory. 
%
% Output:
%   - Saves Figure_S2a, Figure_S2b and Figure_S2c in the ../results
%   - Saves entire worskpace in the ../results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

fprintf('\n Running Run_FigureS2 \n');

clear all;
%% Load Results from Run_Figure2c
load('../results/figure_2c_data');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Plot Gamma Max Results

% === Use same styling params from first plot ===
marker_size = 17;
line_width = 4;
RGB = orderedcolors("gem");

col_1 = RGB(1,:);  % 3 μm color
col_2 = RGB(2,:);  % 5 μm color
col_3 = RGB(3,:);  % 10 μm color

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Figure Params

% === Use same styling params from first plot ===
marker_size = 17;
line_width = 4;
fig_ax_font_size = 40;
RGB = orderedcolors("gem");

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Get Perfectly Incoherent Transition Lines
% Under assumption that the wavefunction dephases entirely between TLF
% transitions

delta_mu = (mu_high - mu_low);

L_even_jump_a = abs((delta_mu.^2 ./ (delta .* w)) .* (N_a) ./ 32);
P_qpp_inco_a  = t_final * gamma_vec_a * L_even_jump_a * 2 / 2;

L_even_jump_b = abs((delta_mu.^2 ./ (delta .* w)) .* (N_b) ./ 32);
P_qpp_inco_b  = t_final * gamma_vec_b * L_even_jump_b * 2 / 2;

L_even_jump_c = abs((delta_mu.^2 ./ (delta .* w)) .* (N_c) ./ 32);
P_qpp_inco_c  = t_final * gamma_vec_c * L_even_jump_c * 2 / 2;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%
% === Create figure ===
currFig = figure();
currFig.Position = [205 190 915 571];
ax = gca;
hold on;

% === Plot error bars in requested color mapping ===
h_3um_num = errorbar(gamma_vec_ghz, P_qpp_final_a, err_bars_a, '-.o', ...
   'Color', col_3, 'MarkerSize', marker_size, 'LineWidth', line_width, ...
   'MarkerFaceColor', col_3, 'MarkerEdgeColor', 'w');

h_5um_num = errorbar(gamma_vec_ghz, P_qpp_final_b, err_bars_b, '-.o', ...
   'Color', col_2, 'MarkerSize', marker_size, 'LineWidth', line_width, ...
   'MarkerFaceColor', col_2, 'MarkerEdgeColor', 'w');

h_10um_num = errorbar(gamma_vec_ghz, P_qpp_final_c, err_bars_c, '-.o', ...
   'Color', col_1, 'MarkerSize', marker_size, 'LineWidth', line_width, ...
   'MarkerFaceColor', col_1, 'MarkerEdgeColor', 'w');

% === Plot Completely Incoherent Theory Result === %
blend_factor = 0.3;
col_1_light = col_1 + (1 - col_1) * blend_factor;
col_2_light = col_2 + (1 - col_2) * blend_factor;
col_3_light = col_3 + (1 - col_3) * blend_factor;
h_3um_inco = plot(gamma_vec_ghz, P_qpp_inco_a, ':', ...
   'Color', col_3_light, 'MarkerSize', marker_size, 'LineWidth', line_width, ...
   'MarkerFaceColor', col_3_light, 'MarkerEdgeColor', 'w');

h_5um_inco = plot(gamma_vec_ghz, P_qpp_inco_b, ':', ...
   'Color', col_2_light, 'MarkerSize', marker_size, 'LineWidth', line_width, ...
   'MarkerFaceColor', col_2_light, 'MarkerEdgeColor', 'w');

h_10um_inco = plot(gamma_vec_ghz, P_qpp_inco_c, ':', ...
   'Color', col_1_light, 'MarkerSize', marker_size, 'LineWidth', line_width, ...
   'MarkerFaceColor', col_1_light, 'MarkerEdgeColor', 'w');

% === Log scale axes ===
set(ax, 'XScale', 'log');
set(ax, 'YScale', 'log');

% === Axis tweaks matching previous plot ===
ax.FontSize = fig_ax_font_size;
ax.LineWidth = 2.3;  % same as fig_border_LW in first plot
ax.TickDir = 'in';
ax.Box = 'off';
ax.XMinorTick = 'on';
ax.YMinorTick = 'off';

% === Axis limits and labels ===
axis tight;
xlim([9.99, 1.01e3]);
xlabel('Two level fluctuator switching rate, \Gamma (GHz)');
%ylabel({'Average Number of', 'Quasiparticle Pairs'}, 'HorizontalAlignment', 'center');
%ylabel('$P_{\rm QPP}$', 'interpreter', 'latex');
ylabel('P_{QPP,1ns}');

ylim_curr = ylim; ylim_curr(2) = 1;
ylim(ylim_curr);

xlim([2*9, 2*1.1e3]);

xl = xlim(ax); yl = ylim(ax);

plot([xl(1) xl(2)], [yl(1) yl(1)], 'k', 'LineWidth', ax.LineWidth) % bottom
plot([xl(1) xl(2)], [yl(2) yl(2)], 'k', 'LineWidth', ax.LineWidth) % top
plot([xl(1) xl(1)], [yl(1) yl(2)], 'k', 'LineWidth', ax.LineWidth) % left
plot([xl(2) xl(2)], [yl(1) yl(2)], 'k', 'LineWidth', ax.LineWidth) % right

% === Legend with requested order ===
hLeg = legend([h_10um_num, h_5um_num, h_3um_num, h_10um_inco, h_5um_inco, h_3um_inco], ...
    {sprintf('$\\mathcal{L} = %.3g \\, \\mu \\rm{m}$', round(1e6*N_c*a,1)), ...
     sprintf('$\\mathcal{L} = %.3g \\, \\mu \\rm{m}$', round(1e6*N_b*a,1)), ...
     sprintf('$\\mathcal{L} = %.3g \\, \\mu \\rm{m}$', round(1e6*N_a*a,1)),...
     sprintf('$\\mathcal{L} = %.3g \\, \\mu \\rm{m}$', round(1e6*N_c*a,1)), ...
     sprintf('$\\mathcal{L} = %.3g \\, \\mu \\rm{m}$', round(1e6*N_b*a,1)), ...
     sprintf('$\\mathcal{L} = %.3g \\, \\mu \\rm{m}$', round(1e6*N_a*a,1))}, ...
    'Interpreter', 'latex', ...
    'Location', 'northwest', ...
    'Box', 'off');

hLeg = legend([h_10um_num, h_5um_num, h_3um_num], ...
    {sprintf('$\\mathcal{L} = %.3g \\, \\mu \\rm{m}$', round(1e6*N_c*a,1)), ...
     sprintf('$\\mathcal{L} = %.3g \\, \\mu \\rm{m}$', round(1e6*N_b*a,1)), ...
     sprintf('$\\mathcal{L} = %.3g \\, \\mu \\rm{m}$', round(1e6*N_a*a,1))}, ...
    'Interpreter', 'latex', ...
    'Location', 'northwest', ...
    'Box', 'off');

%hLeg.Position = [0.6770 0.1677 0.1768 0.1790];
hLeg.Location = 'northwest';
hLeg.FontSize = 33;
%hLeg.Position = [205 190 915 571];
hold off;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Save Figure S2a

% Save as MATLAB .fig image
saveas(gcf, '../results/figure_S2a.fig');

% Save as PNG image
saveas(gcf, '../results/figure_S2a.png');
saveas(gcf, '../results/figure_S2a.svg');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Get F Factors

f_factor_a = P_qpp_final_a(:)./P_qpp_inco_a(:);
f_factor_b = P_qpp_final_b(:)./P_qpp_inco_b(:);
f_factor_c = P_qpp_final_c(:)./P_qpp_inco_c(:);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Plot F Factor Linear Scale

% === Create figure ===
currFig = figure();
currFig.Position = [205 190 915 571];
ax = gca;
hold on;
% === Plot Completely Incoherent Theory Result === %
h_3um_fac = plot(gamma_vec_ghz, f_factor_a, 'o-.', ...
   'Color', col_3*0.8, 'MarkerSize', marker_size, 'LineWidth', line_width, ...
   'MarkerFaceColor', col_3, 'MarkerEdgeColor', 'w');

h_5um_fac = plot(gamma_vec_ghz, f_factor_b, 'o-.', ...
   'Color', col_2*0.8, 'MarkerSize', marker_size, 'LineWidth', line_width, ...
   'MarkerFaceColor', col_2, 'MarkerEdgeColor', 'w');

h_10um_fac = plot(gamma_vec_ghz, f_factor_c, 'o-.', ...
   'Color', col_1*0.8, 'MarkerSize', marker_size, 'LineWidth', line_width, ...
   'MarkerFaceColor', col_1*0.8, 'MarkerEdgeColor', 'w');


% === Log scale axes ===
set(ax, 'XScale', 'log');
 %set(ax, 'YScale', 'log');

% === Axis tweaks matching previous plot ===
ax.FontSize = fig_ax_font_size;
ax.LineWidth = 2.3;  % same as fig_border_LW in first plot
ax.TickDir = 'in';
ax.Box = 'off';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

ax.YAxis.MinorTickValues = 0.1:0.1:1;

% === Axis limits and labels ===
axis tight;

xlabel('Frequency (GHz)');
xlabel('Two level fluctuator switching rate, \Gamma (GHz)');
ylabel('$\mathcal{F} = P_{\rm{QPP,1ns}}^{\rm{num}}/P_{\rm{QPP,1ns}}^{\rm{inco}}$', 'interpreter', 'latex');

xlim([2*9, 2*1.1e3]);

yticks([0:0.2:1]);

ylim([0, 1]);

xl = xlim(ax); yl = ylim(ax);

plot([xl(1) xl(2)], [yl(1) yl(1)], 'k', 'LineWidth', ax.LineWidth) % bottom
plot([xl(1) xl(2)], [yl(2) yl(2)], 'k', 'LineWidth', ax.LineWidth) % top
plot([xl(1) xl(1)], [yl(1) yl(2)], 'k', 'LineWidth', ax.LineWidth) % left
plot([xl(2) xl(2)], [yl(1) yl(2)], 'k', 'LineWidth', ax.LineWidth) % right


% === Legend with requested order ===
hLeg = legend([h_10um_fac, h_5um_fac, h_3um_fac], ...
    {sprintf('$\\mathcal{L} = %.3g \\, \\mu \\rm{m}$', round(1e6*N_c*a,1)), ...
     sprintf('$\\mathcal{L} = %.3g \\, \\mu \\rm{m}$', round(1e6*N_b*a,1)), ...
     sprintf('$\\mathcal{L} = %.3g \\, \\mu \\rm{m}$', round(1e6*N_a*a,1))}, ...
    'Interpreter', 'latex', ...
    'Location', 'southeast', ...
    'Box', 'off');

hLeg.Location = 'northeast';
hold off;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Save Figure S2b

% Save as MATLAB .fig image
saveas(gcf, '../results/figure_S2b.fig');

% Save as PNG image
saveas(gcf, '../results/figure_S2b.png');
saveas(gcf, '../results/figure_S2b.svg');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Call Constants
delta_GHz = 10^(-9)*delta_muev/(2*pi*hbar_mueVs);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% f Factor Linear Scale
delta_GHz = 10^(-9)*delta_muev/(2*pi*hbar_mueVs);

% === Create figure ===
currFig = figure();
currFig.Position = [205   190   938   571];
ax = gca;
hold on;
% === Plot Completely Incoherent Theory Result === %
h_3um_fac = plot(gamma_vec_ghz, f_factor_a(:).*gamma_vec_ghz(:)/(4*delta_GHz), 'o-.', ...
   'Color', col_3*0.8, 'MarkerSize', marker_size, 'LineWidth', line_width, ...
   'MarkerFaceColor', col_3, 'MarkerEdgeColor', 'w');

h_5um_fac = plot(gamma_vec_ghz, f_factor_b(:).*gamma_vec_ghz(:)/(4*delta_GHz), 'o-.', ...
   'Color', col_2*0.8, 'MarkerSize', marker_size, 'LineWidth', line_width, ...
   'MarkerFaceColor', col_2, 'MarkerEdgeColor', 'w');

h_10um_fac = plot(gamma_vec_ghz, f_factor_c(:).*gamma_vec_ghz(:)/(4*delta_GHz), 'o-.', ...
   'Color', col_1*0.8, 'MarkerSize', marker_size, 'LineWidth', line_width, ...
   'MarkerFaceColor', col_1*0.8, 'MarkerEdgeColor', 'w');


% === Log scale axes ===
set(ax, 'XScale', 'log');
 %set(ax, 'YScale', 'log');

% === Axis tweaks matching previous plot ===
ax.FontSize = fig_ax_font_size;
ax.LineWidth = 2.3;  % same as fig_border_LW in first plot
ax.TickDir = 'in';
ax.Box = 'off';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

ax.YAxis.MinorTickValues = 0.1:0.1:1;

% === Axis limits and labels ===
axis tight;


%ylabel({'Average Number of', 'Quasiparticle Pairs'}, 'HorizontalAlignment', 'center');
xlabel('Two level fluctuator switching rate, \Gamma (GHz)');
%ylabel('$\mathcal{F}_{\Gamma} = P_{\rm{QPP, num}}/P_{\rm{QPP, inco}}$', 'interpreter', 'latex');
ylabel('$\mathcal{F} \Gamma/(4 \Delta/h)$', 'interpreter', 'latex');

xlim([2*9, 2*1.1e3]);
ylim([0, 0.8]);

xl = xlim(ax); yl = ylim(ax);
plot([xl(1) xl(2)], [yl(1) yl(1)], 'k', 'LineWidth', ax.LineWidth) % bottom
plot([xl(1) xl(2)], [yl(2) yl(2)], 'k', 'LineWidth', ax.LineWidth) % top
plot([xl(1) xl(1)], [yl(1) yl(2)], 'k', 'LineWidth', ax.LineWidth) % left
plot([xl(2) xl(2)], [yl(1) yl(2)], 'k', 'LineWidth', ax.LineWidth) % right


% === Legend with requested order ===
hLeg = legend([h_10um_fac, h_5um_fac, h_3um_fac], ...
    {sprintf('$\\mathcal{L} = %.3g \\, \\mu \\rm{m}$', round(1e6*N_c*a,1)), ...
     sprintf('$\\mathcal{L} = %.3g \\, \\mu \\rm{m}$', round(1e6*N_b*a,1)), ...
     sprintf('$\\mathcal{L} = %.3g \\, \\mu \\rm{m}$', round(1e6*N_a*a,1))}, ...
    'Interpreter', 'latex', ...
    'Location', 'southeast', ...
    'Box', 'off');

hLeg.Location = 'northeast';
hold off;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Save Figure S2c

% Save as MATLAB .fig image
saveas(gcf, '../results/figure_S2c.fig');

% Save as PNG image
saveas(gcf, '../results/figure_S2c.png');
saveas(gcf, '../results/figure_S2c.svg');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Save Data
close all;
save('../results/figure_S2_data');