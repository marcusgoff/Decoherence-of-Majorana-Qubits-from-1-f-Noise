%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Run_Figure3cd.m
% Purpose: Generate Figure 3c and 3d from the paper
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
%   - Generates Figure 3c and 3d from our paper.
%
% Dependencies:
%   - Results from Run_FigureS1.m must be run and saved in '../results'
%
% Requirements:
%   - MATLAB R2024 or newer
%   - Dependencies: functions in /QPP_Library directory. 
%
% Output:
%   - Saves Figure_3c and Figure_3d in the ./results
%   - Saves entire worskpace in the ./data directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

fprintf('\n Running Run_Figure3cd \n');


%% Load Results from Run_FigureS1
load('../results/figure_S1_data');

warning('off', 'all');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Extract the maxima


% Preallocate
Pqpp_max     = zeros(1, numel(P_qpp_final));
Pqpp_max_idx = zeros(1, numel(P_qpp_final));

for i = 1:numel(P_qpp_final)
    Pqpp_max(i)     = max(P_qpp_final{i});
    Pqpp_max_idx(i) = find(P_qpp_final{i} == Pqpp_max(i), 1);
end


delta_vec_GHz = delta_0_GHz*delta_factors;

Rqpp_max_MHz = 2*Pqpp_max(1:end)*1e3; % MHz

Rqpp_max_error_MHz = zeros(1, numel(err_bars));
for i = 1:numel(err_bars)
    Rqpp_max_error_MHz(i) = 2 * 1e3 * err_bars{i}(Pqpp_max_idx(i)); % MHz
end



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Plot T2^\ast Time - FIGURE 3C


T2_ns = 1e3./Rqpp_max_MHz;
T2_error_ns = T2_ns.*Rqpp_max_error_MHz./Rqpp_max_MHz;

% === Create figure ===
currFig = figure();
currFig.Position = [205 543 439 323]; %[205 543 495 323]; %[205   341   823   505];
ax = gca;

% Technically should divide by time below, but this is 1ns, so I haven't
% put it in. I've scalled by 1000 to get Rqpp in MHz. The 2 is for the two
% chains of the tetron.
% plot(delta_vec(1:end), Rqpp_max_MHz, 'o', ...
%     'Color', col_1, 'MarkerSize', marker_size, 'LineWidth', line_width*0.7, ...
%    'MarkerFaceColor', col_1, 'MarkerEdgeColor', 'w');

errorbar(delta_vec_GHz, T2_ns, T2_error_ns, 'o', ...
   'Color', col_mat(4,:), 'MarkerSize', marker_size, 'LineWidth', line_width*0.7, ...
   'MarkerFaceColor', col_mat(4,:), 'MarkerEdgeColor', 'w'); hold on;

errorbar(delta_vec_GHz(3), T2_ns(3), T2_error_ns(3), 'o', ...
   'Color', col_mat(3,:), 'MarkerSize', marker_size, 'LineWidth', line_width*0.7, ...
   'MarkerFaceColor', col_mat(3,:), 'MarkerEdgeColor', 'w'); hold on;

xlabel('\Delta/h (GHz)'); ylabel('T_2^\ast (ns)');

% === Log scale axes === 
set(ax, 'XScale', 'log');

% === Axis settings ===
ax.FontSize = fig_ax_font_size;
ax.LineWidth = fig_border_LW;
ax.TickDir = 'in';
ax.Box = 'off';
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';

ylim([0, 30]);
xlim([10^0, 10^3]);
ax.XTick = [1, 1e1, 1e2, 1e3];

% === Manually draw box ===
ax = gca;
xl = xlim(ax); yl = ylim(ax);
hold on;
plot([xl(1) xl(2)], [yl(1) yl(1)], 'k', 'LineWidth', ax.LineWidth); % bottom
plot([xl(1) xl(2)], [yl(2) yl(2)], 'k', 'LineWidth', ax.LineWidth); % top
plot([xl(1) xl(1)], [yl(1) yl(2)], 'k', 'LineWidth', ax.LineWidth); % left
plot([xl(2) xl(2)], [yl(1) yl(2)], 'k', 'LineWidth', ax.LineWidth); % right


%% Save Figure 3c
saveas(gcf, '../results/Figure_3c.svg');
saveas(gcf, '../results/Figure_3c.png');
saveas(gcf, '../results/Figure_3c.fig');


%-------------------------------------------------------------------------
%% -----------------------------------------------------------------------
% Plot Figure 3d Inset
% ------------------------------------------------------------------------

% T2 Time for experimental parameters in seconds:
T2_s = T2_ns(3)*1e-9; % seconds


x_point = log10((S_0)); %mu eV
y_point = log10(T2_s);

x_vec = linspace(-13+6, -4+6); 
log_y_int = y_point + 1*x_point; 

currFig = figure(); 
currFig.Position = [205 543 439 323]; 
line_1 = -x_vec + log_y_int;
loglog(10.^(x_vec), 10.^(line_1), 'LineWidth', line_width);


%xlabel('$S_0 ~(\mu\rm{eV})$', 'interpreter', 'latex');
%ylabel('$T_{2}^{\ast}~(\rm s)$', 'interpreter', 'latex');

xlabel('S_0 (\mueV)');
ylabel('T_2^\ast (s)');

ax = gca;
ax.FontSize = fig_ax_font_size;
ax.LineWidth = fig_border_LW;
ax.TickDir = 'in';
ax.Box = 'off';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

axis tight;
ylim([min(ylim), max(ylim)+1e-3]);


ax.YTick = 10.^([-9, -6, -3]);
ax.YAxis.MinorTickValues = 10.^([-8, -7, -5, -4, -3]);

ax.XTick = 10.^[-6, -3, -0 ];
ax.XAxis.MinorTickValues = 10.^[-13:0];

ylim( [3.9811e-10, 0.0033]);
xlim([10^(-5 - 0.3), 10^(0 + 0.3)]);


hold on;
ylim_current = ylim;
xlim_current = xlim;
% Bottom
loglog(xlim_current, [min(ylim_current) min(ylim_current)], 'k', 'LineWidth', ax.LineWidth);
% Left
loglog([min(xlim_current) min(xlim_current)], ylim_current, 'k', 'LineWidth', ax.LineWidth);
% Top
loglog(xlim_current, [max(ylim_current) max(ylim_current)], 'k', 'LineWidth', ax.LineWidth);
% Right
loglog([max(xlim_current) max(xlim_current)], ylim_current, 'k', 'LineWidth', ax.LineWidth);


hold on;
loglog(10.^(x_point), 10.^(y_point), '*', 'Color', col_mat(3,:), ...
    'MarkerSize', 13, 'LineWidth', line_width);



%% Save Figure 3d
saveas(gcf, '../results/Figure_3d.svg');
saveas(gcf, '../results/Figure_3d.png');
saveas(gcf, '../results/Figure_3d.fig');


%% Save workspace
close all;
save('../results/Figure_3cd_workspace');


