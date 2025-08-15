%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Run_All_Figures.m
% Purpose: Computes all data for and plots all numerics results presented in
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
%   Runs the full code for all numerical results presented in our paper and 
%   plots the figures.
% Expected Run Time:
%   >72 hours (tested on a Macbook pro). 
%
% Requirements:
%   - MATLAB R2024 or newer
%   - Dependencies: functions in /QPP_Library directory. 
%
% Output:
%   - Saves all figures as .png, .svg, and .fig files in the ./results
%   - Saves all data corresponding to each figure in the ./data directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note to the Reviewers:
% You may comment out any of the following lines to exclude running that
% figure's matlab script. To comment out a line simply type a "%" character 
% at the start of that line. 

try
    run('Code/Run_Figure2b.m');
    clearvars;
    run('Code/Run_Figure2c.m');
    clearvars;
    run('Code/Run_FigureS1.m');
    clearvars;
    run('Code/Run_Figure3cd.m');
    clearvars;
    run('Code/Run_FigureS2.m');
    clearvars;
    run('Code/Run_FigureS4.m');
catch ME
    disp(getReport(ME));
end


