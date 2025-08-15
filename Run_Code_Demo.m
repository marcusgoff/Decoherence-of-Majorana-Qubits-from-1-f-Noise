%----------------------------------------------------%
% Run Code Demo
%----------------------------------------------------%
% Description:
%   - Demo code which runs in 4 minutes (tested on a Macbook pro). Calculates and plots P_qpp
%     as a function of time for L = 3 microns, and for a TLF rate of 200
%     GHz, for 1 nanosecond (this corresponds to the yellow line in figure
%     2b in the main text). 
%
% Requirements:
%   - MATLAB R2024 or newer
%   - Dependencies: functions in /QPP_Library directory. 
%   - Run_Figure_2c must be run before this script
%
% Output:
%   - Saves Figure_2b in the ./results
%   - Saves entire worskpace in the ./data directory
%
%----------------------------------------------------%


try
    run('Code/code_demo.m');
catch ME
    disp(getReport(ME));
end


