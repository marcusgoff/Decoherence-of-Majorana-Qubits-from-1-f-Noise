%----------------------------------------------------%
% Run all figures 
%----------------------------------------------------%

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


