%----------------------------------------------------%
% Run all figures from a single instance of MATLAB
%----------------------------------------------------%

try
    Run_Figure2b;
    clearvars; 
    Run_Figure2c;
    clearvars; 
    Run_FigureS1;
    clearvars; 
    Run_Figure3cd;
    clearvars; 
    Run_FigureS2;
    clearvars; 
    Run_FigureS4;
catch ME
    disp(getReport(ME));
end