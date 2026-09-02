function Fig7_OV_Resistance_Model
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Written by: Khaphetsi Joseph Mahasa
	% Date: 03 Jul 2026
	% Function: Code to generate Figure 7 in the article.
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% 1. Document Configuration and Setup
    % Vary MOI (according to Liu et al (2022)) initial number of OVs in x0 
    % Figure 8 in the manuscript
    % Scenario 2: High infection rate (beta0 = 2.5e-9; % p3), strong immune response (eta = 0.5; % p12)
    
    global lambda_s K beta0 delta_c h_c gamma0 b omega phi h_m d_m eta h_d psi d_c 
    global x0 tStart tEnd tSpan TFtSpan

    beep off
    clc; close all
    set(0,'DefaultFigureWindowStyle','docked')
    set(0,'DefaultLineLineSmoothing','on');
    set(groot,'DefaultFigureWindowStyle','docked',...
        'defaultLineLineWidth',2,...
        'defaultLineMarkerSize',4,...
        'defaultAxesFontSize',18,...
        'defaultAxesFontWeight','bold',...
        'defaultAxesFontName','Times New Roman');

    mFileWorkingDirectory = pwd;
    mFileParentDirectory = mFileWorkingDirectory;
    
    format shortg
    clk = datestr(datetime('now'), 'dd-mmm-yy-HH');
    PlotsFolder = ['OV_Resistance_PlotsDirectory-', clk];
    
    if ~exist(fullfile(mFileParentDirectory, PlotsFolder), 'dir')
      mkdir(mFileParentDirectory, PlotsFolder);
    end
    
    ResultsFolder = ['OV_Resistance_ResultsDirectory-', clk];
    if ~exist(fullfile(mFileParentDirectory, ResultsFolder), 'dir')
      mkdir(mFileParentDirectory, ResultsFolder);
    end

    PlotsFullPath = fullfile(mFileParentDirectory, PlotsFolder);
    ResultsFullPath = fullfile(mFileParentDirectory, ResultsFolder);
    diary(fullfile(ResultsFullPath, 'OV_Resistance_Results.txt'));
    
%% 2. Define Parameter Set
    lambda_s = 0.075; % p1
    K = 5.3196e8; % p2 
    beta0 = 2.5e-9; % p3
    delta_c = 2.68e-9; % p4
    h_c = 40; % p5
    gamma0 = 0.05; % p6
    b = 3000; % p7
    omega = 1.2; % p8
    phi =0.7e2; % p9
    h_m = 2e7; % p10
    d_m = 0.25; % p11
    eta = 0.5; % p12
    h_d = 2.019e7; % p13
    psi = 0.005; % p14
    d_c = 0.02; % p15 
   
%% 3. Primary Model Solution Configurations
    tStart = 0;
    tEnd = 100;
    MOI = [0.001, 0.01, 0.1, 1, 10];
    PFU = [20, 200, 2000, 20000, 200000];
    tSpan = [0 300]; 
    TFtSpan = tStart:180;

%% 4. Parametric Iteration Loop & Plot Generation  
    cd(PlotsFullPath);
    Styles = {'-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--', '-.', ':'}; 
    options = odeset('Refine', 10, 'RelTol', 1e-4);

    % Run through each experimental PFU initial constraint tier
    for i = 1:length(PFU)
        % Map updated vector conditions dynamically matching each trial segment
        x0 = [2e6, 0, PFU(i), 0, 1e6];  
        
        % Solve using ode23s tracking the anonymous function to map global equations safely
        [t, x] = ode23s(@(t,x) ModelEquations(t,x), tSpan, x0, options);   
        
        FUnT = x(:,1);
        FInfT = x(:,2);
        FV = x(:,3);
        FM = x(:,4);
        FCD = x(:,5);

        % Initialize an independent, isolated figure shell window for the current index run
        figure(i)
        hold on; box on;
        plot(t, FUnT + FInfT, 'LineStyle', Styles{1}, 'LineWidth', 2.5)  % Total tumor
        plot(t, FCD, 'LineStyle', Styles{4}, 'LineWidth', 2.5)          % Immune cells
        plot([0, 300], [1e6, 1e6], 'k--', 'LineWidth', 1.5)             % Threshold line

        % Coordinate layouts layout management rules
        ylim([0 3.5e7]) 
        xlim([tSpan(1) tSpan(2)])
        legend('T(t)+T_{i}(t)', 'C_{d}(t)', 'Threshold', 'location', 'SouthEast', 'Interpreter', 'latex')
        xlabel('Time (days)', 'Interpreter', 'latex')
        ylabel('Population size (\# cells)', 'Interpreter', 'latex')
        title(sprintf('MOI: %g', MOI(i)), 'Interpreter', 'latex')
             
        %%% Inset Figure 1: Free Oncolytic Virus Profile %%%
        axes('position', [.315 .55 .25 .28], 'NextPlot', 'add')
        box on 
        indexOfInterest = (t <= 200) & (t >= 0); 
        plot(t(indexOfInterest), FV(indexOfInterest), 'LineStyle', Styles{2}, 'Color', 'r', 'LineWidth', 2)
        legend('V(t)', 'Interpreter', 'latex', 'Location', 'best')
        ylabel('PFU/cell', 'Interpreter', 'latex')
        xlabel('Time (days)', 'Interpreter', 'latex')
        set(gca, 'YScale', 'log')
        axis tight

        %%% Inset Figure 2: Resistance Cell Profile %%%
        axes('position', [.65 .55 .25 .28], 'NextPlot', 'add')
        box on 
        plot(t(indexOfInterest), FM(indexOfInterest), 'LineStyle', Styles{3}, 'Color', 'g', 'LineWidth', 2)
        legend('M(t)', 'Interpreter', 'latex', 'Location', 'best')
        ylabel('\# cells', 'Interpreter', 'latex')
        xlabel('Time (days)', 'Interpreter', 'latex')
        set(gca, 'YScale', 'log')
        axis tight
        hold off
        
        % Normalize syntax names replacing float dots to prevent compilation save failures
        filename = sprintf('Figure7_MOI_%g.png', MOI(i));
        filename = strrep(filename, '.', 'p'); 
        filename = strrep(filename, 'png', '.png'); 
        saveas(gcf, filename);
    end 
       
    diary off
    cd(mFileWorkingDirectory);
end 

%% ========================================================================
%  SUBFUNCTIONS SECTION
%% ========================================================================

function dx = ModelEquations(~, x) 
    global lambda_s K beta0 delta_c h_c gamma0 b omega phi h_m d_m eta h_d psi d_c 
    dx = zeros(5,1);  
    
    dx(1) = lambda_s*x(1)*(1 - x(1)/K) - beta0*x(1)*x(3) - (delta_c*x(1)*x(5))/(h_c + x(5));
    dx(2) = beta0*x(1)*x(3) - (delta_c*x(2)*x(5))/(h_c + x(5)) - gamma0*x(2);
    dx(3) = b*gamma0*x(2) - omega*x(3);
    dx(4) = (phi*x(2))/(h_m + x(2)) - d_m*x(4); 
    dx(5) = ((eta*(x(1) + x(2)))/(h_d + x(1) + x(2)))*(1 - psi*x(4))*x(5) - d_c*x(5);
end

function [value, isterminal, direction] = myEvents(t, x)
    DerivSol = ModelEquations(t, x);
    value = [DerivSol(1) + 1e-01, DerivSol(2) + 1e-01];
    isterminal = [0, 0];
    direction = [-1, -1];
end
