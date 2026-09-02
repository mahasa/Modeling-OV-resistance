function Fig6_RelationshipBTW_R0_TCP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Khaphetsi Joseph Mahasa
% Date: 03 Jul 2026
% Function: Code to generate Figure 6 in the article.
%
% Description: This script simulates and plots the multi-panel dynamics 
%              of tumor cell populations (uninfected vs. infected) alongside 
%              R_0, A, and total tumor burden alongside Tumor Control Probability (TCP). 
%              It loops through low and high viral infection rates (\beta_0) 
%              and varies tumor proliferation rates (\lambda_s) over a 600-day 
%              timeline, applying automatic log-scaling for high infection rates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic 
    global lambda_s K beta0 delta_c h_c gamma0 b omega phi h_m d_m eta h_d psi d_c 
    global x0 tStart tEnd tSpan TFtSpan PlotsFullPath ResultsFullPath

    beep off
    clc; close all
    
    % Set publication-grade graphics style
    set(groot, 'DefaultFigureWindowStyle', 'docked', ...
        'defaultLineLineWidth', 2, ...
        'defaultLineMarkerSize', 4, ...
        'defaultAxesFontSize', 18, ...
        'defaultAxesFontWeight', 'bold', ...
        'defaultAxesFontName', 'Times New Roman');

    mFileWorkingDirectory = pwd;
    mFileParentDirectory = mFileWorkingDirectory;
    
    % Generate secure timestamp format string for directory naming
    clk = datestr(datetime('now'), 'dd-mmm-yy-HH');
    PlotsFolder = ['\OV_Resistance_PlotsDirectory-', clk];
    ResultsFolder = ['\OV_Resistance_ResultsDirectory-', clk];
    
    % Create paths if they do not exist
    PlotsFullPath = [mFileParentDirectory, PlotsFolder];
    ResultsFullPath = [mFileParentDirectory, ResultsFolder];
    
    if ~exist(PlotsFullPath, 'dir')
        mkdir(PlotsFullPath);
    end
    if ~exist(ResultsFullPath, 'dir')
        mkdir(ResultsFullPath);
    end
    
    diary([ResultsFullPath, '\OV_Resistance_Results.txt']);
    
    %% Define parameters here
    lambda_s = 0.075;   % p1: Tumor cell growth rate
    K        = 10^9;    % p2: Tumor carrying capacity
    delta_c  = 2.68e-9; % p4: T-cell killing rate
    h_c      = 40;      % p5: T-cell influx / recruitment constant
    gamma0   = 0.06;    % p6: Virus clearance rate
    b        = 3000;    % p7: Viral burst size
    omega    = 1.2;     % p8: Infected cell death rate
    phi      = 0.7e2;   % p9: Recruitment rate
    h_m      = 2e7;     % p10: Half-saturation constant for TAMC
    d_m      = 0.25;    % p11: Clearance rate of TAMC
    eta      = 0.1;     % p12: T-cell proliferation rate
    h_d      = 2.019e7; % p13: Half-saturation constant for T-cell proliferation
    psi      = 0.005;   % p14: T-cell suppression constant
    d_c      = 0.02;    % p15: T-cell natural death rate
    
    %% Primary model solution setup
    tStart = 0;
    tEnd = 100;
    x0 = [2e7, 1, 2000, 0, 1e6];   % Initial conditions 
    tSpan = [0 600];               % Time span of the simulation
    TFtSpan = tStart:180;          % Target evaluation grid
    
    %% Array Configuration for Parameter Variations
    beta0_vec = [2.3e-15, 2.3e-12]; % Target infection rate variations
    lambdaSRate = [0.065, 0.070, 0.075, 0.085, 0.095]; 
    Styles = {'-', '--', '-.', ':', '-', '--'};
    
    % Safely change directory for plot saving destination
    if exist(PlotsFullPath, 'dir')
        cd(PlotsFullPath);
    end
    
    %% =========================================================================
    % MAIN PARAMETRIC LOOP OVER BETA0
    % =========================================================================
    for bIdx = 1:length(beta0_vec)
        beta0 = beta0_vec(bIdx); % Update global variable for ODE visibility
        
        % Dynamic Tag strings for identification file naming
        betaTag = num2str(beta0, '%10.0e');
        betaTagClean = strrep(strrep(betaTag, ' ', ''), '+', ''); % Strip spaces and plus signs
        
        % ---------------------------------------------------------------------
        % FIGURE 1: Dynamic Populations (Uninfected vs Infected Tumor Size)
        % ---------------------------------------------------------------------
        fig1 = figure('Name', ['Tumor Dynamics Beta = ', betaTag], 'NumberTitle', 'off'); 
        
        % --- Subplot 1: Uninfected Tumor ---
        subplot(1,2,1)
        legendInfo1 = cell(1, length(lambdaSRate));
        h1 = zeros(1, length(lambdaSRate));
        
        for k = 1:length(lambdaSRate)
            lambda_s = lambdaSRate(k);
            eta_val = 0.1; 
            
            SR = (eta_val * delta_c) / (lambda_s * d_c); 
            ttle = ['\beta = ', betaTag];
            R0 = (beta0 * b * K) / omega; 
            
            options = odeset('NonNegative', 1:5);          
            [t1, x1] = ode23s(@(t,x) ModelEquations(t,x,lambda_s,eta_val), tSpan, x0, options);
            
            T = x1(:,1); 
            
            hold on;
            grid on;
            h1(k) = plot(t1, T, 'LineStyle', Styles{k}, 'LineWidth', 2.5); 
            
            legendInfo1{k} = ['\lambda_s = ', num2str(lambdaSRate(k)), ', A = ', num2str(SR, '%.2e'), ', R_0 = ', num2str(R0, '%.2f')];
        end
        ylabel('Uninfected tumor (log(\# cells))', 'Interpreter', 'latex');
        xlabel('Time (days)', 'Interpreter', 'latex');
        set(gca, 'YScale', 'log'); % Always log scale for Uninfected cells
        title(ttle);
        legend(h1, legendInfo1, 'location', 'NorthEast', 'FontSize', 11);
        legend boxoff;
        hold off;
        
        % --- Subplot 2: Infected Tumor ---
        subplot(1,2,2)
        legendInfo2 = cell(1, length(lambdaSRate));
        hp = zeros(1, length(lambdaSRate));
        
        for k = 1:length(lambdaSRate)
            lambda_s = lambdaSRate(k);
            eta_val = 0.1;
            
            SR = (eta_val * delta_c) / (lambda_s * d_c); 
            
            options = odeset('NonNegative', 1:5);          
            [t1, x1] = ode23s(@(t,x) ModelEquations(t,x,lambda_s,eta_val), tSpan, x0, options);
            
            Ti = x1(:,2); 
            
            hold on;
            grid on;
            hp(k) = plot(t1, Ti, 'LineStyle', Styles{k}, 'LineWidth', 2.5); 
            
            legendInfo2{k} = ['\lambda_s = ', num2str(lambdaSRate(k)), ', A = ', num2str(SR, '%.2e')];
        end
        plot([0,600],[1e0,1e0],'b--') % Baseline threshold line
        
        % Apply log scale conditionally based on requested configuration criteria
        if beta0 == 2.3e-12
            ylabel('Infected tumor size (log(\# cells))', 'Interpreter', 'latex');
            set(gca, 'YScale', 'log');
        else
            ylabel('Infected tumor size (\# cells)', 'Interpreter', 'latex');
        end
        
        xlabel('Time (days)', 'Interpreter', 'latex');
        title(['\beta = ', betaTag]);
        legend(hp, legendInfo2, 'location', 'SouthEast', 'FontSize', 11);
        legend boxoff;
        hold off;
    
        % Save Figure 1 safely per loop condition iteration
        fig1Name = ['Tumor_Dynamics_Breakdown_Beta_', betaTagClean, '.png'];
        print(fig1, fullfile(PlotsFullPath, fig1Name), '-dpng', '-r300');
    
        % ---------------------------------------------------------------------
        % FIGURE 2: Total Tumor Volume vs. Tumor Control Probability (TCP)
        % ---------------------------------------------------------------------
        fig2 = figure('Name', ['Total Tumor vs TCP Beta = ', betaTag], 'NumberTitle', 'off'); 
    
        % --- Subplot 1: Total Tumor (Uninfected + Infected) ---
        subplot(1,2,1) 
        legendInfo3 = cell(1, length(lambdaSRate));
        h3 = zeros(1, length(lambdaSRate));
    
        for k = 1:length(lambdaSRate)
            lambda_s = lambdaSRate(k);
            eta_val = 0.1;
            
            SR = (eta_val * delta_c) / (lambda_s * d_c);
            
            options = odeset('NonNegative', 1:5);          
            [t1, x1] = ode23s(@(t,x) ModelEquations(t,x,lambda_s,eta_val), tSpan, x0, options);
            
            T = x1(:,1); 
            Ti = x1(:,2); 
            
            hold on;
            grid on;
            h3(k) = plot(t1, T + Ti, 'LineStyle', Styles{k}, 'LineWidth', 2.5); 
            
            legendInfo3{k} = ['\lambda_s = ', num2str(lambdaSRate(k)), ',  A = ', num2str(SR, '%.2e')];
        end
        ylabel('Total tumor size (\# cells)', 'Interpreter', 'latex');
        xlabel('Time (days)', 'Interpreter', 'latex');
        title(['\beta = ', betaTag]);
        legend(h3, legendInfo3, 'location', 'SouthEast', 'FontSize', 11);
        legend boxoff;
        hold off;
    
        % --- Subplot 2: Tumor Control Probability (TCP) ---
        subplot(1,2,2) 
        legendInfo4 = cell(1, length(lambdaSRate));
        h0 = zeros(1, length(lambdaSRate));
    
        for k = 1:length(lambdaSRate)
            lambda_s = lambdaSRate(k);
            eta_val = 0.1;
            
            SR = (eta_val * delta_c) / (lambda_s * d_c); 
            
            options = odeset('NonNegative', 1:5);          
            [t1, x1] = ode23s(@(t,x) ModelEquations(t,x,lambda_s,eta_val), tSpan, x0, options);
            
            T = x1(:,1); 
            
            % Calculate TCP (Poisson formula based on remaining uninfected tumor burden)
            TCP = exp(-T);
            
            hold on;
            grid on;
            h0(k) = plot(t1, TCP, 'LineStyle', Styles{k}, 'LineWidth', 2.5); 
            
            legendInfo4{k} = ['\lambda_s = ', num2str(lambdaSRate(k)), ',  A = ', num2str(SR, '%.2e')];
        end
        ylabel('Tumor control probability', 'Interpreter', 'latex');
        xlabel('Time (days)', 'Interpreter', 'latex');
        ylim([0 1]);
        legend(h0, legendInfo4, 'location', 'SouthEast', 'FontSize', 11);
        legend boxoff;
        hold off;
    
        % Save Figure 2 using your requested print specification pattern
        fig2Name = ['BurstInfection_Optimal_OV_Dose_Beta_', betaTagClean, '.png'];
        print(fig2, fullfile(PlotsFullPath, fig2Name), '-dpng', '-r300');
        
    end % End of beta0 loop parameter evaluation block

    %% Clean up logging resources and return to workspace root safely
    diary off;
    toc;
    cd(mFileWorkingDirectory);
end


%% =========================================================================
% MODEL EQUATIONS SUBFUNCTION
% =========================================================================
function dx = ModelEquations(~, x, lambda_s, eta) 
    global K beta0 delta_c h_c gamma0 b omega phi h_m d_m h_d psi d_c
    
    % Preallocate tracking space as a column vector
    dx = zeros(5,1); 
    
    % Core Differential Equations System
    dx(1) = lambda_s*x(1)*(1 - x(1)/K) - beta0*x(1)*x(3) - (delta_c*x(1)*x(5))/(h_c + x(5));
    dx(2) = beta0*x(1)*x(3) - (delta_c*x(2)*x(5))/(h_c + x(5)) - gamma0*x(2);
    dx(3) = b*gamma0*x(2) - omega*x(3);
    dx(4) = (phi*x(2))/(h_m + x(2)) - d_m*x(4);
    dx(5) = ((eta*(x(1) + x(2)))/(h_d + x(1) + x(2)))*(1 - psi*x(4))*x(5) - d_c*x(5);
end

%% =========================================================================
% AUXILIARY MODEL SOLUTION SUBFUNCTION (Using Alternative Conditions Grid)
% =========================================================================
function Solutions = ModelSolution(x02, lambda_s, eta, TFtSpan)
    options = odeset('Events', @(t,x) myEvents(t,x,lambda_s,eta), 'RelTol', 1e-1, 'AbsTol', 1e-1);
    Solutions = ode23s(@(t,x) ModelEquations(t,x,lambda_s,eta), TFtSpan, x02, options);
end

%% =========================================================================
% ODE EVENTS MONITOR SUBFUNCTION
% =========================================================================
function [value, isterminal, direction] = myEvents(t, x, lambda_s, eta)
    DerivSol = ModelEquations(t, x, lambda_s, eta);
    value = [DerivSol(1) + 1e-1, DerivSol(2) + 1e-1];
    isterminal = [0, 0]; % 0: Do not stop simulation when event triggers
    direction = [-1, -1]; 
end

