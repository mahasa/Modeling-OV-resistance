function Fig5_SR_lambdaS_Beta
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Written by: Khaphetsi Joseph Mahasa
    % Date: 03 Jul 2026
    % Function: Code to generate Figure 5 in the article.
    %
    % Description: This script simulates and plots the dual-panel dynamics 
    %              of total tumor size alongside Tumor Control Probability (TCP). 
    %              It iterates across varying viral infection rates (\beta_0) 
    %              and tumor proliferation rates (\lambda_s) over a 300-day 
    %              timeline, using custom event tracking to evaluate stability.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tic 
    global lambda_s K beta0 delta_c h_c gamma0 b omega phi h_m d_m eta h_d psi d_c 
    global x0 tStart tEnd tSpan

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

    mFileWorkingDirectory=pwd;
    mFileParentDirectory=mFileWorkingDirectory;
    
    % Use filesep for cross-platform folder generation
    format shortg
    clk = datestr(datetime('now'),'dd-mmm-yy-HH');
    PlotsFolder=[filesep, 'OV_Resistance_PlotsDirectory-',clk];
    
    if ~exist('OV_Resistance_PlotsDirectory', 'dir')
      mkdir(mFileParentDirectory, PlotsFolder);
    end
    
    ResultsFolder=[filesep, 'OV_Resistance_ResultsDirectory-',clk];
    if ~exist('OV_Resistance_ResultsDirectory', 'dir')
      mkdir(mFileParentDirectory, ResultsFolder);
    end
    
    PlotsFullPath=[mFileParentDirectory,PlotsFolder];
    ResultsFullPath=[mFileParentDirectory,ResultsFolder];
    diary([ResultsFullPath, filesep, 'OV_Resistance_Results.txt']);
    
    %% define parameters here
    lambda_s = 0.075; % p1
    K = 10^9; % p2
    % beta0 is now swept in the loop below
    delta_c = 2.68e-9; % p4
    h_c = 40; % p5
    gamma0 = 0.05; % p6
    b = 3000; % p7
    omega = 1.2; % p8
    phi = 0.7e2; % p9
    h_m = 2e7; % p10
    d_m = 0.25; % p11
    eta = 0.1; % p12
    h_d = 2.019e7; % p13
    psi = 0.005; % p14
    d_c = 0.02; % p15
    
    %% Primary model solution
    tStart=0;
    tEnd=100;
    x0 = [2e7, 0, 2000, 0, 1e6];  % Initial conditions 
    tSpan = [0 300]; % time span of the simulation
    
    %% Evaluating the Model Solutions
    lambdaSRate = [0.065 0.070 0.075 0.085 0.095];  % Suppression rate arrays
    betaValues = [2e-10, 2.3e-12];                  % Infection rates to compare

    % Changing directory to save plots
    cd(PlotsFullPath);
    
    % Line Styles
    Styles={'-','--','-.',':','-','--','-.',':','-','--','-.',':'};
     
    %% Outer Loop: Iterate through each beta0 value
    for bIdx = 1:length(betaValues)
        beta0 = betaValues(bIdx);
        betaStr = num2str(beta0,'%10.1e');
        
        % =====================================================================
        % Plot Figure A: Total Tumor vs TCP
        % =====================================================================
        figure('Name', ['Tumor_Dynamics_Beta_' betaStr]); 
        
        % Subplot 1: Total tumor
        subplot(1,2,1) 
        legendInfo1 = cell(1, length(lambdaSRate));
        h3 = zeros(1, length(lambdaSRate));
        for k = 1: length(lambdaSRate)
            lambda_s = lambdaSRate(k);
            eta = 0.1;
            SR = eta*delta_c/(lambdaSRate(k)*d_c);
            ttle=['\beta = ', num2str(beta0,'%10.1e')];
            options = odeset('NonNegative',[1,2,3,4,5]);          
            [t1, x1] = ode23s(@(t,x)ModelEquations(t,x,lambda_s,eta),tSpan,x0,options);
            T = x1(:,1); Ti = x1(:,2); 
            
            hold on;
            h3(k) = plot(t1,T+Ti,'LineStyle',Styles{k},'LineWidth',2.5); 
            grid on;
            if k==1, ylabel('Total tumor size (\# cells)','Interpreter', 'latex'); end
            xlabel('Time (days)','Interpreter', 'latex')
            legendInfo1{k}=['\lambda_s = ', num2str(lambdaSRate(k)), ',  A = ', num2str(SR)];
        end
        title(ttle);
        legend(h3, legendInfo1,'location','SouthEast','FontSize',12)
        legend boxoff
        
        % Subplot 2: TCP
        subplot(1,2,2) 
        legendInfo2 = cell(1, length(lambdaSRate));
        h0 = zeros(1, length(lambdaSRate));
        for k = 1: length(lambdaSRate)
            lambda_s = lambdaSRate(k);
            eta = 0.1;
            SR = eta*delta_c/(lambdaSRate(k)*d_c);
            options = odeset('NonNegative',[1,2,3,4,5]);          
            [t1, x1] = ode23s(@(t,x)ModelEquations(t,x,lambda_s,eta),tSpan,x0,options);
            T = x1(:,1); 
            TCP = exp(-T);
            
            hold on;
            h0(k) = plot(t1,TCP,'LineStyle',Styles{k},'LineWidth',2.5); 
            grid on;
            if k==1, ylabel('Tumor control probability','Interpreter', 'latex'); end
            xlabel('Time (days)','Interpreter', 'latex')
            legendInfo2{k}=['\lambda_s = ', num2str(lambdaSRate(k)), ',  A = ', num2str(SR)];
        end
        legend(h0, legendInfo2,'location','SouthEast','FontSize',12)
        ylim([0 1]);
        legend boxoff
        
        % Save Figure A 
        fig_tumor = gcf; 
        fileNameA = ['Fig5_Tumor_TCP_beta_' betaStr '.png'];
        print(fig_tumor, fullfile(PlotsFullPath, fileNameA), '-dpng', '-r300');
        
        % The editable MATLAB file
        saveas(fig_tumor, fullfile(PlotsFullPath, ['Fig5_Tumor_TCP_beta_' betaStr '.fig']));
        
        
        % =====================================================================
        % Plot Figure B: TAMCs vs CD8+ T cells
        % =====================================================================
        figure('Name', ['Immune_Dynamics_Beta_' betaStr]); 
        
        % Subplot 1: TAMCs
        subplot(1,2,1)   
        legendInfo3 = cell(1, length(lambdaSRate));
        h5 = zeros(1, length(lambdaSRate));
        for k = 1: length(lambdaSRate)
            lambda_s = lambdaSRate(k);
            eta = 0.1;
            SR = eta*delta_c/(lambdaSRate(k)*d_c);
            ttle=['\beta = ', num2str(beta0,'%10.1e')];
            options = odeset('NonNegative',[1,2,3,4,5]);
            [t1, x1] = ode23s(@(t,x)ModelEquations(t,x,lambda_s,eta),tSpan,x0,options);
            M = x1(:,4); 
            
            hold on;
            h5(k) = plot(t1,M,'LineStyle',Styles{k},'LineWidth',2.5); 
            grid on;
            if k==1, ylabel('TAMC size (\# cells)','Interpreter', 'latex'); end
            xlabel('Time (days)','Interpreter', 'latex')
            legendInfo3{k}=['\lambda_s = ', num2str(lambdaSRate(k)), ',  A = ', num2str(SR)];
        end
        title(ttle);
        legend(h5,legendInfo3,'location','NorthEast','FontSize',12)
        legend boxoff
        
        % Subplot 2: CD8+ T cells
        subplot(1,2,2) 
        legendInfo4 = cell(1, length(lambdaSRate));
        h6 = zeros(1, length(lambdaSRate));
        for k = 1: length(lambdaSRate)
            lambda_s = lambdaSRate(k);
            eta = 0.1;
            SR = eta*delta_c/(lambdaSRate(k)*d_c);
            ttle=['\beta = ', num2str(beta0,'%10.1e')];
            options = odeset('NonNegative',[1,2,3,4,5]);         
            [t1, x1] = ode23s(@(t,x)ModelEquations(t,x,lambda_s,eta),tSpan,x0,options);
            Cd = x1(:,5); 
            
            hold on;
            h6(k) = plot(t1,Cd,'LineStyle',Styles{k},'LineWidth',2.5); 
            grid on;
            if k==1, ylabel('CD$8^+$ T cell size (\# cells)','Interpreter', 'latex'); end
            xlabel('Time (days)','Interpreter', 'latex')
            legendInfo4{k}=['\lambda_s = ', num2str(lambdaSRate(k)), ',  A = ', num2str(SR)];
        end
        title(ttle);
        legend(h6, legendInfo4,'location','NorthEast','FontSize',12)
        ylim([0, inf]);
        legend boxoff
        
        % Save Figure B 
        fig_immune = gcf;
        fileNameB = ['Fig5_TAMC_CD8_beta_' betaStr '.png'];
        print(fig_immune, fullfile(PlotsFullPath, fileNameB), '-dpng', '-r300');
        
        % The editable MATLAB file
        saveas(fig_immune, fullfile(PlotsFullPath, ['Fig5_TAMC_CD8_beta_' betaStr '.fig']));
    end

    %% Wrap up
    diary off
    toc
    cd(mFileWorkingDirectory);
end

%% ========================================================================
%  Local Functions (Keep Outside Main Function Scope)
%  ========================================================================

function dx = ModelEquations(~,x,lambda_s,eta) 
    global K beta0 delta_c h_c gamma0 b omega phi h_m d_m h_d psi d_c 
    dx = zeros(5,1);  
    dx(1) = lambda_s*x(1)*(1-x(1)./K) - beta0*x(1)*x(3) - (delta_c*x(1)*x(5))./(h_c + x(5));
    dx(2) = beta0*x(1)*x(3) - (delta_c*x(2)*x(5))./(h_c + x(5)) - gamma0*x(2);
    dx(3) = b*gamma0*x(2) - omega*x(3);
    dx(4) = (phi*x(2))./(h_m + x(2)) - d_m*x(4);
    dx(5) = ((eta*(x(1)+x(2)))./(h_d + x(1)+x(2)))*(1-psi*x(4))*x(5) -  d_c*x(5);
end
