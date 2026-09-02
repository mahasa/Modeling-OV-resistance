function Fig8_OV_Resistance_Model_Varying_Beta
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Written by: Khaphetsi Joseph Mahasa
    % Date: 03 Jul 2026
    % Function: Code to generate Figure 8 in the article.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Important notes %%%%%%%%%%%
    % Vary MOI (according to Liu et al (2022)) initial number of OVs in x0 
    % Figure 8 in the manuscript
    % Scenario 2: High infection rate (beta0 = 0.5e-6; % p3), strong immune response (eta = 0.5; % p12)
    %%
    
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

    mFileWorkingDirectory=pwd;
    mFileParentDirectory=mFileWorkingDirectory;
    
    format shortg
    clk = datestr(datetime('now'),'dd-mmm-yy-HH');
    PlotsFolder=['OV_Resistance_PlotsDirectory-',clk];
    
    if ~exist(fullfile(mFileParentDirectory, PlotsFolder), 'dir')
      mkdir(mFileParentDirectory, PlotsFolder);
    end
    
    ResultsFolder=['OV_Resistance_ResultsDirectory-',clk];
    if ~exist(fullfile(mFileParentDirectory, ResultsFolder), 'dir')
      mkdir(mFileParentDirectory, ResultsFolder);
    end

    PlotsFullPath=fullfile(mFileParentDirectory,PlotsFolder);
    ResultsFullPath=fullfile(mFileParentDirectory,ResultsFolder);
    diary(fullfile(ResultsFullPath,'OV_Resistance_Results.txt'));
    
    %% define parameters here
    lambda_s = 0.075; 
    K = 5.3196e8;  
    beta0 = 0.5e-6; 
    delta_c = 2.68e-9; 
    h_c = 40; 
    gamma0 = 0.05; 
    b = 30; 
    omega = 1.2; 
    phi = 10e8; 
    h_m = 2e7; 
    d_m = 0.25; 
    eta = 0.5; 
    h_d = 2.019e7; 
    psi = 0; 
    d_c = 0.02; 
   
    
    %% Primary model solution
    tStart=0;
    tEnd=100;
    MOI = [0.001, 0.01, 0.1, 1, 10];
    PFU = [20, 200, 2000, 20000, 200000];
    tSpan = [0 300]; 
    TFtSpan=tStart:180;

    %% -----------------------------------------------------------------------
    %                           Generate the figures
    % ------------------------------------------------------------------------  
    cd(PlotsFullPath);
    Styles={'-','--','-.',':','-','--','-.',':','-','--','-.',':'}; 
    
    options = odeset('Refine',10, 'RelTol',1e-4);

    % START FOR LOOP FOR EACH PFU CONDITION
    for i = 1:length(PFU)
        % Dynamic initial condition vector setup
        x0 = [2e6, 0, PFU(i), 0, 1e6];  
        
        [t, x] = ode23s(@ModelEquations, tSpan, x0, options);   
        
        FUnT = x(:,1);
        FInfT = x(:,2);
        FV = x(:,3);
        FM = x(:,4);
        FCD = x(:,5);

        % Create a brand new distinct figure window for this loop iteration
        figure(i)
        hold on; box on;
        plot(t, FUnT + FInfT, 'LineStyle',Styles{1},'LineWidth',2.5)  
        plot(t, FCD, 'LineStyle',Styles{4},'LineWidth',2.5)  
        plot([0,300],[1e6,1e6],'k--') 

        ylim([0 9e6])  
        xlim([tSpan(1) tSpan(2)])
        legend('T(t)+T_{i}(t)', 'C_{d}(t)','Threshold','location','SouthEast','Interpreter', 'latex')
        xlabel('Time (days)', 'Interpreter', 'latex')
        ylabel('Population size (\# cells)', 'Interpreter', 'latex')
        
        % Dynamic titles using index trackers
        title(sprintf('MOI: %g', MOI(i)), 'Interpreter', 'latex')
             
        %%% inset figure 1 %%%
        axes('position',[.315 .55 .25 .28], 'NextPlot', 'add')
        box on 
        indexOfInterest = (t <= 200) & (t >= 0); 
        plot(t(indexOfInterest), FV(indexOfInterest),'LineStyle',Styles{2},'Color','r', 'LineWidth', 2)
        legend('V(t)','Interpreter', 'latex')
        ylabel('PFU/cell', 'Interpreter', 'latex')
        xlabel('Time (days)', 'Interpreter', 'latex')
        set(gca,'YScale','log')
        axis tight

        %%% inset figure 2 %%%
        axes('position',[.65 .55 .25 .28], 'NextPlot', 'add')
        box on 
        plot(t(indexOfInterest), FM(indexOfInterest),'LineStyle',Styles{3},'Color','g', 'LineWidth', 2)
        legend('M(t)','Interpreter', 'latex')
        ylabel('\# cells', 'Interpreter', 'latex')
        xlabel('Time (days)', 'Interpreter', 'latex')
        set(gca,'YScale','log')
        axis tight
        hold off
        
        % Save each individual image cleanly matching the exact MOI run
        filename = sprintf('Figure_MOI_%g.png', MOI(i));
        % Replacing decimal dots with 'p' in filename if desired to prevent extension bugs
        filename = strrep(filename, '.', 'p'); 
        filename = strrep(filename, 'png', '.png'); % correct double dot cleanup
        saveas(gcf, filename);
    end 
    % END FOR LOOP
       
    diary off
    cd(mFileWorkingDirectory);
end 

%% ========================================================================
%  SUBFUNCTIONS SECTION
%% ========================================================================

function dx = ModelEquations(~, x) 
    global lambda_s K beta0 delta_c h_c gamma0 b omega phi h_m d_m eta h_d psi d_c 
    dx = zeros(5,1);  
    
    dx(1) = lambda_s*x(1)*(1-x(1)/K) - beta0*x(1)*x(3) - (delta_c*x(1)*x(5))/(h_c + x(5));
    dx(2) = beta0*x(1)*x(3) - (delta_c*x(2)*x(5))/(h_c + x(5)) - gamma0*x(2);
    dx(3) = b*gamma0*x(2) - omega*x(3);
    dx(4) = (phi*x(2))/(h_m + x(2)) - d_m*x(4); 
    dx(5) = ((eta*(x(1)+x(2)))/(h_d + x(1)+x(2)))*(1-psi*x(4))*x(5) - d_c*x(5);
end


function [value, isterminal, direction] = myEvents(t, x)
    DerivSol = ModelEquations(t, x);
    value = [DerivSol(1) + 1e-01, DerivSol(2) + 1e-01];
    isterminal = [0,0];
    direction = [-1, -1];
end
