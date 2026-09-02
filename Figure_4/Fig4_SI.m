function Fig4_SI
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Written by: Khaphetsi Joseph Mahasa
	% Date: 03 Jul 2026
	% Function: Code to generate Figure 4 in the article.
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
    
    % This is to creat subfolders for the plots and results. The folders' 
    % names include the current date and time to avoid overwritting them
    format shortg
    clk = datestr(datetime('now'),'dd-mmm-yy-HH');
    PlotsFolder=['\OV_Resistance_PlotsDirectory-',clk];
    
    % This is to avoid "creating" a folder that was exists already
    if ~exist('OV_Resistance_PlotsDirectory', 'dir')
      mkdir(mFileParentDirectory, PlotsFolder);
    end
    
    % Same here: The results folder name containes the current date. 
    ResultsFolder=['\OV_Resistance_ResultsDirectory-',clk];
    if ~exist('OV_Resistance_ResultsDirectory', 'dir')
      mkdir(mFileParentDirectory, ResultsFolder);
    end
    
    % The code's results go here. Anything that we would like to report
    PlotsFullPath=[mFileParentDirectory,PlotsFolder];
    ResultsFullPath=[mFileParentDirectory,ResultsFolder];
    diary([ResultsFullPath,'\OV_Resistance_Results.txt']);
    
    %% define parameters here
    lambda_s = 0.045; % p1
    K = 10^9;  % p2 
    beta0 = 2.3e-10; % 1e-3; p3
    delta_c = 2.68e-9; % p4
    h_c = 40; % p5
    gamma0 = 0.06; % p6
    b = 3000; % p7
    omega = 1.2; % p8
    phi =0.7e2; % p9
    h_m = 2e7; % p10
    d_m = 0.25; % p11
    eta = 0.1; % p12
    h_d = 2.019e7; % p13
    psi = 0.0;  % p14
    d_c = 0.02; % p15 mu2=9.12;

    
    
    %% Primary model solution
    tStart=0;
    tEnd=100;
    x0 = [ 2e7 0 2000 0 1e6];  % Initial conditions 
    x02 = [ 2e7 0 0 0 1e6];  % Initial conditions 
    tSpan = [0 300]; %time span of the simulation
    TFtSpan=tStart:180;
    
   
    function dx = ModelEquations(t,x,phi,psi) 

    %% Model equations 
    dx = [0; 0; 0; 0; 0];  %% preallocation of vector that keeps track of variables 
    %% systems of differential equations
    dx(1) = lambda_s*x(1)*(1-x(1)./K) - beta0*x(1)*x(3) - (delta_c*x(1)*x(5))./(h_c + x(5));
    dx(2) = beta0*x(1)*x(3) - (delta_c*x(2)*x(5))./(h_c + x(5)) - gamma0*x(2);
    dx(3) = b*gamma0*x(2) - omega*x(3);
    dx(4) = (phi*x(2))./(h_m + x(2)) - d_m*x(4);
    dx(5) = ((eta*(x(1)+x(2)))./(h_d + x(1)+x(2)))*(1-psi*x(4))*x(5) -  d_c*x(5);

    end


%% Model Solution
    function Solutions = ModelSolution(tStart,x02,phi,psi)
        options = odeset('Events',@(t,x) myEvents(t,x,phi,psi),'RelTol',1e-01,'AbsTol',1e-01);
        %options0 = odeset('RelTol',1e-03,'AbsTol',1e-03);
        Solutions = ode23s(@(t,x)ModelEquations(t,x,phi,psi),TFtSpan,x02,options);
        
    end

%% Events
    function [value,isterminal,direction] = myEvents(t,x,phi,psi)
        DerivSol=ModelEquations(t,x,phi,psi);
        value = [DerivSol(1) + 1e-01,DerivSol(2) + 1e-01];
        isterminal = [0,0];
        direction = [-1,-1];
    end

%% Evaluating the Model Solutions
    phiRate = [0 50 500 1e4 2e4]; % 
    psiRate = [0 5e-6 3e-5 1e-4 0.08];  % Supression rate of CD8+ T cells  

    %% -----------------------------------------------------------------------
    %                           Generate the figures
    % ------------------------------------------------------------------------  
   % Changing directory for saving the plots
    cd(PlotsFullPath);
    
    
    
    % Line Styles
    Styles={'-','--','-.',':','-','--','-.',':','-','--','-.',':'};
    
  for k=1:1
      % Figures' legends
      if k == 1
          lgd='Without virotherapy';
      end
  end
    
    
  options = odeset('Refine',10, 'RelTol',1e-4);
  [t,x] = ode23s(@(t,x)ModelEquations(t,x,phi,psi), tSpan, x0, options);   % Solution
 
  
  
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Plot the dynamics
%==========================================================================

figure; % Total Tumor and TCP 
    j=0;  % accumulator;
    subplot(1,2,1) % Total tumor (Uninfected and infected tumor cells)
    for k = 1: length(psiRate)
        psi = psiRate(k);
        phi = phiRate(k);
        j=j+1;
        SI = (psiRate(k)*phi*lambda_s)/(d_m*d_c);
        ttle=['I = ', num2str(SI,'%10.0e\n')];
        options = odeset('NonNegative',[1,2,3,4,5]);
        [t,x] = ode23s(@(t,x)ModelEquations(t,x,phi,psi), tSpan, x02, options);   % Solution: Without virotherapy
        FUnT = x(:,1); FInfT = x(:,2); %FV = x(:,3); %FM = x(:,4); %FCD = x(:,5);            
        [t1, x1] = ode23s(@(t,x)ModelEquations(t,x,phi,psi),tSpan,x0,options);
        T = x1(:,1); Ti = x1(:,2); V = x1(:,3); M = x1(:,4); Cd = x1(:,5); % 
        hold all
        h = plot(t,FUnT + FInfT,'LineStyle',Styles{1},'LineWidth',2.5); hold all grid
        h3(k) = plot(t1,T+Ti,'LineStyle',Styles{k},'LineWidth',2.5); hold all grid
        if k==1; ylabel('Total tumor size (\# cells)','Interpreter', 'latex'); end
        xlabel('Time (days)','Interpreter', 'latex')
        if j==1, hold on, end  % only set hold after first plot...
        legendInfo{j}=['\psi = ', num2str(psiRate(k)), ',  \phi = ', num2str(phiRate(k)), ',  I = ', num2str(SI)];
    end
    legend(h3, legendInfo,'location','SouthEast','FontSize',14)
    hold off;
    legend boxoff
    subplot(1,2,2) % TCP
    for k = 1: length(psiRate)
        psi = psiRate(k);
        phi = phiRate(k);
        j=j+1;
        options = odeset('NonNegative',[1,2,3,4,5]);          
        [t1, x1] = ode23s(@(t,x)ModelEquations(t,x,phi,psi),tSpan,x0,options);
        T = x1(:,1); Ti = x1(:,2); V = x1(:,3); M = x1(:,4); Cd = x1(:,5); % 
        % Calculate TCP for Tumor 
        TCP = exp(-T);
        hold all
        h0(k) = plot(t1,TCP,'LineStyle',Styles{k},'LineWidth',2.5); hold all grid
        if k==2; ylabel('Tumor control probability','Interpreter', 'latex'); end
        xlabel('Time (days)','Interpreter', 'latex')
        if j==1, hold on, end  % only set hold after first plot...
        legendInfo{j}=['\psi = ', num2str(psiRate(k)), ',  \phi = ', num2str(phiRate(k)), ',  I = ', num2str(SI)];
    end
    legend(h0, legendInfo,'location','NorthEast','FontSize',14)
    hold off;
    ylim([0 1]);
    legend boxoff
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; % TAMCs amd CD8+ T cells  
    j=0;  % accumulator;
    subplot(1,2,1)   % TAMCs
    for k = 1: length(psiRate)
        psi = psiRate(k);
        phi = phiRate(k);
        j=j+1;
        SI = (psiRate(k)*phi*lambda_s)/(d_m*d_c);
        ttle=['I = ', num2str(SI,'%10.0e\n')];
        options = odeset('NonNegative',[1,2,3,4,5]);           
        [t1, x1] = ode23s(@(t,x)ModelEquations(t,x,phi,psi),tSpan,x0,options);
        T = x1(:,1); Ti = x1(:,2); V = x1(:,3); M = x1(:,4); Cd = x1(:,5); % 
        hold all
        h5(k) = plot(t1,M,'LineStyle',Styles{k},'LineWidth',2.5); hold all grid
        if k==1; ylabel('TAMC size (\# cells)','Interpreter', 'latex'); end
        xlabel('Time (days)','Interpreter', 'latex')
        if j==1, hold on, end  % only set hold after first plot...

        legendInfo{j}=['\psi = ', num2str(psiRate(k)), ',  \phi = ', num2str(phiRate(k)), ',  I = ', num2str(SI)];
    end
    legend(h5,legendInfo,'location','NorthEast','FontSize',14)
    hold off;
    legend boxoff
    subplot(1,2,2) % CD8+ T cells
    for k = 1: length(psiRate)
        psi = psiRate(k);
        phi = phiRate(k);
        j=j+1;
        SI = (psiRate(k)*phi*lambda_s)/(d_m*d_c);
        ttle=['I = ', num2str(SI,'%10.0e\n')];
        options = odeset('NonNegative',[1,2,3,4,5]);         
        [t1, x1] = ode23s(@(t,x)ModelEquations(t,x,phi,psi),tSpan,x0,options);
        T = x1(:,1); Ti = x1(:,2); V = x1(:,3); M = x1(:,4); Cd = x1(:,5); % 
        hold all
        h6(k) = plot(t1,Cd,'LineStyle',Styles{k},'LineWidth',2.5); hold all grid
        if k==1; ylabel('CD$8^+$ T cell size (\# cells)','Interpreter', 'latex'); end
        xlabel('Time (days)','Interpreter', 'latex')
        if j==1, hold on, end  % only set hold after first plot...
        legendInfo{j}=['\psi = ', num2str(psiRate(k)), ',  \phi = ', num2str(phiRate(k)), ',  I = ', num2str(SI)];
    end
    legend(h6,legendInfo,'location','NorthEast','FontSize',14)
    hold off;
    ylim([0,inf])
    legend boxoff    
    
diary off
toc
%% returning to the working directory
cd(mFileWorkingDirectory);
end
