function Fig3_SynergyHeatmaps()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Written by: Khaphetsi Joseph Mahasa
	% Date: 03 Jul 2026
	% Function: Code to generate Figure 3 in the article and Supplementary Fig. S1.
	%
	% Description: This script simulates global sensitivity landscapes
	% of the synergy quantity evaluated at t = 300 days.
	%
	% Visual Outputs (Saved to 'Fig3-Plots'):
	%   1. BurstInfection_Tumor_Efficacy_2D.png : Flat pcolor heatmap
	%   2. BurstInfection_Tumor_Efficacy_3D.png : 3D Surface mesh map
	%   3. Tumor_Synergy_TimeCourse.png         : Trajectories & Synergy Index
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % --- Fixed Model Parameters ---
    K = 10^9;        
    beta0 = 2.3e-10; 
    h_c = 40;        
    b = 3000;        
    omega = 1.2;     
    phi = 0.7e2;     
    h_m = 2e7;       
    d_m = 0.25;      
    h_d = 2.019e7;   
    psi = 5e-4;       % Active feedback metric
    d_c = 0.02;      
    
    lambda_s = 0.045;   
    eta = 0.1;          

    % FIXED: Initial conditions matched perfectly to 2e7
    y0 = [2e7; 0; 2000; 0; 1e6]; 
    Z0_minus = y0(1) + y0(2); 

    opt = odeset('AbsTol',1e-9,'RelTol',1e-6);
    tspan1 = linspace(0,300,1000); 
    
    %% Parameters that vary for heatmaps:    
    delta_c_sweep = logspace(-8, -3, 11);        
    gamma0_sweep  = linspace(0.001605, 2.667, 11);   

    Z1 = zeros(length(gamma0_sweep), length(delta_c_sweep));

    % 1. Run Surface Parameter Sweep Mesh 
    for ii = 1:length(delta_c_sweep)
        for kk = 1:length(gamma0_sweep)
            % ODEs
            odefcn = @(t,x)[lambda_s*x(1)*(1-x(1)./K) - beta0*x(1)*x(3) - (delta_c_sweep(ii)*x(1)*x(5))./(h_c + x(5));
               beta0*x(1)*x(3) - (delta_c_sweep(ii)*x(2)*x(5))./(h_c + x(5)) - gamma0_sweep(kk)*x(2);
               b*gamma0_sweep(kk)*x(2) - omega*x(3);
               (phi*x(2))/(h_m + x(2)) - d_m*x(4);
               ((eta*(x(1)+x(2)))/(h_d + x(1)+x(2)))*(1 - psi*x(4))*x(5) -  d_c*x(5)];
            
            [~, x_temp] = ode23s(odefcn, tspan1, y0, opt);
            T_final = x_temp(end, 1);
            Z1(kk, ii) = T_final;
        end
    end

    %% 2. Generate Heatmaps
    folderName = 'Fig3-Plots';
    if ~exist(folderName, 'dir'), mkdir(folderName); end

    x_contour = gamma0_sweep';
    y_contour = delta_c_sweep';
    [Y,X] = meshgrid(y_contour,x_contour);

    % 2D plot
    fig=figure(1); clf; ax = axes;
    set(ax, 'FontSize', 18, 'FontName', 'Arial');
    pcolor(Y,X,Z1); hold on; shading interp;
    cb1 = colorbar; cb1.LineWidth = 1.5; ax1=gca; set(gca,'XScale','log'); 
    xlabel('\delta_{c} (days^{-1})'), ylabel('\gamma (days^{-1})'); zlabel('Tumor Cell Count'); 
    set(gca,'linewidth', 2,'fontsize',16,'TickDir','out');
    h = axes(fig,'visible','off'); bottom = min(min(Z1)); top = max(max(Z1));
    caxis(ax1,[bottom, top]); cb1.Position = [0 0 0 0];                   
    c = colorbar(h,'Position',[0.92 0.23 0.022 0.7]); c.LineWidth = 1.5; c.FontSize = 18;
    caxis(get(ax1, 'CLim'));                    
    print(fig, fullfile(folderName, 'BurstInfection_Tumor_Efficacy_2D.png'), '-dpng', '-r300');

    % 3D plot
    fig2=figure(2); clf; ax2 = axes;
    set(ax2, 'FontSize', 18, 'FontName', 'Arial');
    surf(Y,X,Z1); hold on; shading interp;
    cb2 = colorbar; cb2.LineWidth = 1.5; ax3=gca; set(gca,'XScale','log'); 
    xlabel('\delta_{c} (days^{-1})'), ylabel('\gamma (days^{-1})'); zlabel('Tumor Cell Count'); 
    set(gca,'linewidth', 2,'fontsize',16,'TickDir','out'); view(3);                                    
    h2 = axes(fig2,'visible','off'); caxis(ax3,[bottom, top]); cb2.Position = [0 0 0 0];                   
    c3 = colorbar(h2,'Position',[0.92 0.23 0.022 0.7]); c3.LineWidth = 1.5; c3.FontSize = 18;
    caxis(get(ax3, 'CLim'));                    
    print(fig2, fullfile(folderName, 'BurstInfection_Tumor_Efficacy_3D.png'), '-dpng', '-r300');

    %% 3. FIXED TIME-COURSE TRAJECTORIES (Identical to Benchmark)
    % Forces the benchmark baseline single values instead of matrix index averages
    plot_delta = 2.68e-9; 
    plot_gamma = 0.06;   
    
    odefcn_single = @(t,x)[lambda_s*x(1)*(1-x(1)./K) - beta0*x(1)*x(3) - (plot_delta*x(1)*x(5))./(h_c + x(5));
                   beta0*x(1)*x(3) - (plot_delta*x(2)*x(5))./(h_c + x(5)) - plot_gamma*x(2);
                   b*plot_gamma*x(2) - omega*x(3);
                   (phi*x(2))/(h_m + x(2)) - d_m*x(4);
                   ((eta*(x(1)+x(2)))/(h_d + x(1)+x(2)))*(1 - psi*x(4))*x(5) -  d_c*x(5)];
               
    [t, x_trajectory] = ode23s(odefcn_single, tspan1, y0, opt);
    
    T  = x_trajectory(:, 1);
    Ti = x_trajectory(:, 2);
    
    Z_minus = (K * Z0_minus) ./ (Z0_minus + (K - Z0_minus) * exp(-lambda_s * t));
    S_t = (T + Ti) ./ Z_minus;
    
    fig3 = figure('Color', 'w');
    
    % Subplot 1: Cell Populations
    subplot(2, 1, 1);
    plot(t, T, 'b-', 'LineWidth', 2); hold on;
    plot(t, Ti, 'r--', 'LineWidth', 2);
    plot(t, Z_minus, 'k:', 'LineWidth', 2);
    grid on;
    xlabel('Time (t)');
    ylabel('Cell Count');
    legend('Uninfected (T)', 'Infected (T_i)', 'No Treatment Tumor (Z^-)', 'Interpreter', 'latex');
    title('Tumor Dynamics Evolution Under Treatment');

    % Subplot 2: Synergy Index S(t)
    subplot(2, 1, 2);
    plot(t, S_t, 'm-', 'LineWidth', 2);
    grid on;
    xlabel('Time (t)');
    ylabel('Synergy Index $(S_{\delta_c}^{\gamma}(t))$', 'Interpreter', 'latex');
    title('Treatment-Efficacy Ratio Over Time');
    
    print(fig3, fullfile(folderName, 'Tumor_Synergy_TimeCourse.png'), '-dpng', '-r300');
end
