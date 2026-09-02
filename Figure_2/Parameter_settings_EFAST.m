%% PARAMETER INITIALIZATION

% PARAMETER BASELINE VALUES
lambdaS = 0.075; % 
K = 10^9;
beta0 = 2.3e-9; % 
deltaC = 2.68e-9; % 
hC = 40; % 0.5;
gamma0 = 0.06;
b =30;
omega = 1.2;
phi = 0.7e2;
hM = 2e7;
dM = 0.25;
eta = 0.1;
hD =  2.019e7;
psi = 0.0005;
dC =  0.02;
dummy=1; %25




% Redefine parameters for circuit topologies with NAR and PAR, so that all
% three synthesis same steady state conc of FadR in absence of inducer:
    %params = p2;
    
    params = [lambdaS,K,beta0,deltaC,hC,gamma0,b,omega,phi,hM,...
            dM,eta,hD,psi,dC,dummy];
    p_var_pos = [1,2,3,4,5,6,17,18];

    % Parameter labels:
    % Parameter Labels 
    efast_var={'\lambda_s','K','\beta','\delta_c','h_c','\gamma','b',...
    '\omega','\phi','h_m','d_m','\eta','h_d','\psi','d_c','dummy'};%,

    % Define min and max vectors of parameters:
    pmin = params * 0.1; % 0.1
    pmax = params * 5; % 10
    pmin(end-1) = 0.1; % for param sT
    pmax(end-1) = 0.9; % for param sT



%% TIME SPAN OF THE SIMULATION
t_end = 400; % length of the simulations
tspan=(0:1:t_end);   % time points where the output is calculated
time_points = 200; 
%time_points=[9 13 70 200]; % time points of interest for the model analysis

% INITIAL CONDITION FOR THE ODE MODEL 

y0=[1e6,0,1e6,0,1e6];

% Variables Labels
y_var_label={'T(t)','T_i(t)','V(t)','M(t)','C_d(t)','Tumor population'};
