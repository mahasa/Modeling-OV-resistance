%% This ODE represents the CTC model 
function dy=ODEmodel(t,y,X,run_num)

    %% PARAMETERS %%
    Parameter_settings_EFAST;

    lambdaS=X(run_num,1);
    K=X(run_num,2);
    beta0=X(run_num,3);
    deltaC=X(run_num,4);
    hC=X(run_num,5);
    gamma0=X(run_num,6);
    b=X(run_num,7);
    omega=X(run_num,8);
    phi=X(run_num,9);
    hM=X(run_num,10);
    dM=X(run_num,11);
    eta=X(run_num,12);
    hD=X(run_num,13);
    psi=X(run_num,14);
    dC=X(run_num,15);
    dummy=X(run_num,16);
    
    
    dy=zeros(5,1);
         
    % Difine ODEs here
    dy(1) = lambdaS*y(1)*(1-y(1)./K) - beta0*y(1)*y(3) - (deltaC*y(1)*y(5))./(hC + y(5));
    dy(2) = beta0*y(1)*y(3) - (deltaC*y(2)*y(5))./(hC + y(5)) - gamma0*y(2);
    dy(3) = b*gamma0*y(2) - omega*y(3);
    dy(4) = (phi*y(2))./(hM + y(2)) - dM*y(4);
    dy(5) = ((eta*(y(1)+y(2)))./(hD + y(1)+y(2)))*(1-psi*y(4))*y(5) -  dC*y(5);

         

