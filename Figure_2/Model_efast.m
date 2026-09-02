% First order and total effect indices for a given
% model computed with Extended Fourier Amplitude
% Sensitivity Test (EFAST).
% 
clc;clear;close all
%% INPUT
NR = 1; %: no. of search curves - RESAMPLING
k = 15 + 1; % # of input factors (parameters varied) + dummy parameter
NS = 80; %65; % # of samples per search curve
wantedN=NS*k*NR; % wanted no. of sample points

% OUTPUT
% SI[] : first order sensitivity indices
% STI[] : total effect sensitivity indices
% Other used variables/constants:
% OM[] : vector of k frequencies
% OMi : frequency for the group of interest
% OMCI[] : set of freq. used for the compl. group
% X[] : parameter combination rank matrix
% AC[],BC[]: fourier coefficients
% FI[] : random phase shift
% V : total output variance (for each curve)
% VI : partial var. of par. i (for each curve)
% VCI : part. var. of the compl. set of par...
% AV : total variance in the time domain
% AVI : partial variance of par. i
% AVCI : part. var. of the compl. set of par.
% Y[] : model output

MI = 4; %: maximum number of fourier coefficients
% that may be retained in calculating the partial
% variances without interferences between the
% assigned frequencies

%% PARAMETERS AND ODE SETTINGS (they are included in the following file)
Parameter_settings_EFAST;


% Computation of the frequency for the group
% of interest OMi and the # of sample points NS (here N=NS)
OMi = floor(((wantedN/NR)-1)/(2*MI)/k);
NS = 2*MI*OMi+1;
if(NS*NR < 65)
    fprintf(['Error: sample size must be >= ' ...
    '65 per factor.\n']);
    return;
end


%% Pre-allocation of the output matrix Y
%% Y will save only the points of interest specified in
%% the vector time_points
Y(NS,length(time_points),length(y0),length(pmin),NR)=0;  % pre-allocation

% Loop over k parameters (input factors)
tic
for i=1:k % i=# of replications (or blocks)
    % Algorithm for selecting the set of frequencies.
    % OMci(i), i=1:k-1, contains the set of frequencies
    % to be used by the complementary group.
    OMci = SETFREQ(k,OMi/2/MI,i);   
    % Loop over the NR search curves.
    for L=1:NR
        % Setting the vector of frequencies OM
        % for the k parameters
        cj = 1;
        for j=1:k
            if(j==i)
                % For the parameter (factor) of interest
                OM(i) = OMi;
            else
                % For the complementary group.
                OM(j) = OMci(cj);
                cj = cj+1;
            end
        end
        % Setting the relation between the scalar
        % variable S and the coordinates
        % {X(1),X(2),...X(k)} of each sample point.
        FI = rand(1,k)*2*pi; % random phase shift
        S_VEC = pi*(2*(1:NS)-NS-1)/NS;
        OM_VEC = OM(1:k);
        FI_MAT = FI(ones(NS,1),1:k)';
        ANGLE = OM_VEC'*S_VEC+FI_MAT;
        
        X(:,:,i,L) = 0.5+asin(sin(ANGLE'))/pi; % between 0 and 1
        
        % Transform distributions from standard
        % uniform to general.
        X(:,:,i,L) = parameterdist(X(:,:,i,L),pmax,pmin,0,1,NS,'unif'); %%this is what assigns 'our' values rather than 0:1 dist
        % Do the NS model evaluations.
        for run_num=1:NS
            [i run_num L] % keeps track of [parameter run NR]
            % ODE system file
            f=@ODE_efast;
            % ODE solver call 
            options = odeset('Refine',10, 'RelTol',1e-4);
            [t,y]=ode23s(@(t,y)f(t,y,X(:,:,i,L),run_num),tspan,y0,options);
            %[t,y]=ode15s(@(t,y)f(t,y,X(:,:,i,L),run_num),tspan,y0,[]);
            % It saves only the output at the time points of interest
            Y(run_num,:,:,i,L)=y(time_points+1,:);
        end %run_num=1:NS
    end % L=1:NR
end % i=1:k
%save Model_efast.mat;
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % CALCULATE Si AND STi for each resample (1,2,...,NR) [ranges]
% [Si,Sti,rangeSi,rangeSti] = efast_sd(Y,OMi,MI,time_points,1:length(y0));
% 
% % Calculate Coeff. of Var. for Si and STi for Viral load (variable 4). See
% % online Supplement A.5 for details.
% [CVsi CVsti] = CVmethod(Si,rangeSi,Sti,rangeSti,1);
% 
% % T-test on Si and STi for Viral load (variable 4)
% s_KM = efast_ttest(Si,rangeSi,Sti,rangeSti,1:length(time_points),efast_var,2,y_var_label,0.05);


% % PLOTTING RESULTS
% figure(); clf
% m_BarPlot_Sensitivities
% %title({['eFAST (GSA) to find sensitivity of: ',char(y_var_label)]; ['Circuit: ', char(arch)]})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE Si AND STi for each resample (1,2,...,NR) [ranges]
% Primary tumor
[Si1,Sti1,rangeSi1,rangeSti1] = efast_sd(Y,OMi,MI,time_points,1:4)
% Calculate Coeff. of Var. for Si and STi for Uninfected tumor (variable 1).
[CVsi1 CVsti1]=CVmethod(Si1, rangeSi1,Sti1,rangeSti1 ,1)


% Circulating tumor
[Si2,Sti2,rangeSi2,rangeSti2] = efast_sd(Y,OMi,MI,time_points,1:4)
% Calculate Coeff. of Var. for Si and STi for Infected tumor (variable 2). 
[CVsi2 CVsti2]=CVmethod(Si2, rangeSi2,Sti2,rangeSti2,2)

% Total tumor cell population
TotalTumor_Si = Si1 + Si2;
TotalTumor_Sti = Sti1 + Sti2;
TotalTumor_rangeSi = rangeSi1 + rangeSi2;
TotalTumor_rangeSti = rangeSti1 + rangeSti2;
%[TotalTumor_CVsi TotalTumor_CVSti] = [CVsi1 CVsi2 CVsti1 CVsti2];
AllTumor_Si = TotalTumor_Si(:,:,1) + TotalTumor_Si(:,:,2);
AllTumor_Sti = TotalTumor_Sti(:,:,1) + TotalTumor_Sti(:,:,2);
AllTumor_rangeSi = TotalTumor_rangeSi(:,:,1) + TotalTumor_rangeSi(:,:,2);
AllTumor_rangeSti = TotalTumor_rangeSti(:,:,1) + TotalTumor_rangeSti(:,:,2);


% T-test on Si and STi for CTC population (variable 2)
%s_Tumor = efast_ttest(TotalTumor_Si,TotalTumor_rangeSi,TotalTumor_Sti,TotalTumor_rangeSti,1:length(time_points),efast_var,VarOI,y_var_label,0.05)

s_Tumor = efast_ttest(AllTumor_Si,AllTumor_rangeSi,AllTumor_Sti,AllTumor_rangeSti,1:length(time_points),efast_var,1,y_var_label(6),0.01)

% Preparing results to plot average and 1 standard deviation from search
% repeats:
% ... first order sensitivities:
Si = squeeze(s_Tumor.rangeSi);
Si_avg = mean(Si,2);
Si_std = std(Si,0,2);
% ... total order sensitivities:
Sti = squeeze(s_Tumor.rangeSti);
Sti_avg = mean(Sti,2);
Sti_std = std(Sti,0,2);

Days = ["at day= 200"]; % Time points

% colors
colr = autumn(length(Si(1,:)));
colr = colr(end:-1:1,:);
%colr = winter(length(Si(1,:)));
%colr = colr(end:-1:1,:);
%Plot Si 
figure(1),
% % for i = 1:length(Si(1,:))
% %     hold on
% %     bar(AllTumor_Si(i,1),'FaceColor',colr(i,:),'EdgeColor','k','LineWidth',1)
% %     hold off
% % end
bar(AllTumor_Si(:,1),'FaceColor',[0.517647058823529,0.568627450980392,0.705882352941177],'EdgeColor','k','LineWidth',1)
hold on
% errorbar(1:length(Sti_avg),Sti_avg,Sti_std,'k.')
for i = 1:length(Si(1,:))
    plot(1:length(Si_avg),Si(:,i),'ko','MarkerSize',6,'LineWidth',1)
end
errorbar(1:length(Si_avg),Si_avg,Si_std,'k.')
hold off
set(gca,'FontSize',12,'FontName','Arial','FontWeight','Bold','LineWidth',1.0 )
xlim([0 17])
set(gca,'XTickLabel',efast_var,'XTick',[1:length(efast_var)])
xtickangle(45)
xlabel('Model parameters')
ylabel('eFAST sensitivity index')
% legend('Total order S_{Ti}','First order S_{i}','Location','SouthEast')
legend('First-order index (S_{i})','Location','NorthEast')
ylim([0 1])

% Indentifying parameters with sensitivities significantly different from
% dummy parameter:
p_sigDiff = find(s_Tumor.p_Si < 0.01);
hold on
text(p_sigDiff,Si_avg(p_sigDiff) + Si_std(p_sigDiff),'*','FontSize',20,'VerticalAlignment','bottom','HorizontalAlignment','center');
hold off


%Plot Sti              
figure(2), bar(AllTumor_Sti(:,1),'FaceColor',[0.8,0.8,0.8],'EdgeColor','k','LineWidth',1)
hold on
% errorbar(1:length(Sti_avg),Sti_avg,Sti_std,'k.')
for i = 1:length(Sti(1,:))
    plot(1:length(Sti_avg),Sti(:,i),'ko','MarkerSize',6,'LineWidth',1)
end
errorbar(1:length(Sti_avg),Sti_avg,Sti_std,'k.')
hold off
set(gca,'FontSize',12,'FontName','Arial','FontWeight','Bold','LineWidth',1.0 )
xlim([0 17])
set(gca,'XTickLabel',efast_var,'XTick',[1:length(efast_var)])
xtickangle(45)
colormap(jet(6));
xlabel('Model parameters')
ylabel('eFAST sensitivity index')
% legend('Total order S_{Ti}','First order S_{i}','Location','SouthEast')
legend('Total order index','Location','NorthEast')
ylim([0 1])

% Indentifying parameters with sensitivities significantly different from
% dummy parameter:
p_sigDiff = find(s_Tumor.p_Sti < 0.01);
hold on
text(p_sigDiff,Sti_avg(p_sigDiff) + Sti_std(p_sigDiff),'*','FontSize',20,'VerticalAlignment','bottom','HorizontalAlignment','center');
hold off

% colormap(jet(6));
% %im_hatch = applyhatch_pluscolor(gcf,'|x.+\/',0,[1 1 0 1 0 0],[],200,3,2);
% applyhatch_pluscolor(gcf,'|x.\+c',0,[1 1 0 1 0 0],jet(6),200,3,2);


% %Plot Si
% figure, bar(AllTumor_Si,'stacked','linewidth',1.5)
% %bar(AllTumor_Si,'FaceColor',['b' 'flat'],'EdgeColor','k','LineWidth',1)
% %figure, bar(Sti(:,:,4))
% %title('PRCCs for parameters at day 2,5 and 7');
% set(gca,'FontSize',12,'FontName','Arial','FontWeight','Bold','LineWidth',1.0 )
% xlim([0 27])
% set(gca,'XTickLabel',efast_var,'YScale','log')
% xlabel('Model parameters')
% ylabel('Measured sensitivity')
% % legend('Total order S_{Ti}','First order S_{i}','Location','SouthEast')
% legend('First order S_{i}','Location','SouthEast')
% ylim([0 1])
% set(gca, 'xtick',[1:length(efast_var)], 'XTickLabel', efast_var,'XTickLabelRotation',45); 
% colormap gray
% set(gca,'Linewidth',1.5);

% % Plot the Sti
% figure, bar(AllTumor_Sti,'stacked','linewidth',1.5)
% %title('PRCCs for parameters at day 2,5 and 7');
% xlabel('Parameter'),ylabel('eFAST total-order index (S_{ti})');
% legend(Days(1:end),'location','NorthWest','FontSize',10);
% set(gca,'FontSize',12,'FontName','Arial','FontWeight','Bold','LineWidth',1.0 )
% xlim([0 25])
% set(gca, 'xtick',[1:length(efast_var)], 'XTickLabel', efast_var,'XTickLabelRotation',45); 
% colormap gray
% set(gca,'Linewidth',1.5);

% % Plot the Sti
% figure, bar(AllTumor_Sti,'stacked','linewidth',1.5)
% %title('PRCCs for parameters at day 2,5 and 7');
% xlabel('Parameter'),ylabel('eFAST total-order index (S_{ti})');
% legend(Days(1:end),'location','NorthWest','FontSize',10);
% set(gca,'FontSize',12,'FontName','Arial','FontWeight','Bold','LineWidth',1.0 )
% xlim([0 25])
% %set(gca,'XTickLabel',efast_var,'XTick',[1:length(efast_var)])
% set(gca,'XTickLabel',efast_var,'XTick',[1:length(efast_var)])
% xtickangle(45)
% % colormap(jet(6));
% % %im_hatch = applyhatch_pluscolor(gcf,'|x.+\/',0,[1 1 0 1 0 0],[],200,3,2);
% % applyhatch_pluscolor(gcf,'|x.\+c',0,[1 1 0 1 0 0],jet(6),200,3,2);
% 
