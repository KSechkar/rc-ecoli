%% Initial preparations

% decipher=containers.Map('KeyType','int','ValueType', 'char');
% decipher(1)='regulated';
% decipher(2)='fixed';
% decipher(3)='optimal';

Delta=0.001; % difference in growth rate not considered important
Max_iter=100; % max. number of iterations
sim = advanced_simulator;
sim.tf=7; % time frame of single integration (h)

close all

nutrients=[0.1,0.5,2]; %[0.05,0.1,0.25,0.5,1,2,10];
% storing growth rates
l_map = zeros(3,size(nutrients,2));

% storing ribosome transcription rates
F_r_map = zeros(3,size(nutrients,2));

% storing ribosome concentrations
R_map = zeros(3,size(nutrients,2));

fixed_F_rs= 0.01:0.01:0.6; %0.05:0.05:1; % range of fixed ribosome transcription coefficients
ls=zeros(size(nutrients,2),size(fixed_F_rs,2));
Rs=zeros(size(nutrients,2),size(fixed_F_rs,2));

for j=1:size(nutrients,2)
    sim.parameters('s')=sim.parameters('K_nus')*nutrients(j);

    % Run with regulated ribosome transcription
    sim.parameters('is_fixed_F_r')=0; % F_r regulated now
    [l, F_r, R] = get_steady_for_sweep(sim,Delta,Max_iter);
    l_map(1,j)=l;
    F_r_map(1,j)=F_r;
    R_map(1,j)=R;
    sim.init_conditions=advanced_init_conds(sim.parameters); % reset intial conditions % REDUCED
    %sim.init_conditions=full_init_conds(sim.parameters); % reset intial conditions % FULL
    
    % Run with fixed ribosome transcription
    sim.parameters('is_fixed_F_r')=1; % F_r fixed now
    sim.parameters('fixed_F_r')=0.6;
    [l, F_r, R] = get_steady_for_sweep(sim,Delta,Max_iter);
    l_map(2,j)=l;
    F_r_map(2,j)=F_r;
    R_map(2,j)=R;
    sim.init_conditions=advanced_init_conds(sim.parameters); % reset intial conditions % REDUCED
    %sim.init_conditions=full_init_conds(sim.parameters); % reset intial conditions % FULL

    % Run with a range of fixed ribosome transcriptions
    sim.parameters('is_fixed_F_r')=1; % F_r fixed now
    for i=1:size(fixed_F_rs,2)
        disp(fixed_F_rs(i))
        sim.parameters('fixed_F_r')=fixed_F_rs(i);
        [l, F_r, R] = get_steady_for_sweep(sim,Delta,Max_iter);
        ls(j,i)=l;
        Rs(j,i)=R;
        sim.init_conditions=advanced_init_conds(sim.parameters); % reset intial conditions % REDUCED
        %sim.init_conditions=full_init_conds(sim.parameters); % reset intial conditions % FULL
    end
    % finding optimal ribosome content
    [maxl,maxlpos]=max(ls(j,:));
    l_map(3,j)=maxl;
    F_r_map(3,j)=fixed_F_rs(maxlpos);
    R_map(3,j)=Rs(j,maxlpos);
    disp(j)
end

%% Fit regulated ribosome content to growth law
% fit_reg = polyfit(l_map(1,:),R_map(1),1);
% fitted_reg = polyval(fit_reg,R_map(1));
% 
% fit_fix = polyfit(l_map(2,:),R_map(2),1);
% fitted_fix = polyval(fit_reg,R_map(2));
% 
% fit_opt = polyfit(l_map(3,:),R_map(3),1);
% fitted_opt = polyval(fit_reg,R_map(3));


%% Plot
Fig = figure('Position',[0 0 540 480]);
set(Fig, 'defaultLineLineWidth', 1)
set(Fig, 'defaultAxesFontSize', 12)
hold on
for j=1:size(nutrients,2)
    linecolour=0.75*(1-j/size(nutrients,2));
    plot(fixed_F_rs,ls(j,:),'Color',[linecolour, linecolour,linecolour])

    plot(F_r_map(1,j),l_map(1,j),'o','Color','r')
    plot(F_r_map(2,j),l_map(2,j),'o','Color','b')
    plot(F_r_map(3,j),l_map(3,j),'o','Color','g')
end

xlabel('Ribosome transcription activation (-)');
ylabel('SS growth rate (/h)');
title('Growth dependence on ribosome regulation');
hold off

Fig2 = figure('Position',[0 0 540 480]);
set(Fig2, 'defaultLineLineWidth', 2)
set(Fig2, 'defaultAxesFontSize', 16)
hold on
for j=1:size(nutrients,2)
    linecolour=0.75*(1-j/size(nutrients,2));
    plot(l_map(1,j),R_map(1,j),'x','Color','r')
    plot(l_map(2,j),R_map(2,j),'x','Color','b')
    plot(l_map(3,j),R_map(3,j),'x','Color','g')
end

xlabel('SS growth rate (/h)');
ylabel('Ribosome concentration (nM)');
title('Growth dependence on ribosome regulation');
hold off

%% Function for getting steady state
function [l, F_r, R]=get_steady_for_sweep(sim,Delta,Max_iter)
    % iterate integrations until steady state reached
    lRBm_old=zeros(1,3); % start with growth rate assumed to be zero
    for i=1:Max_iter
        % evaluate growth rate
        sim = sim.simulate_model;
           
        % evaluate growth
        par=sim.parameters;
        m_a = sim.x(end,1);
        m_r = sim.x(end,2);
        p_a = sim.x(end,3);
        R = sim.x(end,4);
        tc = sim.x(end,5);
        tu = sim.x(end,6);
        Bm = sim.x(end,7);

        e=sim.coll.e(par,tc);
        kcmh=par('kcm').*par('h');
        k_a=sim.coll.k(e,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
        k_r=sim.coll.k(e,par('k+_r'),par('k-_r'),par('n_r'),kcmh);
        D=1+(m_a./k_a+m_r./k_r)./(1-par('phi_q')); % denominator in ribosome competition cobj.optalculations
        B=R.*(1-1./D);

        l=sim.coll.l(par,e,B); % growth rate!
        
        T = tc./tu;
        F_r = sim.coll.F_r(par,T) ; % ribosome regulation
    
        if(norm([l,R,Bm]-lRBm_old)<Delta) % if SS reached, quit
            break
        else % if not, continue integrating
            lRBm_old=[l,R,Bm];
            sim.init_conditions('m_a')=m_a;
            sim.init_conditions('m_r')=m_r;
            % proteins
            sim.init_conditions('p_a')=p_a;
            sim.init_conditions('R')=R;
            % tRNAs
            sim.init_conditions('tc')=tc;
            sim.init_conditions('tu')=tu;
        end
    end
    if(i==Max_iter)
        disp('Warning! SS not reached yet')
%     else
%         disp(['SS reached in ',num2str(i),' iterations'])
    end
end