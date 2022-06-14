%% Initial preparations

addpath(genpath('..'))

Delta=0.01; % difference in growth rate not considered important
Max_iter=100; % max. number of iterations
sim = advanced_simulator;
sim.tf=100; % time frame of single integration (h)
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9);

close all

nutrients=flip([0.05,0.1,0.25,0.5,0.75,1]);
% storing growth rates
l_map = zeros(3,size(nutrients,2));

% storing ribosome transcription rates
F_r_map = zeros(3,size(nutrients,2));

% storing ribosome concentrations
R_map = zeros(3,size(nutrients,2));

fixed_F_rs= 0.01:0.01:1; % range of fixed ribosome transcription coefficients
ls=zeros(size(nutrients,2),size(fixed_F_rs,2));
Rs=zeros(size(nutrients,2),size(fixed_F_rs,2));

for j=1:size(nutrients,2)
    sim=sim.set_default_parameters(); % reset initial consitions and parameters
    sim=sim.set_default_init_conditions(); % reset initial consitions and parameters
    sim.init_conditions('s')=nutrients(j); % set nutrient quality

    % Run with regulated ribosome transcription
    sim.parameters('is_fixed_F_r')=0; % F_r regulated now
    [l, F_r, R] = get_steady_for_sweep(sim,Delta,Max_iter);
    l_map(1,j)=l;
    F_r_map(1,j)=F_r;
    R_map(1,j)=R;

    sim=sim.set_default_parameters(); % reset initial consitions and parameters
    sim=sim.set_default_init_conditions(); % reset initial consitions and parameters
    sim.init_conditions('s')=nutrients(j); % set nutrient quality    
    
    % Run with fixed ribosome transcription
    sim.parameters('is_fixed_F_r')=1; % F_r fixed now
    sim.parameters('fixed_F_r')=0.6;
    [l, F_r, R] = get_steady_for_sweep(sim,Delta,Max_iter);
    l_map(2,j)=l;
    F_r_map(2,j)=F_r;
    R_map(2,j)=R;
    
    sim=sim.set_default_parameters(); % reset initial consitions and parameters
    sim=sim.set_default_init_conditions(); % reset initial consitions and parameters
    sim.init_conditions('s')=nutrients(j); % set nutrient quality

    % Run with a range of fixed ribosome transcriptions
    for i=1:size(fixed_F_rs,2)
        sim=sim.set_default_parameters(); % reset initial consitions and parameters
        sim=sim.set_default_init_conditions(); % reset initial consitions and parameters
        sim.init_conditions('s')=nutrients(j); % set nutrient quality
        sim.parameters('is_fixed_F_r')=1; % F_r fixed now
        
        disp(fixed_F_rs(i))
        sim.parameters('fixed_F_r')=fixed_F_rs(i);
        [l, F_r, R] = get_steady_for_sweep(sim,Delta,Max_iter);
        ls(j,i)=l;
        Rs(j,i)=R;
    end
    % finding optimal ribosome content
    [maxl,maxlpos]=max(ls(j,:));
    l_map(3,j)=maxl;
    F_r_map(3,j)=fixed_F_rs(maxlpos);
    R_map(3,j)=Rs(j,maxlpos);
    disp(j)
end

%% Plot
Fig = figure('Position',[0 0 316 240]);
set(Fig, 'defaultAxesFontSize', 9)
set(Fig, 'defaultLineLineWidth', 1)
hold on
for j=1:size(nutrients,2)
    plot(F_r_map(1,j),l_map(1,j),'o','Color','r')
    plot(F_r_map(3,j),l_map(3,j),'o','Color',[0,0.75,0])

    linecolour=0.75*(1-nutrients(j));
    plot(fixed_F_rs,ls(j,:),'Color',[linecolour, linecolour,linecolour])
end

xlabel({'F_r, ribosomal gene transcription','regulation func.'});
ylabel('Growth rate \lambda, 1/h');
ylim([0 2.5])
grid on
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