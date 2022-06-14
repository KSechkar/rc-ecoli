%% CLEAR

addpath(genpath('..'))

clear all
close all

%% LOAD 55 years of E.coli measurements!
% rib fractions vs growth rates
[dataset,captions,~] = xlsread('data/55_yrs.csv');
ydata_rib=zeros(size(dataset,1),2); % intialise outputs array
for i = 1:size(dataset,1)
    % outputs
    ydata_rib(i,1) = dataset(i,1); % growth rate (1/h)
    ydata_rib(i,2) = dataset(i,2); % ribosome mass fraction
end
% disp(ydata_rib)

% translation elongation rates vs growth rates
[dataset,captions,~] = xlsread('data/55_yrs_elongation.csv');
ydata_el=zeros(size(dataset,1),2); % intialise outputs array
for i = 1:size(dataset,1)
    % outputs
    ydata_el(i,1) = dataset(i,1); % growth rate (1/h)
    ydata_el(i,2) = dataset(i,2); % translation elongation rate
end
% disp(ydata_el)

%% SET UP SIMULATOR
% params for finding steady state
Delta = 0.01; % threshold below which the changes in l and phi_r assumed negligible
Max_iter=100; % max. no iterations for finding steady state
tf=5;
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9);
sim=advanced_simulator; % initialise simulator

%% GET model predictions
% set up simulator

nutr_quals=logspace(log10(0.01),log10(1),8); % get a wide range of nutrient qualities

l_map=zeros(size(nutr_quals));
phir_map=zeros(size(nutr_quals));
el_map=zeros(size(nutr_quals));
for j=1:size(nutr_quals,2)
    disp(nutr_quals(j))
    % reset simulator
    sim.parameters=advanced_params();
    sim.init_conditions=advanced_init_conds(sim.parameters); % reset intial conditions
    
    % set nutrient quality
    sim.init_conditions('s')=nutr_quals(j);
    sim=sim.push_het();

    % Run
    sim.parameters('is_fixed_F_r')=0; % F_r regulated now
    [l, R, e] = get_steady_for_55(sim,Delta,Max_iter);
    l_map(j)=l;
    phir_map(j)=R.*sim.parameters('n_r')./sim.parameters('M');
    el_map(j)=e./3600;
end

%% PLOT and compare
Frib = figure('Position',[0 0 316 240]);
set(Frib, 'defaultAxesFontSize', 9)
hold on
plot(ydata_rib(:,1),ydata_rib(:,2),'.','Color','b')
plot(l_map,phir_map,'-','Color','r')
plot([0.5 0.5],[0 0.3],':','Color','k')
ylabel('\phi_r, ribosomal mass fraction');
xlabel('\lambda, growth rate [1/h]');
legend('Experimental data','Model predictions','Location','southeast')
grid on
hold off

Fel = figure('Position',[0 0 316 240]);
set(Fel, 'defaultAxesFontSize', 9)
hold on
plot(ydata_el(:,1),ydata_el(:,2),'.','Color','b')
plot(l_map,el_map,'-','Color','r')
plot([0.5 0.5],[6 20],':','Color','k')
ylabel('\epsilon, translation elong. rate [aa/s]');
xlabel('\lambda, growth rate [1/h]');
legend('Experimental data','Model predictions','Location','southeast')
grid on
hold off

%% Function for getting steady states
function [l, R, e]=get_steady_for_55(sim,Delta,Max_iter)
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
        Bcm = sim.x(end,7);

        e=sim.coll.e(par,tc);
        kcmh=par('kcm').*par('h');
        k_a=sim.coll.k(e,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
        k_r=sim.coll.k(e,par('k+_r'),par('k-_r'),par('n_r'),kcmh);
        D=1+(m_a./k_a+m_r./k_r)./(1-par('phi_q')); % denominator in ribosome competition cobj.optalculations
        B=R.*(1-1./D);

        l=sim.coll.l(par,e,B); % growth rate!
    
        if(norm([l,R,Bcm]-lRBm_old)<Delta) % if SS reached, quit
            break
        else % if not, continue integrating
            lRBm_old=[l,R,Bcm];
            sim.init_conditions('m_a')=m_a;
            sim.init_conditions('m_r')=m_r;
            % proteins
            sim.init_conditions('p_a')=p_a;
            sim.init_conditions('R')=R;
            % tRNAs
            sim.init_conditions('tc')=tc;
            sim.init_conditions('tu')=tu;
            % inactivated ribosomes
            sim.init_conditions('Bcm')=Bcm;
        end
    end
    if(i==Max_iter)
        disp('Warning! SS not reached yet')
    end
end
