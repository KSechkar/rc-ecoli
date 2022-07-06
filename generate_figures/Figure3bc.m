%% Figure3bc
% Compare model predictions of the cell's growth rate, ribosomal mass
% fraction and translation elongation rate with experimental data from the 
% last 55 years of measurements, compiled by Chure et al. 2022

%% CLEAR

addpath(genpath('..'))

clear all
close all

%% LOAD experimental data (55 years of E.coli measurements!)
% ribosomal mass fractions vs growth rates
[dataset,captions,~] = xlsread('data/55_yrs_ribosomes.csv');
data_rib(:,1) = dataset(:,1); % het. prot. mass fraction
data_rib(:,2) = dataset(:,2); % (growth rate w/ het. prot):(growth rate w/out)

% translation elongation rates vs growth rates
[dataset,captions,~] = xlsread('data/55_yrs_elongation.csv');
data_el(:,1) = dataset(:,1); % growth rate (1/h)
data_el(:,2) = dataset(:,2); % translation elongation rate

%% SET UP the simulator

sim=cell_simulator; % initialise simulator

% parameters for getting steady state
sim.tf = 5; % single integraton step timeframe
Delta = 0.01; % threshold that determines if we're in steady state
Max_iter = 100; % maximum no. iterations (checking if SS reached over first 750 h)

sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9); % more lenient integration tolerances for speed

%% SPECIFY range of nutrient qualit for which we run the simulation

nutr_quals=logspace(log10(0.01),log10(1),8);

%% INITIALISE the arrays where model predictions will be stored

l_map=zeros(size(nutr_quals));
phir_map=zeros(size(nutr_quals));
el_map=zeros(size(nutr_quals));

%% GET model predictions

for j=1:size(nutr_quals,2)
    disp(nutr_quals(j))

    % reset simulator
    sim.parameters=cell_params();
    sim.init_conditions=cell_init_conds(sim.parameters); % reset intial conditions
    
    % set nutrient quality
    sim.init_conditions('s')=nutr_quals(j);
    sim=sim.push_het();

    % Run
    sim.parameters('is_fixed_F_r')=0; % F_r regulated now
    ss = get_steady(sim,Delta,Max_iter);
    [l, e, phi_r] = get_lephir(sim,ss);
    l_map(j)=l;
    phir_map(j)=phi_r;
    el_map(j)=e./3600; % record in aa/s, NOT aa/h
end

%% Figure 3b - growth rates vs rib. mass fractions

Frib = figure('Position',[0 0 316 240]);
set(Frib, 'defaultAxesFontSize', 9)
hold on
plot(data_rib(:,1),data_rib(:,2),'.','Color','b')
plot(l_map,phir_map,'-','Color','r')
plot([0.5 0.5],[0 0.3],':','Color','k')
ylabel('\phi_r, ribosomal mass fraction');
xlabel('\lambda, growth rate [1/h]');
legend('Experimental data','Model predictions','Location','southeast')
grid on
hold off

%% Figure 3c - growth rates vs translation elongation rates

Fel = figure('Position',[0 0 316 240]);
set(Fel, 'defaultAxesFontSize', 9)
hold on
plot(data_el(:,1),data_el(:,2),'.','Color','b')
plot(l_map,el_map,'-','Color','r')
plot([0.5 0.5],[6 20],':','Color','k')
ylabel('\epsilon, translation elong. rate [aa/s]');
xlabel('\lambda, growth rate [1/h]');
legend('Experimental data','Model predictions','Location','southeast')
grid on
hold off

%% FUNCTION for getting growth rate, translation elongation rate and rib. mass fraction from the system's steady state

function [l, e, phi_r]=get_lephir(sim,ss)
    % evaluate growth
    par=sim.parameters;
    m_a = ss(1); % metabolic gene mRNA
    m_r = ss(2); % ribosomal mRNA
    p_a = ss(3); % metabolic proteins
    R = ss(4); % operational ribosomes
    tc = ss(5); % charged tRNAs
    tu = ss(6); % uncharged tRNAs
    Bm = ss(7); % inactivated ribosomes
    s=ss(8); % nutrient quality
    h=ss(9); % chloramphenicol conc.
    ss_het=ss(10:9+2*sim.num_het); % heterologous genes

    e=sim.form.e(par,tc); % translation elongation ratio
    k_a=sim.form.k(e,par('k+_a'),par('k-_a'),par('n_a'),par('kcm').*h); % ribosome-mRNA dissociation constant (metab. genes)
    k_r=sim.form.k(e,par('k+_r'),par('k-_r'),par('n_r'),par('kcm').*h); % ribosome-mRNA dissociation constant (rib. genes)
    D=1+(m_a./k_a+m_r./k_r)./(1-par('phi_q')); % denominator in ribosome competition calculations
    B=R.*(1-1./D); % actively translating ribosomes (inc. those translating housekeeping genes)

    l=sim.form.l(par,e,B); % growth rate!
    e=sim.form.e(par,tc); % translation elongation rate!

    % rib. mass fraction
    phi_r=R.*par('n_r')./par('M');
end
