%% CLEAR all variables

addpath(genpath('..'))

close all
clear all

%% LOAD experimental data (55 years of E.coli measurements!)
% rib fractions vs growth rates
[dataset,captions,~] = xlsread('data/55_yrs_heterologous.csv');
ydata_xtra(:,1) = dataset(:,1); % het. prot. mass fraction
ydata_xtra(:,2) = dataset(:,3); % (growth rate w/ het. prot):(growth rate w/out)

%% SET UP the simulators
Delta=0.1; % difference in the state not considered important
Max_iter=100; % max. number of iterations
tf=20; % time frame of single integration (h)

sim={}; % simulators for different conditions will be stored here

%% DEFINE conditions
num_conds=1; % number of conditions to check
for cond=1:num_conds
    sim=advanced_simulator; % initialise simulation
    sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9);
    
    % amend default parameters
    a_total=300;
    if(cond==1)
        sim=sim.load_heterologous_and_external('two_constit','no_ext');
        sim.het.parameters('a_xtra1')=200;
        sim.het.parameters('a_xtra2')=100;
    end
    sim=sim.push_het;
end

%% SET UP the approximate estimator
approx=heterologous_approx;

%% SPECIFY which heterologous gene expression ranges you consider
plasmid_concs = logspace(0,log10(2500),50); % concentrations of plasmid w/ heterologous genes

% for different plasmid concs, we will get and record (for 2 different conditions)...
sss = {}; % system steady states
for cond=1:num_conds
    sss=zeros(size(sim.init_conditions,1),size(plasmid_concs,2));
end
ls = zeros(1,size(plasmid_concs,2)); % growth rates
es = zeros(1,size(plasmid_concs,2)); % elongation rates
phi_hets = zeros(1,size(plasmid_concs,2)); % heterologous protein fractions
t_sss = zeros(1,size(plasmid_concs,2)); % times in which steady state reached 
phi_xtra1s = zeros(1,size(plasmid_concs,2)); % heterologous protein fractions

% comparing them with values without het. prot. exp.
ss0=zeros((size(sim.init_conditions,1)+size(sim.het.init_conditions,1)),1);
l0=0;
e0=0;
phi_het0=0;
t_ss0=0;
phi_xtra10=0;

% approximate estimates
l0_approx=0;
l0_full_approx=0;
phi_het0_approx=0;
phi_xtra10_approx=0;
ls_approx = zeros(1,size(plasmid_concs,2)); % growth rates
phi_hets_approx = zeros(1,size(plasmid_concs,2)); % heterologous protein fractions
phi_xtra1s_approx = zeros(1,size(plasmid_concs,2)); % output protein fractions

%% RUN the model

% values w/out het. prot.
for i=1:sim.num_het
    sim.het.parameters(['c_',sim.het.names{i}])=0;
end

sim = sim.push_het(); % push initial condition to main object
sstss=get_steady(sim,Delta,Max_iter); ss=sstss{1}; t_ss=sstss{2}; % integrate + unpack outcome

[l,e,phi_het,phi_xtra1]=get_lephihetphixtra1(sim,ss); % get desired values
% record
ss0 = ss;
l0=l;
e0=e;
phi_het0=phi_het;
t_ss0=t_ss;
phi_xtra10=phi_xtra1;

% get approximate estimates
l0_approx=approx.ss_l(0,ss0,ss0,e0,sim);
l0_full_approx=approx.ss_l_full(0,ss0,ss0,e0,sim);
ls_approx=approx.ss_l(plasmid_concs,sss,ss0,e0,sim);
ls_full_approx=approx.ss_l_full(plasmid_concs,sss,ss0,e0,sim);
phi_hets_approx=approx.ss_phi_het(plasmid_concs,sss,ss0,e0,sim);

% find for xtra1
phi_xtra1s_approx=approx.ss_phi(plasmid_concs,sss, ... % steady state of the system
            1, ... % gene of interest (GOI) ss regulation function value
            plasmid_concs, ... % GOI concentration
            sim.parameters('a_xtra1'), ... % GOI transcription rate
            sim.parameters('k+_xtra1'), ... % GOI mRNA-ribosome binding rate
            sim.parameters('k-_xtra1'), ... % GOI mRNA-ribosome unbinding rate
            sim.parameters('n_xtra1'), ...  % GOI length in aa
            ss0, ... % steady state without heterologous gene expression
            e0, ... % steady state translation elongation rate without heterologous gene expression
            sim); % simulator (required to get parameters and regulatory functions for all genes)

% get actual model predictions
for i=1:size(plasmid_concs,2)
    % disp(plasmid_concs(i))
    for j=1:sim.num_het
        sim.het.parameters(['c_',sim.het.names{j}])=plasmid_concs(i);
    end
    sim = sim.push_het(); % reset initial condition
    sstss=get_steady(sim,Delta,Max_iter); ss=sstss{1}; t_ss=sstss{2}; % integrate + unpack outcome
    [l,e,phi_het,phi_xtra1]=get_lephihetphixtra1(sim,ss); % get desired values

    % record
    sss(:,i) = ss;
    ls(i)=l;
    es(i)=e;
    phi_hets(i)=phi_het;
    phi_xtra1s(i)=phi_xtra1;
end

%% PLOT: het. prot. mass fraction as a function of plasmid conc.
colours={'r', 'g'};
colours_approx={'#A2142F', '#77AC30'}; % colours for plotting approximations

Fig3c = figure('Position',[0 0 360 320]);
set(Fig3c, 'defaultAxesFontSize', 9)
set(Fig3c, 'defaultLineLineWidth', 1)
hold on

% plot overall het. prot. mass frac.
plot([0,plasmid_concs],[phi_het0,phi_hets],[colours{1},'-']) % plot model predictions
plot([0,plasmid_concs],[phi_het0_approx,phi_hets_approx], ...
    'Color',colours_approx{1},'LineStyle','--') % plot approximation

% plot output prot.mass frac.
plot([0,plasmid_concs],[phi_xtra10,phi_xtra1s],[colours{2},'-']) % plot model predictions
plot([0,plasmid_concs],[phi_xtra10_approx,phi_xtra1s_approx], ...
    'Color',colours_approx{2},'LineStyle','--') % plot approximation

grid on
xlim([0 2500])

xlabel('c, plasmid concentration [nM]');
ylabel('\phi_x, protein mass fraction');

legend('Total het. protein','Total het. protein (approx.)',...
    'Output protein','Output protein (approx.)','Location','southeast')
hold off

%% Function for getting desired steady state values
function [l,e,phi_het,phi_xtra1]=get_lephihetphixtra1(sim,ss)
    % PAREMETER MAP TO VARIABLE PAR
    par = sim.parameters;
    
    % STATE VECTOR TO SINGLE VARIABLES
    m_a = ss(1);
    m_r = ss(2);
    p_a = ss(3);
    R = ss(4);
    tc = ss(5);
    tu = ss(6);
    Bcm = ss(7);
    s = ss(8);
    h = ss(9);
    ss_het=ss(10 : (9+2*sim.num_het) ).';

    % USEFUL PRE-CALCULATIONS
    % translation elongation rate
    e=sim.coll.e(par,tc);

    % ribosome inactivation rate due to chloramphenicol
    kcmh=par('kcm').*h;

    % ribosome dissociation constants
    k_a=sim.coll.k(e,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
    k_r=sim.coll.k(e,par('k+_r'),par('k-_r'),par('n_r'),kcmh);
    % heterologous genes
    if(sim.num_het>0)
        k_het=zeros(sim.num_het,1);
        for i=1:sim.num_het
            k_het(i)=sim.coll.k(e,...
            sim.parameters(['k+_',sim.het.names{i}]),...
            sim.parameters(['k-_',sim.het.names{i}]),...
            sim.parameters(['n_',sim.het.names{i}]),...
            kcmh);
        end
    end

    T=tc./tu; % ratio of charged to uncharged tRNAs
    D=1+(m_a./k_a+m_r./k_r+sum(ss_het(1:sim.num_het)./k_het))./...
        (1-par('phi_q')); % denominator in ribosome competition csim.optalculations
    B=R.*(1-1./D); % actively translating ribosomes - INCLUDING Q

    % tRNA charging and synthesis
    nu=sim.coll.nu(par,tu,s);
    psi=sim.coll.psi(par,T);

    % growth rate
    l=sim.coll.l(par,e,B);

    % heterologous protein mass fraction
    phi_het=0;
    for i=1:sim.num_het
        phi_het=phi_het+ss_het(sim.num_het+i).*sim.parameters(['n_',sim.het.names{i}])./sim.parameters('M');
    end

    % output protein mass fraction
    phi_xtra1=ss_het(sim.num_het+1).*sim.parameters(['n_',sim.het.names{1}])./sim.parameters('M');
end
