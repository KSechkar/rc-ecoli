%% Figure6.m
% Analyse the two-gene circuit for plasmid copy number -indeoendent gene
% expression: investigate the output protein mass fraction a a function of
% plasmid concentration

%% CLEAR all variables

addpath(genpath('..'))

close all
clear all

%% SET UP the simulator

sim=cell_simulator; % initialise simulator

% parameters for getting steady state
sim.tf = 20; % single integraton step timeframe
Delta = 0.1; % threshold that determines if we're in steady state
Max_iter = 100; % maximum no. iterations (checking if SS reached over first 750 h)

sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9); % more lenient integration tolerances for speed


%% DEFINE gene expression paramaters

a_total=300; % total het. gene transcription rate per plasmid

sim=sim.load_heterologous_and_external('two_constit','no_ext'); % load the het. gene exp. module
sim.het.parameters('a_xtra1')=200; % outptut gene transcription rate
sim.het.parameters('a_xtra2')=100; % controller gene transcription rate

plasmid_concs = logspace(0,log10(2500),50); % concentrations of plasmid w/ heterologous genes

%% SET UP the approximate estimator
approx=heterologous_approx;

%% INITIALISE arrays where obtained values of phys. variables will be stored

% values without het. prot. exp.
ss0=zeros((size(sim.init_conditions,1)+size(sim.het.init_conditions,1)),1); % steady states of the system
l0=0; % growth rate
e0=0;  % elongation rate
phi_het0=0; % heterologous protein fraction
phi_xtra10 = zeros(1,size(plasmid_concs,2)); % output protein mass fraction

% numerically obtained values
sss=zeros(size(sim.init_conditions,1),size(plasmid_concs,2)); % steady states of the system
ls = zeros(1,size(plasmid_concs,2)); % growth rates
es = zeros(1,size(plasmid_concs,2)); % elongation rates
phi_hets = zeros(1,size(plasmid_concs,2)); % heterologous protein fractions
phi_xtra1s = zeros(1,size(plasmid_concs,2)); % output protein mass fractions

% approximate estimates
l0_approx=0;
phi_het0_approx=0;
phi_xtra10_approx=0;
ls_approx = zeros(1,size(plasmid_concs,2)); % growth rates
phi_hets_approx = zeros(1,size(plasmid_concs,2)); % total heterologous protein mass fractions
phi_xtra1s_approx = zeros(1,size(plasmid_concs,2)); % output protein mass fractions

%% RUN withput heterologous gene expression
% values w/out het. prot.
for i=1:sim.num_het
    sim.het.parameters(['c_',sim.het.names{i}])=0;
end

sim = sim.push_het(); % push initial condition to main object
ss=get_steady(sim,Delta,Max_iter);

[l,e,phi_het,phi_xtra1]=get_lephihetphixtra1(sim,ss); % get values of variables

% record
ss0 = ss;
l0=l;
e0=e;
phi_het0=phi_het;
phi_xtra10=0;

%% GET approximate estimates
l0_approx=approx.ss_l(0,ss0,ss0,e0,sim);
ls_approx=approx.ss_l(plasmid_concs,sss,ss0,e0,sim);
phi_hets_approx=approx.ss_phi_het(plasmid_concs,sss,ss0,e0,sim); % total

% output gene mass fractions
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

%% OBTAIN numerical predictions

% get actual model predictions
for i=1:size(plasmid_concs,2)
    % disp(plasmid_concs(i))
    for j=1:sim.num_het
        sim.het.parameters(['c_',sim.het.names{j}])=plasmid_concs(i);
    end
    sim = sim.push_het(); % reset initial condition
    ss=get_steady(sim,Delta,Max_iter);
    [l,e,phi_het,phi_xtra1]=get_lephihetphixtra1(sim,ss); % get desired values

    % record
    sss(:,i) = ss;
    ls(i)=l;
    es(i)=e;
    phi_hets(i)=phi_het;
    phi_xtra1s(i)=phi_xtra1;
end

%% Figure 6: plot het. prot. mass fraction as a function of plasmid conc.
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

%% FUNCTION for getting growth rate, translation elongation rate, het. prot. and output prot. mass fraction from the system's steady state
function [l,e,phi_het,phi_xtra1]=get_lephihetphixtra1(sim,ss)
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

    % heterologous protein mass fraction
    phi_het=0;
    for i=1:sim.num_het
        phi_het=phi_het+ss_het(sim.num_het+i).*sim.parameters(['n_',sim.het.names{i}])./sim.parameters('M');
    end

    % output protein mass fraction
    phi_xtra1=ss_het(sim.num_het+1).*sim.parameters(['n_',sim.het.names{1}])./sim.parameters('M');
end
