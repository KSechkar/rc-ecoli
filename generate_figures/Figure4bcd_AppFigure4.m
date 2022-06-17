%% Figure4bcd_AppFigure4.m
% Analyse the effects of heterologous gene expression: investigate growth
% rate and het. prot. mas fraction (Figure 4), as well as charged/uncharged 
% tRNA conc. ratio, traslation el. rate, tRNA aminoacylation rate and  
% ribosomal genetranscription regulation (Appendix Figure 4)as a function of 
% het. gene transcription. Determine which heterologous protein mass fraction 
% leads to the highest protein production rates in a population of bacterai
% when the culturing starts. 

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

a_xtra=1000; % transcription rate of the heterologous gene PER PLASMID

sim=sim.load_heterologous_and_external('one_constit','no_ext'); % load heterologous gene expression module
sim.het.parameters('a_xtra')=a_xtra;
sim=sim.push_het;

plasmid_concs = logspace(0,3,100); % range of plasmid concentrations

%% SET UP the approximate estimator
approx=heterologous_approx;

%% INITIALISE arrays where obtained values of phys. variables will be stored

% values without het. prot. exp.
ss0=zeros((size(sim.init_conditions,1)+size(sim.het.init_conditions,1)),1); % steady states of the system
l0=0; % growth rates
e0=0;  % elongation rates
phi_het0=0; % heterologous protein fractions

% numerically obtained values
sss=zeros(size(sim.init_conditions,1),size(plasmid_concs,2)); % steady states of the system
ls = zeros(1,size(plasmid_concs,2)); % growth rates
es = zeros(1,size(plasmid_concs,2)); % elongation rates
phi_hets = zeros(1,size(plasmid_concs,2)); % heterologous protein fractions

% approximate estimates
l0_approx=0;
phi_het0_approx=0;
ls_approx = zeros(1,size(plasmid_concs,2)); % growth rates
phi_hets_approx = zeros(1,size(plasmid_concs,2)); % heterologous protein fractions

%% RUN withput heterologous gene expression
% values w/out het. prot.
for i=1:sim.num_het
    sim.het.parameters(['c_',sim.het.names{i}])=0;
end

sim = sim.push_het(); % push initial condition to main object
ss=get_steady(sim,Delta,Max_iter);

[l,e,phi_het]=get_lephihet(sim,ss); % get values of variables

% record
ss0 = ss;
l0=l;
e0=e;
phi_het0=phi_het;

% find k_xtra^NB (mRNA-ribosome dissoc. const for het. gene - will be needed later)
par=sim.het.parameters;
kxNB=sim.form.k(e,par('k+_xtra'),par('k-_xtra'),par('n_xtra'),0);

%% GET approximate estimates
l0_approx=approx.ss_l(0,ss0,ss0,e0,sim);
ls_approx=approx.ss_l(plasmid_concs,sss,ss0,e0,sim);
phi_hets_approx=approx.ss_phi_het(plasmid_concs,sss,ss0,e0,sim);

%% OBTAIN numerical predictions

% get actual model predictions
for i=1:size(plasmid_concs,2)
    % disp(plasmid_concs(i))
    for j=1:sim.num_het
        sim.het.parameters(['c_',sim.het.names{j}])=plasmid_concs(i);
    end
    sim = sim.push_het(); % reset initial condition
    ss=get_steady(sim,Delta,Max_iter);
    [l,e,phi_het]=get_lephihet(sim,ss); % get desired values

    % record
    sss(:,i) = ss;
    ls(i)=l;
    es(i)=e;
    phi_hets(i)=phi_het;
end

%% Figure 4b: cell growth as a function of burden

Fig4b = figure('Position',[0 0 316 280]);
set(Fig4b, 'defaultAxesFontSize', 9)
set(Fig4b, 'defaultLineLineWidth', 1)
hold on

plot([0,plasmid_concs*a_xtra]./kxNB,[l0,ls],['r','-']) % plot model predictions
plot([0,plasmid_concs*a_xtra]./kxNB,[l0_approx,ls_approx],...
    'Color',approx_colour,'LineStyle','--') % plot approximation

xlabel({'\xi, translational burden [1/h]'});
ylabel('\lambda, growth rate [1/h]');

grid on
hold off

%% Figure 4c: het. prot. mass fraction as a function of burden

Fig4c = figure('Position',[0 0 316 280]);
set(Fig4c, 'defaultAxesFontSize', 9)
set(Fig4c, 'defaultLineLineWidth', 1)

hold on

plot([0,plasmid_concs*a_xtra]./kxNB,[phi_het0,phi_hets],['r','-']) % plot model predictions
plot([0,plasmid_concs*a_xtra]./kxNB,[phi_het0_approx,phi_hets_approx], ...
    'Color',approx_colour,'LineStyle','--') % plot approximation

xlabel('\xi, translational burden [1/h]');
ylabel('\phi_x, het. prot. mass fraction');
grid on
hold off

%% Figure 4d - total heterologous protein production rate at t=0 in a population of cells (mu_het)
delta=0.25;
Fig4d = figure('Position',[0 0 345 288]);
hold on

set(Fig4d, 'defaultAxesFontSize', 9)
set(Fig4d, 'defaultLineLineWidth', 1)

mu_het0=0; % initialise
mu_het0_approx=0; % initialise


% make delta the same dimension as growth rate
deltas=delta*ones(size(ls));

% actual value of mu_het
mu_hets=(ls-deltas).*phi_hets.*sim.parameters('M');

% approximation of mu_het
mu_hets_approx=approx.ss_mu_het(plasmid_concs,sss,ss0,e0,l0,sim,delta);

plot([0,phi_hets],[mu_het0,mu_hets],['r','-']) % plot model predictions
plot([0,phi_hets_approx],[mu_het0_approx,mu_hets_approx], ...
    'Color',approx_colour,'LineStyle','--') % plot approximation

% draw a line to mark phi_het that maximises production (NUMERICAL RESULT)
[~,index_of_max_mu]=max(mu_hets); % get index
phi_het_max=phi_hets(index_of_max_mu); % get the value
mu_max_approx=mu_hets(index_of_max_mu); % find the corresponding mu
% get the value
plot([phi_het_max phi_het_max],[0 mu_max_approx], ...
    'Color',[0,0.8,0],'LineStyle','-','LineWidth',0.75) % plot

% draw a line to mark phi_het that maximises production  (ANALYTICAL PREDICTION)
phi_het_max_approx=approx.phi_het_max(l0,sim,delta); %calculate
mu_max_approx=sim.parameters('M').*phi_het_max_approx.*...
    (l0.*(1-phi_het_max_approx./(1-sim.parameters('phi_q')))-delta); % find corresponding mu
plot([phi_het_max_approx phi_het_max_approx],[0 mu_max_approx], ...
    'Color',[0,0.8,0],'LineStyle','--','LineWidth',0.75) % plot

xlabel('\phi_{poi}, prot. of interest mass fraction');
ylabel({'\mu, total prot. production','rate constant [aa/(h \cdot cell)]'});

xlim([0 0.25])
ylim([0 15*10^7])

yticks((0:2.5:15)*10^7)

grid on
hold off 

%% Appendix Figure 4 - metabolic changes induced by burden

Fig_app2 = figure('Position',[0 0 616 560]);

% T
subplot(2,2,3)
hold on

T0=ss0(5)/ss0(6);
Ts=sss(5,:)./sss(6,:);

plot([0,plasmid_concs*a_xtra],100*[T0,Ts]/T0,['r','-'])


xlabel('c_x\alpha_x, het. gene transcription rate [nM/h]');

ylabel({'T, proportional to 1/ppGpp conc.','[% of value with no burden]'});
ylim([85 115])


yticks(85:5:115)
xticks(0:2*10^5:10^6)
grid on
grid minor
hold off

% Elongation rates
subplot(2,2,1)
hold on

plot([0,plasmid_concs*a_xtra],100*[e0,es]/e0,['r','-']) % plot model predictions

xlabel('c_x\alpha_x, het. gene transcription rate [nM/h]');

ylabel({'\epsilon, translation elongation rate','[% of value with no burden]'});

yticks(85:5:115)
xticks(0:2*10^5:10^6)
grid on
grid minor
hold off

% Ribosome/tRNA transcription
subplot(2,2,2)
hold on

F_r0=sim.form.F_r(sim.parameters,T0);
F_rs=sim.form.F_r(sim.parameters,Ts);

plot([0,plasmid_concs*a_xtra],100*[F_r0,F_rs]/F_r0,['r','-'])

xlabel('c_x\alpha_x, het. gene transcription rate [nM/h]');
ylabel({'F_r, rib. gene transc. reg. func.','[% of value with no burden]'});
ylim([85 115])


yticks(85:5:115)
xticks(0:2*10^5:10^6)
grid on
grid minor
hold off

% tRNA charging
subplot(2,2,4)
hold on

nu0=sim.form.nu(sim.parameters,ss0(6),ss0(8));
nus=sim.form.nu(sim.parameters,sss(6,:),sss(8));

plot([0,plasmid_concs*a_xtra],100*[nu0,nus]/nu0,['r','-'])

xlabel('c_x\alpha_x, het. gene transcription rate [nM/h]');
ylabel({'\nu, tRNA charging rate.','[% of value with no burden]'});
ylim([85 115])

yticks(85:5:115)
xticks(0:2*10^5:10^6)
grid on
grid minor
hold off

%% FUNCTION for getting growth rate, translation elongation rate and het. prot. mass fraction from the system's steady state
function [l,e,phi_het]=get_lephihet(sim,ss)
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
end
