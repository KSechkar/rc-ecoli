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
    sim{cond}=advanced_simulator; % initialise simulation
    sim{cond}.opt = odeset('reltol',1.e-6,'abstol',1.e-9);
    
    % amend default parameters
    a_total=1000;
    if(cond==1)
        sim{cond}=sim{cond}.load_heterologous_and_external('one_constit','no_ext');
        sim{cond}.het.parameters('a_xtra')=a_total;
    end
    sim{cond}=sim{cond}.push_het;
end

%% SET UP the approximate estimator
approx=heterologous_approx;

%% SPECIFY which heterologous gene expression ranges you consider
plasmid_concs = logspace(0,3,100); % concentrations of plasmid w/ heterologous genes

% for different plasmid concs, we will get and record (for 2 different conditions)...
sss = {}; % system steady states
for cond=1:num_conds
    sss{cond}=zeros(size(sim{cond}.init_conditions,1),size(plasmid_concs,2));
end
ls = {zeros(1,size(plasmid_concs,2)), zeros(1,size(plasmid_concs,2))}; % growth rates
es = {zeros(1,size(plasmid_concs,2)), zeros(1,size(plasmid_concs,2))}; % elongation rates
phi_hets = {zeros(1,size(plasmid_concs,2)), zeros(1,size(plasmid_concs,2))}; % heterologous protein fractions
t_sss = {zeros(1,size(plasmid_concs,2)), zeros(1,size(plasmid_concs,2))}; % times in which steady state reached 

% comparing them with values without het. prot. exp.
ss0={zeros((size(sim{cond}.init_conditions,1)+size(sim{cond}.het.init_conditions,1)),1),...
    zeros((size(sim{cond}.init_conditions,1)+size(sim{cond}.het.init_conditions,1)),1)};
l0={0,0};
e0={0,0};
phi_het0={0,0};
t_ss0={0,0};

% approximate estimates
l0_approx={0,0};
l0_full_approx={0,0};
phi_het0_approx={0,0};
ls_approx = {zeros(1,size(plasmid_concs,2)), zeros(1,size(plasmid_concs,2))}; % growth rates
ls_full_approx = {zeros(1,size(plasmid_concs,2)), zeros(1,size(plasmid_concs,2))}; % growth rates with more accurate approximation
phi_hets_approx = {zeros(1,size(plasmid_concs,2)), zeros(1,size(plasmid_concs,2))}; % heterologous protein fractions

%% RUN for 2 conditions
for cond=1:num_conds
    disp(cond)
    % values w/out het. prot.
    for i=1:sim{cond}.num_het
        sim{cond}.het.parameters(['c_',sim{cond}.het.names{i}])=0;
    end

    sim{cond} = sim{cond}.push_het(); % push initial condition to main object
    sstss=get_steady(sim{cond},Delta,Max_iter); ss=sstss{1}; t_ss=sstss{2}; % integrate + unpack outcome
    
    [l,e,phi_het]=get_lephihet(sim{cond},ss); % get desired values
    % record
    ss0{cond} = ss;
    l0{cond}=l;
    e0{cond}=e;
    phi_het0{cond}=phi_het;
    t_ss0{cond}=t_ss;

    % find k_x^NB
    par=sim{cond}.het.parameters;
    kxNB=sim{cond}.coll.k(e,par('k+_xtra'),par('k-_xtra'),par('n_xtra'),0);

    % get approximate estimates
    l0_approx{cond}=approx.ss_l(0,ss0{cond},ss0{cond},e0{cond},sim{cond});
    l0_full_approx{cond}=approx.ss_l_full(0,ss0{cond},ss0{cond},e0{cond},sim{cond});
    ls_approx{cond}=approx.ss_l(plasmid_concs,sss{cond},ss0{cond},e0{cond},sim{cond});
    ls_full_approx{cond}=approx.ss_l_full(plasmid_concs,sss{cond},ss0{cond},e0{cond},sim{cond});
    phi_hets_approx{cond}=approx.ss_phi_het(plasmid_concs,sss{cond},ss0{cond},e0{cond},sim{cond});

    % get actual model predictions
    for i=1:size(plasmid_concs,2)
        % disp(plasmid_concs(i))
        for j=1:sim{cond}.num_het
            sim{cond}.het.parameters(['c_',sim{cond}.het.names{j}])=plasmid_concs(i);
        end
        sim{cond} = sim{cond}.push_het(); % reset initial condition
        sstss=get_steady(sim{cond},Delta,Max_iter); ss=sstss{1}; t_ss=sstss{2}; % integrate + unpack outcome
        [l,e,phi_het]=get_lephihet(sim{cond},ss); % get desired values
    
        % record
        sss{cond}(:,i) = ss;
        ls{cond}(i)=l;
        es{cond}(i)=e;
        phi_hets{cond}(i)=phi_het;
    end
end

%% Figure 4a: growth rate vs het. prot. mass frac.
colours={'r', 'g'};
colours_approx={'#A2142F', '#77AC30'}; % colours for plotting approximations
Fig = figure('Position',[0 0 316 280]);
set(Fig, 'defaultAxesFontSize', 9)
hold on
for cond=1:num_conds
    plot([0,phi_hets{cond}],[l0{cond},ls{cond}]/l0{cond},[colours{cond},'-'],'LineWidth',1) % plot model predictions
    plot([0,phi_hets_approx{cond}],[l0_approx{cond},ls_approx{cond}]/l0_approx{cond},...
        'Color',colours_approx{cond},'LineStyle','--','LineWidth',1) % plot approximate estimates
end
plot(ydata_xtra(:,1),ydata_xtra(:,2),'b.') % plot experimental data

plot([0.25 0.25],[0 1.05],':','Color','k')
xlabel('\phi_x, het. prot. mass fraction');
ylabel('\lambda:\lambda^{NB}, relative growth rate');
xlim([0 0.41])
ylim([0 1.05])

xticks([0:0.125:0.41,0.41])
yticks([0:0.2:1,1.05])
grid on
hold off

%% Figure 4b: cell growth as a function of burden
absolutes=true;

Fig4b = figure('Position',[0 0 316 280]);
set(Fig4b, 'defaultAxesFontSize', 9)
set(Fig4b, 'defaultLineLineWidth', 1)
hold on
for cond=1:num_conds
    if(absolutes)
        plot([0,plasmid_concs*a_total]./kxNB,[l0{cond},ls{cond}],[colours{cond},'-']) % plot model predictions
        plot([0,plasmid_concs*a_total]./kxNB,[l0_approx{cond},ls_approx{cond}],...
            'Color',colours_approx{cond},'LineStyle','--') % plot approximation
    else
        plot([0,plasmid_concs*a_total]./kxNB,[l0{cond},ls{cond}]/l0{cond},[colours{cond},'-']) % plot model predictions
        plot([0,plasmid_concs*a_total]./kxNB,[l0_approx{cond},ls_approx{cond}]/l0_approx{cond},...
            'Color',colours_approx{cond},'LineStyle','--') % plot approximation
    end
end
xlabel({'\xi, translational burden [1/h]'});
if(absolutes)
    ylabel('\lambda, growth rate [1/h]');
else
    ylabel('l with xtra:l without');
    ylim([0 1])
end
%title('Growth dependence on het. prot. exp.');
%xlim([0 0.2])
%ylim([0 1.05])
%legend('nutr. qual.=1','nutr. qual.=0.2')
grid on
hold off

%% Figure 4c: het. prot. mass fraction as a function of burden

Fig4c = figure('Position',[0 0 316 280]);
set(Fig4c, 'defaultAxesFontSize', 9)
set(Fig4c, 'defaultLineLineWidth', 1)

hold on
for cond=1:num_conds
    plot([0,plasmid_concs*a_total]./kxNB,[phi_het0{cond},phi_hets{cond}],[colours{cond},'-']) % plot model predictions
    plot([0,plasmid_concs*a_total]./kxNB,[phi_het0_approx{cond},phi_hets_approx{cond}], ...
        'Color',colours_approx{cond},'LineStyle','--') % plot approximation
end
xlabel('\xi, translational burden [1/h]');
ylabel('\phi_x, het. prot. mass fraction');
grid on
hold off

%% Figure 4d - total heterologous protein production rate
delta=0.25;
Fig4d = figure('Position',[0 0 345 288]);
hold on

set(Fig4d, 'defaultAxesFontSize', 9)
set(Fig4d, 'defaultLineLineWidth', 1)

mu_het0={0, 0}; % initialise
mu_hets={}; % initialise
mu_het0_approx={0, 0}; % initialise
mu_hets_approx={}; % initialise
phi_het_maxs={}; % initialise

for cond=1:num_conds
    % make delta the same dimension as growth rate
    deltas=delta*ones(size(ls{cond}));

    % actual value
    mu_hets{cond}=(ls{cond}-deltas).*phi_hets{cond}.*sim{cond}.parameters('M');

    % approximation
    mu_hets_approx{cond}=approx.ss_mu_het(plasmid_concs,sss{cond},ss0{cond},e0{cond},l0{cond},sim{cond},delta);

    plot([0,phi_hets{cond}],[mu_het0{cond},mu_hets{cond}],[colours{cond},'-']) % plot model predictions
    plot([0,phi_hets_approx{cond}],[mu_het0_approx{cond},mu_hets_approx{cond}], ...
        'Color',colours_approx{cond},'LineStyle','--') % plot approximation
    
    % draw a line for phi_het maximising production (NUMERICAL RESULT)
    [~,index_of_max_mu]=max(mu_hets{cond}); % get index
    phi_het_max{cond}=phi_hets{cond}(index_of_max_mu); % get the value
    mu_max_approx{cond}=mu_hets{cond}(index_of_max_mu); % find the corresponding mu
    % get the value
    plot([phi_het_max{cond} phi_het_max{cond}],[0 mu_max_approx{cond}], ...
        'Color',[0,0.8,0],'LineStyle','-','LineWidth',0.75) % plot

    % draw a line for phi_het maximising production (ANALYTICAL PREDICTION)
    phi_het_max_approx{cond}=0.5.*(1-delta./l0{cond}).*(1-sim{cond}.parameters('phi_q')); %calculate
    mu_max_approx{cond}=sim{cond}.parameters('M').*phi_het_max_approx{cond}.*...
        (l0{cond}.*(1-phi_het_max_approx{cond}./(1-sim{cond}.parameters('phi_q')))-delta); % find corresponding mu
    plot([phi_het_max_approx{cond} phi_het_max_approx{cond}],[0 mu_max_approx{cond}], ...
        'Color',[0,0.8,0],'LineStyle','--','LineWidth',0.75) % plot

end


xlabel('\phi_{poi}, prot. of interest mass fraction');
ylabel({'\mu, total prot. production','rate constant [aa/(h \cdot cell)]'});

xlim([0 0.25])
ylim([0 15*10^7])

yticks((0:2.5:15)*10^7)

grid on
hold off 

%% Appendix Figure 4 - metabolic changes induced by burden
absolutes = false;

Fig_app2 = figure('Position',[0 0 616 560]);

% T
subplot(2,2,3)
hold on
T0={};
Ts={};
for cond=1:num_conds
    T0{cond}=ss0{cond}(5)/ss0{cond}(6);
    Ts{cond}=sss{cond}(5,:)./sss{cond}(6,:);
    if(absolutes)
        plot([0,plasmid_concs*a_total],[T0{cond},Ts{cond}],[colours{cond},'-']) % plot model predictions
    else
        plot([0,plasmid_concs*a_total],100*[T0{cond},Ts{cond}]/T0{cond},[colours{cond},'-'])
    end
end
xlabel('c_x\alpha_x, het. gene transcription rate [nM/h]');
if(absolutes)
    ylabel('T');
else
    ylabel({'T, proportional to 1/ppGpp conc.','[% of value with no burden]'});
    ylim([85 115])
end

yticks(85:5:115)
xticks(0:2*10^5:10^6)
grid on
grid minor
hold off

% Elongation rates
subplot(2,2,1)
hold on
for cond=1:num_conds
    if(absolutes)
        plot([0,plasmid_concs*a_total],[e0{cond},es{cond}],[colours{cond},'-']) % plot model predictions
    else
        plot([0,plasmid_concs*a_total],100*[e0{cond},es{cond}]/e0{cond},[colours{cond},'-']) % plot model predictions
    end
end
xlabel('c_x\alpha_x, het. gene transcription rate [nM/h]');
if(absolutes)
    ylabel('e, /h');
else
    ylabel({'\epsilon, translation elongation rate','[% of value with no burden]'});
    ylim([85 115])
end

yticks(85:5:115)
xticks(0:2*10^5:10^6)
grid on
grid minor
hold off

% Ribosome/tRNA transcription
subplot(2,2,2)
hold on
F_r0={};
F_rs={};
for cond=1:num_conds
    F_r0{cond}=sim{cond}.coll.F_r(sim{cond}.parameters,T0{cond});
    F_rs{cond}=sim{cond}.coll.F_r(sim{cond}.parameters,Ts{cond});
    if(absolutes)
        plot([0,plasmid_concs*a_total],[F_r0{cond},F_rs{cond}],[colours{cond},'-']) % plot model predictions
    else
        plot([0,plasmid_concs*a_total],100*[F_r0{cond},F_rs{cond}]/F_r0{cond},[colours{cond},'-'])
    end
end
xlabel('c_x\alpha_x, het. gene transcription rate [nM/h]');
if(absolutes)
    ylabel('F_r');
else
    ylabel({'F_r, rib. gene transc. reg. func.','[% of value with no burden]'});
    ylim([85 115])
end

yticks(85:5:115)
xticks(0:2*10^5:10^6)
grid on
grid minor
hold off

% tRNA charging
subplot(2,2,4)
hold on
nu0={};
nus={};
for cond=1:num_conds
    nu0{cond}=sim{cond}.coll.nu(sim{cond}.parameters,ss0{cond}(6),ss0{cond}(8));
    nus{cond}=sim{cond}.coll.nu(sim{cond}.parameters,sss{cond}(6,:),sss{cond}(8));
    if(absolutes)
        plot([0,plasmid_concs*a_total],[nu0{cond},nus{cond}],[colours{cond},'-']) % plot model predictions
    else
        plot([0,plasmid_concs*a_total],100*[nu0{cond},nus{cond}]/nu0{cond},[colours{cond},'-'])
    end
end
xlabel('c_x\alpha_x, het. gene transcription rate [nM/h]');
if(absolutes)
    ylabel('phi_r');
else
    ylabel({'\nu, tRNA charging rate.','[% of value with no burden]'});
    ylim([85 115])
end

yticks(85:5:115)
xticks(0:2*10^5:10^6)
grid on
grid minor
hold off

%% Function for getting desired steady state values
function [l,e,phi_het]=get_lephihet(sim,ss)
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
end
