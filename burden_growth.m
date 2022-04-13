%% CLEAR all variables
close all
clear all

%% LOAD experimental data (55 years of E.coli measurements!)
% rib fractions vs growth rates
[dataset,captions,~] = xlsread('data/55_yrs_heterologous.csv');
ydata_xtra(:,1) = dataset(:,1); % het. prot. mass fraction
ydata_xtra(:,2) = dataset(:,3); % (growth rate w/ het. prot):(growth rate w/out)

%% SET UP the simulator
Delta=0.1; % difference in the state not considered important
Max_iter=100; % max. number of iterations
sim = advanced_simulator;
sim.tf=20; % time frame of single integration (h)
sim.parameters('nutr_qual')=0.2; % nutrient conc.

%% SET UP the approximate estimator
approx=heterologous_approx;

%% SPECIFY which heterologous gene expression ranges you consider
plasmid_concs = logspace(0,1.6,25); % concentrations of plasmid w/ heterologous genes

% for different plasmid concs, we will get and record (for 2 different conditions)...
sss = {zeros((size(sim.init_conditions,1)+size(sim.het.init_conditions,1)),size(plasmid_concs,2)),...
    zeros((size(sim.init_conditions,1)+size(sim.het.init_conditions,1)),size(plasmid_concs,2))}; % system steady states
ls = {zeros(1,size(plasmid_concs,2)), zeros(1,size(plasmid_concs,2))}; % growth rates
es = {zeros(1,size(plasmid_concs,2)), zeros(1,size(plasmid_concs,2))}; % elongation rates
phi_hets = {zeros(1,size(plasmid_concs,2)), zeros(1,size(plasmid_concs,2))}; % heterologous protein fractions

% comparing them with values without het. prot. exp.
ss0={zeros((size(sim.init_conditions,1)+size(sim.het.init_conditions,1)),1),...
    zeros((size(sim.init_conditions,1)+size(sim.het.init_conditions,1)),1)};
l0={0,0};
e0={0,0};
phi_het0={0,0};

% approximate estimates
l0_approx={0,0};
l0_full_approx={0,0};
phi_het0_approx={0,0};
ls_approx = {zeros(1,size(plasmid_concs,2)), zeros(1,size(plasmid_concs,2))}; % growth rates
ls_full_approx = {zeros(1,size(plasmid_concs,2)), zeros(1,size(plasmid_concs,2))}; % growth rates with more accurate approximation
phi_hets_approx = {zeros(1,size(plasmid_concs,2)), zeros(1,size(plasmid_concs,2))}; % heterologous protein fractions

%% RUN for 2 conditions
for cond=1:1
    % set desired conditions
    if(cond==1)
        sim.parameters('nutr_qual')=0.5;
    else
        sim.parameters('nutr_qual')=0.2;
    end

    % values w/out het. prot.
    sim.het.parameters('c_xtra')=0;
    ss=get_steady(sim,Delta,Max_iter); % integrate
    [l,e,phi_het]=get_lephihet(sim,ss); % get desired values
    % record
    ss0{cond} = ss;
    l0{cond}=l;
    e0{cond}=e;
    phi_het0{cond}=phi_het;

    % get approximate estimates
    l0_approx{cond}=approx.ss_l(0,ss0{cond},e0{cond},sim);
    l0_full_approx{cond}=approx.ss_l_full(0,ss0{cond},e0{cond},sim);
    ls_approx{cond}=approx.ss_l(plasmid_concs,ss0{cond},e0{cond},sim);
    ls_full_approx{cond}=approx.ss_l_full(plasmid_concs,ss0{cond},e0{cond},sim);
    phi_hets_approx{cond}=approx.ss_phi_het(plasmid_concs,ss0{cond},e0{cond},sim);

    % get actual model predictions
    for i=1:size(plasmid_concs,2)
        % disp(plasmid_concs(i))
        sim.het.parameters('c_xtra')=plasmid_concs(i);
        ss=get_steady(sim,Delta,Max_iter); % integrate
        [l,e,phi_het]=get_lephihet(sim,ss); % get desired values
    
        % record
        sss{cond}(:,i) = ss;
        ls{cond}(i)=l;
        es{cond}(i)=e;
        phi_hets{cond}(i)=phi_het;
    end
end

%% PLOT and COMPARE with experimental data
colours={'r'};
colours_approx={'#A2142F'}; % colours for plotting approximations
Fig = figure('Position',[0 0 540 480]);
hold on
for cond=1:size(colours,2)
    plot([0,phi_hets{cond}],[l0{cond},ls{cond}]/l0{cond},[colours{cond},'-']) % plot model predictions
    plot([0,phi_hets_approx{cond}],[l0_approx{cond},ls_approx{cond}]/l0_approx{cond},...
        'Color',colours_approx{cond},'LineStyle','--') % plot approximate estimates
    plot([0,phi_hets_approx{cond}],[l0_full_approx{cond},ls_full_approx{cond}]/l0_full_approx{cond},...
        'Color',colours_approx{cond},'LineStyle',':') % plot approximate estimates
end
plot(ydata_xtra(:,1),ydata_xtra(:,2),'b.') % plot experimental data
xlabel('Het prot mass fraction');
ylabel('l with xtra:l without');
title('Growth dependence on het. prot. exp. vs 55 years of studies');
xlim([0 0.45])
ylim([0 1.05])
legend('Model pred.','Approximation', 'Better approx')
hold off

%% DECIDE whether to plot absolutes or ratios
absolutes = true;

%% PLOT growh rate ratios vs mass fractions, including experimental data
Fig1 = figure('Position',[0 0 540 480]);
hold on
for cond=1:size(colours,2)
    if(absolutes)
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[l0{cond},ls{cond}],[colours{cond},'-']) % plot model predictions
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[l0_approx{cond},ls_approx{cond}],...
            'Color',colours_approx{cond},'LineStyle','--') % plot approximation
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[l0_full_approx{cond},ls_full_approx{cond}],...
            'Color',colours_approx{cond},'LineStyle',':') % plot approximation
    else
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[l0{cond},ls{cond}]/l0{cond},[colours{cond},'-']) % plot model predictions
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[l0_approx{cond},ls_approx{cond}]/l0_approx{cond},...
            'Color',colours_approx{cond},'LineStyle','--') % plot approximation
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[l0_full_approx{cond},ls_full_approx{cond}]/l0_full_approx{cond},...
            'Color',colours_approx{cond},'LineStyle',':') % plot approximation
    end
end
%plot(ydata_xtra(:,1),ydata_xtra(:,2),'b.') % plot experimental data
xlabel('Het prot transcription rate (nM/h)');
if(absolutes)
    ylabel('Growth rate, /h');
else
    ylabel('l with xtra:l without');
    ylim([0 1])
end
title('Growth dependence on het. prot. exp.');
%xlim([0 0.2])
%ylim([0 1.05])
%legend('nutr. qual.=1','nutr. qual.=0.2')
hold off

%% PLOT cellular characteristics vs mass fractions
Fig2 = figure('Position',[0 0 1200 600]);

% Charged tRNA levels
subplot(2,3,1)
hold on
for cond=1:size(colours,2)
    if(absolutes)
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[ss0{cond}(5),sss{cond}(5,:)],[colours{cond},'-']) % plot model predictions
    else
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[ss0{cond}(5),sss{cond}(5,:)]/ss0{cond}(5),[colours{cond},'-'])
    end
end
xlabel('Het prot transcription rate (nM/h)');
if(absolutes)
    ylabel('t^c, nM');
else
    ylabel('t^c with xtra:t^c without');
    ylim([0 2.5])
end
title('t^c vs het. prot. exp.');
%xlim([0 0.2])
hold off

% Unharged tRNA levels
subplot(2,3,2)
hold on
for cond=1:size(colours,2)
    if(absolutes)
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[ss0{cond}(6),sss{cond}(6,:)],[colours{cond},'-'])
    else
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[ss0{cond}(6),sss{cond}(6,:)]/ss0{cond}(6),[colours{cond},'-'])
    end
end
xlabel('Het prot transcription rate (nM/h)');
if(absolutes)
    ylabel('t^u, nM');
else
    ylabel('t^u with xtra:t^u without');
    ylim([0 2.5])
end
title('t^u vs het. prot. exp.');
%xlim([0 0.2])
hold off

% T
subplot(2,3,3)
hold on
T0={};
Ts={};
for cond=1:size(colours,2)
    T0{cond}=ss0{cond}(5)/ss0{cond}(6);
    Ts{cond}=sss{cond}(5,:)./sss{cond}(6,:);
    if(absolutes)
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[T0{cond},Ts{cond}],[colours{cond},'-']) % plot model predictions
    else
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[T0{cond},Ts{cond}]/T0{cond},[colours{cond},'-'])
    end
end
xlabel('Het prot transcription rate (nM/h)');
if(absolutes)
    ylabel('T');
else
    ylabel('T with xtra:T without');
    ylim([0 1])
end
title('1/ppGpp levels vs het. prot. exp.');
%xlim([0 0.2])
hold off

% Elongation rates
subplot(2,3,4)
hold on
for cond=1:size(colours,2)
    if(absolutes)
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[e0{cond},es{cond}],[colours{cond},'-']) % plot model predictions
    else
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[e0{cond},es{cond}]/e0{cond},[colours{cond},'-']) % plot model predictions
    end
end
xlabel('Het prot transcription rate (nM/h)');
if(absolutes)
    ylabel('e, /h');
else
    ylabel('e with xtra:e without');
    ylim([0 1.25])
end
title('Translation rate vs het. prot. exp.');
%xlim([0 0.2])
hold off

% Ribosome/tRNA transcription
subplot(2,3,5)
hold on
F_r0={};
F_rs={};
for cond=1:size(colours,2)
    F_r0{cond}=sim.coll.F_r(sim.parameters,T0{cond});
    F_rs{cond}=sim.coll.F_r(sim.parameters,Ts{cond});
    if(absolutes)
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[F_r0{cond},F_rs{cond}],[colours{cond},'-']) % plot model predictions
    else
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[F_r0{cond},F_rs{cond}]/F_r0{cond},[colours{cond},'-'])
    end
end
xlabel('Het prot transcription rate (nM/h)');
if(absolutes)
    ylabel('F_r');
else
    ylabel('F_r with xtra:F_r without');
    ylim([0 1])
end
title('ppGpp signal vs het. prot. exp.');
%xlim([0 0.2])
hold off

% tRNA charging
subplot(2,3,6)
hold on
nu0={};
nus={};
for cond=1:size(colours,2)
    nu0{cond}=sim.coll.nu(sim.parameters,ss0{cond}(6));
    nus{cond}=sim.coll.nu(sim.parameters,sss{cond}(6,:));
    if(absolutes)
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[nu0{cond},nus{cond}],[colours{cond},'-']) % plot model predictions
    else
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[nu0{cond},nus{cond}]/nu0{cond},[colours{cond},'-'])
    end
end
xlabel('Het prot transcription rate (nM/h)');
if(absolutes)
    ylabel('phi_r');
else
    ylabel('nu with xtra:nu without');
    ylim([0 1.05])
end
title('tRNA charging per p_a vs het. prot. exp.');
%xlim([0 0.2])
hold off

%% PLOT R:p_a
Fig3 = figure('Position',[0 0 540 480]);
hold on
Rpa0={};
Rpas={};
% approximations
Pras0_approx={};
Rpas_approx={};
for cond=1:size(colours,2)
    % actual predictions
    Rpa0{cond}=ss0{cond}(4)/ss0{cond}(3);
    Rpas{cond}=sss{cond}(4,:)./sss{cond}(3,:);
    % approximation
    Rpa0_approx{cond}=(approx.ss_phi(0,F_r0{cond},sim.parameters('c_r'),sim.parameters('a_r'),...
        sim.parameters('k+_r'),sim.parameters('k-_r'),sim.parameters('n_r'),ss0{cond},e0{cond},sim)./sim.parameters('n_r'))./...
        (approx.ss_phi(0,1,sim.parameters('c_a'),sim.parameters('a_a'),...
        sim.parameters('k+_a'),sim.parameters('k-_a'),sim.parameters('n_a'),ss0{cond},e0{cond},sim)./sim.parameters('n_a'));
    Rpas_approx{cond}=(approx.ss_phi(plasmid_concs,F_r0{cond},sim.parameters('c_r'),sim.parameters('a_r'),...
        sim.parameters('k+_r'),sim.parameters('k-_r'),sim.parameters('n_r'),ss0{cond},e0{cond},sim)./sim.parameters('n_r'))./...
        (approx.ss_phi(plasmid_concs,1,sim.parameters('c_a'),sim.parameters('a_a'),...
        sim.parameters('k+_a'),sim.parameters('k-_a'),sim.parameters('n_a'),ss0{cond},e0{cond},sim)./sim.parameters('n_a')); % NOTE: F_r assumed UNCHANGING!
    if(absolutes)
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[Rpa0{cond},Rpas{cond}],[colours{cond},'-']) % plot model predictions
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[Rpa0_approx{cond},Rpas_approx{cond}],...
            'Color',colours_approx{cond},'LineStyle','--') % plot approximation
    else
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[Rpa0{cond},Rpas{cond}]/Rpa0{cond},[colours{cond},'-'])
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[Rpa0_approx{cond},Rpas_approx{cond}]/Rpa0_approx{cond},...
            'Color',colours_approx{cond},'LineStyle','--') % plot approximation
    end
end
xlabel('Het prot transcription rate (nM/h)');
if(absolutes)
    ylabel('R/p_a');
else
    ylabel('R/p_a with xtra:R/p_a without');
    ylim([0 1.05])
end
ylabel('R:p_a w/:w/out het. prot. ratio');
title('Ribosomes vs met. prot.');
%xlim([0 0.45])
%legend('nutr. qual.=1','nutr. qual.=0.2')
hold off

%% PLOT mRNA levels
Fig4 = figure('Position',[0 0 1500 400]);

% Heterologous
if(absolutes)
    subplot(1,3,1)
    hold on
    for cond=1:size(colours,2)
        % find approximation
        m_het0_approx=approx.ss_m(0,1,0,sim.het.parameters('a_xtra'),sim.het.parameters('b_xtra'),...
            ss0{cond},e0{cond},sim);
        m_hets_approx=approx.ss_m(plasmid_concs,1,plasmid_concs,sim.het.parameters('a_xtra'),sim.het.parameters('b_xtra'),...
            ss0{cond},e0{cond},sim);
        % plot
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[0,sss{cond}(8,:)],[colours{cond},'-']) % plot model predictions
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[m_het0_approx,m_hets_approx],...
            'Color',colours_approx{cond},'LineStyle','--')
    end
    xlabel('Het prot transcription rate (nM/h)');
    ylabel('m_{xtra}, nM');
    title('Het. mRNA levels');
    hold off
end

% Ribosomal
subplot(1,3,2)
hold on
for cond=1:size(colours,2)
    % get approximations
    m_r0_approx=approx.ss_m(0,F_r0{cond},sim.parameters('c_r'),sim.parameters('a_r'),sim.parameters('b_r'),...
        ss0{cond},e0{cond},sim);
    m_rs_approx=approx.ss_m(plasmid_concs,F_r0{cond},sim.parameters('c_r'),sim.parameters('a_r'),sim.parameters('b_r'),...
        ss0{cond},e0{cond},sim); % NOTE: F_r assumed UNCHANGING!
    if(absolutes)
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[ss0{cond}(2),sss{cond}(2,:)],[colours{cond},'-']) % plot model predictions
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[m_r0_approx,m_rs_approx],...
            'Color',colours_approx{cond},'LineStyle','--')
    else
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[ss0{cond}(2),sss{cond}(2,:)]/ss0{cond}(2),[colours{cond},'-']) % plot model predictions
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[m_r0_approx,m_rs_approx]/m_r0_approx,...
            'Color',colours_approx{cond},'LineStyle','--')
    end
end
xlabel('Het prot transcription rate (nM/h)');
if(absolutes)
    ylabel('m_r, nM');
else
    ylabel('m_r with xtra:m_r without');
    ylim([0 1.25])
end
title('Ribosomal mRNA levels');
hold off

% Metabolic
subplot(1,3,3)
hold on
for cond=1:size(colours,2)
    % get approximations
    m_a0_approx=approx.ss_m(0,1,sim.parameters('c_a'),sim.parameters('a_a'),sim.parameters('b_a'),ss0{cond},e0{cond},sim);
    m_as_approx=approx.ss_m(plasmid_concs,1,sim.parameters('c_a'),sim.parameters('a_a'),sim.parameters('b_a'),ss0{cond},e0{cond},sim);
    if(absolutes)
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[ss0{cond}(1),sss{cond}(1,:)],[colours{cond},'-']) % plot model predictions
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[m_a0_approx,m_as_approx],...
            'Color',colours_approx{cond},'LineStyle','--')
    else
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[ss0{cond}(1),sss{cond}(1,:)]/ss0{cond}(1),[colours{cond},'-'])
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[m_a0_approx,m_as_approx]/m_a0_approx,...
            'Color',colours_approx{cond},'LineStyle','--')
    end
end
xlabel('Het prot transcription rate (nM/h)');
if(absolutes)
    ylabel('m_a, nM');
else
    ylabel('m_a with xtra:m_a without');
    ylim([0 1.2])
end
title('Metabolic mRNA levels');
hold off

%% PLOT actual mass fractions
Fig5 = figure('Position',[0 0 1500 400]);

% Heterologous
if(absolutes)
    subplot(1,3,1)
    hold on
    for cond=1:size(colours,2)
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[phi_het0{cond},phi_hets{cond}],[colours{cond},'-']) % plot model predictions
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[phi_het0_approx{cond},phi_hets_approx{cond}], ...
            'Color',colours_approx{cond},'LineStyle','--') % plot approximation
    end
    xlabel('Het prot transcription rate (nM/h)');
    ylabel('phi_{xtra}');
    title('Het. prot. mass frac');
    hold off
end

% Ribosomal
subplot(1,3,2)
hold on
phi_rs={};
phi_r0={};
% approximates
for cond=1:size(colours,2)
    phi_r0{cond}=ss0{cond}(4).*sim.parameters('n_r')./sim.parameters('M');
    phi_rs{cond}=sss{cond}(4,:).*sim.parameters('n_r')./sim.parameters('M');
    % get approximates
    phi_r0_approx{cond}=approx.ss_phi(0,F_r0{cond},sim.parameters('c_r'),sim.parameters('a_r'),...
        sim.parameters('k+_r'),sim.parameters('k-_r'),sim.parameters('n_r'),ss0{cond},e0{cond},sim);
    phi_rs_approx{cond}=approx.ss_phi(plasmid_concs,F_r0{cond},sim.parameters('c_r'),sim.parameters('a_r'),...
        sim.parameters('k+_r'),sim.parameters('k-_r'),sim.parameters('n_r'),ss0{cond},e0{cond},sim); % NOTE: F_r assumed UNCHANGING!
    if(absolutes)
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[phi_r0{cond},phi_rs{cond}],[colours{cond},'-']) % plot model predictions
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[phi_r0_approx{cond},phi_rs_approx{cond}], ...
            'Color',colours_approx{cond},'LineStyle','--') % plot approximation
    else
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[phi_r0{cond},phi_rs{cond}]/phi_r0{cond},[colours{cond},'-']) % plot model predictions
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[phi_r0_approx{cond},phi_rs_approx{cond}]/phi_r0_approx{cond}, ...
            'Color',colours_approx{cond},'LineStyle','--') % plot approximation
    end
end
xlabel('Het prot transcription rate (nM/h)');
if(absolutes)
    ylabel('phi_r');
else
    ylabel('phi_r with xtra:phi_r without');
    ylim([0 1])
end
title('Ribosome mass frac');
hold off

% Metabolic
subplot(1,3,3)
hold on
phi_as={};
phi_a0={};
for cond=1:size(colours,2)
    phi_a0{cond}=ss0{cond}(3).*sim.parameters('n_a')./sim.parameters('M');
    phi_as{cond}=sss{cond}(3,:).*sim.parameters('n_a')./sim.parameters('M');
    % get approximations
    phi_a0_approx{cond}=approx.ss_phi(0,1,sim.parameters('c_a'),sim.parameters('a_a'),...
        sim.parameters('k+_a'),sim.parameters('k-_a'),sim.parameters('n_a'),ss0{cond},e0{cond},sim);
    phi_as_approx{cond}=approx.ss_phi(plasmid_concs,1,sim.parameters('c_a'),sim.parameters('a_a'),...
        sim.parameters('k+_a'),sim.parameters('k-_a'),sim.parameters('n_a'),ss0{cond},e0{cond},sim);
    if(absolutes)
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[phi_a0{cond},phi_as{cond}],[colours{cond},'-']) % plot model predictions
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[phi_a0_approx{cond},phi_as_approx{cond}], ...
            'Color',colours_approx{cond},'LineStyle','--') % plot approximation
    else
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[phi_a0{cond},phi_as{cond}]/phi_a0{cond},[colours{cond},'-'])
        plot([0,plasmid_concs*sim.het.parameters('a_xtra')],[phi_a0_approx{cond},phi_as_approx{cond}]/phi_a0_approx{cond}, ...
            'Color',colours_approx{cond},'LineStyle','--') % plot approximation
    end
end
xlabel('Het prot transcription rate (nM/h)');
if(absolutes)
    ylabel('phi_a');
else
    ylabel('phi_a with xtra:phi_a without');
    ylim([0 1])
end
title('Met. prot. mass frac');
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

    % USEFUL PRE-CALCULATIONS
    % translation elongation rate
    e=sim.coll.e(par,tc);

    % ribosome inactivation rate due to chloramphenicol
    kcmh=par('kcm').*par('h');

    % ribosome dissociation constants
    k_a=sim.coll.k(e,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
    k_r=sim.coll.k(e,par('k+_r'),par('k-_r'),par('n_r'),kcmh);

    % consider the contribution of heterologous genes to resource competition if there are any
    mk_het=0;
    if(sim.het.num_genes>0)
        ss_het=ss(8:end); % current state of heterologous genes
        sim.het = sim.het.find_current_ks(e,kcmh); % ribosome dissociation constants for heterolgous genes
        for i=1:sim.het.num_genes
            mk_het = mk_het + ss_het(i)./sim.het.current_ks(i); % get sum of mx/kx for heterologous genes
        end
    end

    T=tc./tu; % ratio of charged to uncharged tRNAs
    D=1+(m_a./k_a+m_r./k_r+mk_het)./(1-par('phi_q')); % denominator in ribosome competition cobj.optalculations
    B=R.*(1-1./D); % actively translating ribosomes - INCLUDING Q

    % tRNA charging and synthesis
    nu=sim.coll.nu(par,tu);
    psi=sim.coll.psi(par,T);

    % growth rate
    l=sim.coll.l(par,e,B);

    % heterologous protein mass fraction
    phi_het=0;
    for i=1:sim.het.num_genes
        phi_het=phi_het+ss_het(sim.het.num_genes+i).*sim.het.parameters(['n_',sim.het.names{i}])./sim.parameters('M');
    end
end