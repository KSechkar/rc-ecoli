%% Figure3d.m
% Plot growth rates for different fixed values of the ribosomal gene
% transcription function, compare with the growth produced by Flux-Parity
% Regulation

%% CLEAR parameters, add paths of all files

addpath(genpath('..'))

clear
close all

%% SET UP the simulator

sim=cell_simulator; % initialise simulator

% parameters for getting steady state
sim.tf = 10; % single integraton step timeframe
Delta = 0.01; % threshold that determines if we're in steady state
Max_iter = 75; % maximum no. iterations (checking if SS reached over first 750 h)

sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9); % more lenient integration tolerances for speed

%% DEFINE nutrient condiitons to explore

% vector of nurtient qualities
nutrients=flip([0.05,0.1,0.25,0.5,0.75,1]);

% range of fixed ribosome transcription regulation function values
fixed_F_rs= 0.01:0.01:1;

%% INITIALISE arrays that will store differen growth rates
% storing growth rates - optimal and flux-parity
l_map = zeros(2,size(nutrients,2));

% storing ribosome transcription regulation function values - optimal and flux-parity
F_r_map = zeros(2,size(nutrients,2));

% storing growth rates - for different fixed regulation function values
ls=zeros(size(nutrients,2),size(fixed_F_rs,2));

%% RUN simulations
for j=1:size(nutrients,2)
    % RESET
    sim=sim.set_default_parameters(); % reset initial consitions and parameters
    sim=sim.set_default_init_conditions(); % reset initial consitions and parameters
    sim.init_conditions('s')=nutrients(j); % set nutrient quality

    % Run with regulated ribosome transcription
    sim.parameters('is_fixed_F_r')=0; % F_r regulated now

%     [l, F_r, ~] = get_steady_for_sweep(sim,Delta,Max_iter);
    ss=get_steady(sim,Delta,Max_iter);
    [l,F_r]=get_lFr(sim,ss);
    l_map(1,j)=l;
    F_r_map(1,j)=F_r;
    
    % RESET
    sim=sim.set_default_parameters(); % reset initial consitions and parameters
    sim=sim.set_default_init_conditions(); % reset initial consitions and parameters
    sim.init_conditions('s')=nutrients(j); % set nutrient quality    

    % Run with a range of fixed F_r values
    for i=1:size(fixed_F_rs,2)
        sim=sim.set_default_parameters(); % reset initial consitions and parameters
        sim=sim.set_default_init_conditions(); % reset initial consitions and parameters
        sim.init_conditions('s')=nutrients(j); % set nutrient quality
        sim.parameters('is_fixed_F_r')=1; % F_r fixed now
        
        disp(fixed_F_rs(i))
        sim.parameters('fixed_F_r')=fixed_F_rs(i);
%         [l, F_r, ~] = get_steady_for_sweep(sim,Delta,Max_iter);
        ss=get_steady(sim,Delta,Max_iter);
        [l,F_r]=get_lFr(sim,ss);
        ls(j,i)=l;
    end

    % finding optimal ribosome content
    [maxl,maxlpos]=max(ls(j,:));
    l_map(2,j)=maxl;
    F_r_map(2,j)=fixed_F_rs(maxlpos);
    disp(j)
end

%% Plot
Fig = figure('Position',[0 0 316 240]);
set(Fig, 'defaultAxesFontSize', 9)
set(Fig, 'defaultLineLineWidth', 1)
hold on
for j=1:size(nutrients,2)
    linecolour=0.75*(1-nutrients(j));
    plot(fixed_F_rs,ls(j,:),'Color',[linecolour, linecolour,linecolour]) % different fixed F_r values
    
    plot(F_r_map(1,j),l_map(1,j),'o','Color','r') % flux-parity regulation
    plot(F_r_map(2,j),l_map(2,j),'o','Color',[0,0.75,0]) % optimal regulation
end

xlabel({'F_r, ribosomal gene transcription','regulation func.'});
ylabel('Growth rate \lambda, 1/h');
ylim([0 2.5])
grid on
hold off

%% FUNCTION for getting growth rate and F_r from the system's steady state
function [l,F_r]=get_lFr(sim,ss)
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

    e=sim.form.e(par,tc); % translation elongation ratio
    k_a=sim.form.k(e,par('k+_a'),par('k-_a'),par('n_a'),par('kcm').*h); % ribosome-mRNA dissociation constant (metab. genes)
    k_r=sim.form.k(e,par('k+_r'),par('k-_r'),par('n_r'),par('kcm').*h); % ribosome-mRNA dissociation constant (rib. genes)
    D=1+(m_a./k_a+m_r./k_r)./(1-par('phi_q')); % denominator in ribosome competition calculations
    B=R.*(1-1./D); % actively translating ribosomes (inc. those translating housekeeping genes)

    l=sim.form.l(par,e,B); % growth rate!
    
    T = tc./tu; % ratio of charged to uncharged tRNAs
    F_r = sim.form.F_r(par,T) ; % ribosome regulation
end