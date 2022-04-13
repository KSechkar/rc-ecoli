%% main_code.m
%%% OCTOBER 11, 2021

clear
close all 

%% Create the object 's' from the class 'simulator'
sim = advanced_simulator;
sim.parameters('using_nutr_qual') = 1;
sim.parameters('nutr_qual') = 0.5;
sim.het.parameters('c_xtra')=10;
sim.tf = 20;

%% Execute the simulation
sim = sim.simulate_model;

%% Plot the simulation
sim.plot_simulation;

%% Plot other things - concentrations
Fm = figure('Position',[0 0 1000 800]);
set(Fm, 'defaultLineLineWidth', 2)
set(Fm, 'defaultAxesFontSize', 14)

% mRNA levels
subplot(2,2,1)
hold on
plot(sim.t,sim.x(:,1),'Color',[0.9290, 0.6940, 0.1250]);
plot(sim.t,sim.x(:,2),'Color',[0.4940, 0.1840, 0.5560]);
xlabel('Time (h)');
ylabel('mRNA concentration (nM)');
legend('a','r');
title('mRNA levels');
hold off

% protein levels
subplot(2,2,2)
hold on
plot(sim.t,sim.x(:,3),'Color',[0.9290, 0.6940, 0.1250]);
plot(sim.t,sim.x(:,4),'Color',[0.4940, 0.1840, 0.5560]);
xlabel('Time (h)');
ylabel('mRNA concentration (nM)');
legend('a','r');
title('protein levels');
hold off

% tRNA levels
subplot(2,2,3)
plot(sim.t,sim.x(:,5),sim.t,sim.x(:,6))
xlabel('Time (h)');
ylabel('tRNA concentration (nM)');
legend('charged','uncharged');
title('tRNA charging');

%%
% Rates and other variables
es=zeros(size(sim.t));
nus=zeros(size(sim.t));
phis=zeros(size(sim.t));
ls=zeros(size(sim.t));
F_rs=zeros(size(sim.t));
ks=zeros(size(sim.t));
Bs=zeros(size(sim.t));
Ts=zeros(size(sim.t));
Ds=zeros(size(sim.t));

whatnot=zeros(size(sim.t));

par=sim.parameters;
for i=1:size(sim.t)
    m_a = sim.x(i,1);
    m_r = sim.x(i,2);
    p_a = sim.x(i,3);
    R = sim.x(i,4);
    tc = sim.x(i,5);
    tu = sim.x(i,6);

    e=sim.coll.e(par,tc);
    % ribosome dissociation constants
    kcmh=par('kcm').*par('h');
    k_a=sim.coll.k(e,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
    k_r=sim.coll.k(e,par('k+_r'),par('k-_r'),par('n_r'),kcmh);

    % commonly used expressions
    T=tc./tu; % ratio of charged to uncharged tRNAs
    D=1+(m_a./k_a+m_r./k_r)./(1-par('phi_q')); % denominator in ribosome competition calculations
    %B=b_r+b_a+b_q; % actively translating ribosomes % full only!
    B=R.*(1-1./D);

    % tRNA charging and synthesis
    nu=sim.coll.nu(par,tu);
    psi=sim.coll.psi(par,T);

    % growth rate
    l=sim.coll.l(par,e,B);

    es(i)=e;
    nus(i)=nu;
    phis(i)=psi;
    ls(i)=l;

    F_rs(i)=sim.coll.F_r(par,T);

    ks(i)=k_a;

    Bs(i)=B;
    Ts(i)=T;
    Ds(i)=D;
end

Fr = figure('Position',[0 0 1500 800]);
set(Fr, 'defaultLineLineWidth', 1)
set(Fr, 'defaultAxesFontSize', 12)

% elongation rate
subplot(2,3,1)
plot(sim.t,es)
xlabel('Time (h)');
ylabel('Elongation rate');
title('Elongation rate');
ylim([0 80000]);

% tRNA charging rate
subplot(2,3,2)
plot(sim.t,nus)
xlabel('Time (h)');
ylabel('tRNA charging rate');
title('tRNA charging');

% tRNA synthesis rate
subplot(2,3,3)
plot(sim.t,phis)
xlabel('Time (h)');
ylabel('tRNA synthesis rate');
title('tRNA synthesis');

% growth rate
subplot(2,3,4)
plot(sim.t,ls)
xlabel('Time (h)');
ylabel('l, 1/h');
title('Growth rate');
ylim([0 2]);

% regulation
subplot(2,3,5)
hold on
plot(sim.t,F_rs,'Color',[0.4940, 0.1840, 0.5560])
hold off
xlabel('Time (h)');
ylabel('Hill function values');
title('Regulation of rib. transc.');
legend('R')
ylim([0 1]);


% ribosome affinity
subplot(2,3,6)
plot(sim.t,ks)
xlabel('Time (h)');
ylabel('k_a');
title('ribosome affinity');

%% PRE-CALCULATED VALUES
Fr = figure('Position',[0 0 1000 800]);
set(Fr, 'defaultLineLineWidth', 1)
set(Fr, 'defaultAxesFontSize', 12)

% tRNA ratios
subplot(2,2,1)
plot(sim.t,log(Ts))
xlabel('Time (h)');
ylabel('ln(T)');
title('tRNA ratios');

% active ribosomes
subplot(2,2,2)
hold on
plot(sim.t,Bs)
plot(sim.t,sim.x(:,4),'Color',[0.4940, 0.1840, 0.5560])
hold off
xlabel('Time (h)');
ylabel('Concentration');
title('Active and overall ribosomes');
legend('B','R')

% common denominator
subplot(2,2,3)
hold on
plot(sim.t,Ds)
xlabel('Time (h)');
ylabel('D');
title('Competition denominator');





