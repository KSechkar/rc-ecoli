%% AppFigure1.m
% Generate Appendix Figure 1 - example simulation of the system

%% CLEAR parameters, add paths of all files

addpath(genpath('..'))

clear
close all 

%% SIMULATE the model
sim = cell_simulator;
sim.tf=4;
sim.init_conditions('s')=0.5;

sim = sim.simulate_model;

%% Open Figure
Faa = figure('Position',[0 0 632 560]);
set(Faa, 'defaultLineLineWidth', 1.25)
set(Faa, 'defaultAxesFontSize', 9)

%% Figure A1a - cell composition
% Calculate masses of different proteins in amino acids
amino_q=sim.parameters('M')*ones(size(sim.t))*sim.parameters('phi_q');
amino_a=sim.x(:,3)*sim.parameters('n_a');
amino_r=sim.x(:,4)*sim.parameters('n_r');
amino_Bcm=sim.x(:,7)*sim.parameters('n_r');
amino_het=zeros(size(sim.t)); % initialise masses of heterologous proteins
if(sim.num_het~=0)
    x_het=sim.x(:,10: (9+2*sim.num_het) );
    for i=1:sim.num_het
        amino_het=sum((amino_het+x_het(:,(sim.num_het+i))...
            *sim.parameters(['n_',sim.het.names{i}])),2);
    end
end

%PLOT
subplot(2,2,1)
hold all
patch([sim.t; flip(sim.t)],[amino_q+amino_a+amino_r+amino_Bcm+amino_het; flip(amino_a+amino_r+amino_Bcm+amino_het)], [0.75, 0.75, 0.75])
patch([sim.t; flip(sim.t)],[amino_a+amino_r+amino_Bcm+amino_het; flip(amino_r+amino_Bcm)], [0.9290, 0.6940, 0.1250])
patch([sim.t; flip(sim.t)],[amino_r+amino_Bcm+amino_het; flip(amino_Bcm+amino_het)], [0.4940, 0.1840, 0.5560])
patch([sim.t; flip(sim.t)],[amino_Bcm+amino_het; flip(amino_het)], [0.4660 0.6740 0.1880])
patch([sim.t; flip(sim.t)],[amino_het; zeros(size(sim.t))], [0 0.4470 0.7410])

xlabel('t, time [h]');
ylabel('n_i p_i, protein mass [aa]');
legend('housekeeping','metabolic','ribosomes');
hold off

%% Figure A1b - mRNA concentrations
subplot(2,2,2)
hold all
plot(sim.t,sim.x(:,1),'Color',[0.9290, 0.6940, 0.1250]);
plot(sim.t,sim.x(:,2),'Color',[0.4940, 0.1840, 0.5560]);

grid on

xlabel('t, time [h]');
ylabel('m_i, mRNA concentration [nM]');
legend('metabolic','ribosomal','Location','east');
hold off

%% Figure A1c - protein concentrations
subplot(2,2,3)
hold all
plot(sim.t,sim.x(:,3),'Color',[0.9290, 0.6940, 0.1250]);
plot(sim.t,sim.x(:,4),'Color',[0.4940, 0.1840, 0.5560]);

grid on

xlabel('t, time [h]');
ylabel('p_i, protein concentration [nM]');
legend('metabolic','ribosomal','Location','east');
hold off

%% Figure A1d - tRNA concentrations
subplot(2,2,4)
hold all
plot(sim.t,sim.x(:,6),sim.t,sim.x(:,5))

grid on

xlabel('t, time [h]');
ylabel('tRNA concentration [nM]');
legend('charged','uncharged');
hold off