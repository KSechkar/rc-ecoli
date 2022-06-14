%% Figure9ab.m
% Generate Figure 9a,b - parameters and reduced leakiness

addpath(genpath('..'))

clear
close all 

%% SET UP and RUN the simulators for both cases

% define how many resacalings are to be performed
rescales=logspace(-2,2,7);

for cond=1:size(rescales,2)
    sim{cond}=advanced_simulator;

    % disturbance signal parameters
    sim{cond}=sim{cond}.load_heterologous_and_external('integral_dirdist','pulse_inducer');
    sim{cond}.ext.input_func_parameters('inducer_base_level')=1;
    sim{cond}.ext.input_func_parameters('pulse_value_prop')=0;
    sim{cond}.ext.input_func_parameters('pulse_start_time')=0;
    sim{cond}.ext.input_func_parameters('pulse_duration')=21;
    sim{cond}.het.parameters('a_dist')=300;
    
    % integral controller parameters
    sim{cond}.het.parameters('K_dna-sens')=4000;
    sim{cond}.het.parameters('eta_dna-sens')=1;
    sim{cond}.het.parameters('kb_anti')=300;
    sim{cond}.het.parameters('a_sens')=15;
    sim{cond}.het.parameters('a_anti')=2000;
    sim{cond}.het.parameters('a_act')=1000;
    
    % to mitigate leakiness, scale the synthesis and binding rates of aniilator and actuator
    sim{cond}.het.parameters('a_anti')=sim{cond}.het.parameters('a_anti').*rescales(cond);
    sim{cond}.het.parameters('a_act')=sim{cond}.het.parameters('a_act').*rescales(cond);
    sim{cond}.het.parameters('kb_anti')=sim{cond}.het.parameters('kb_anti').*rescales(cond);
       
    % push amended parameter values
    sim{cond}=sim{cond}.push_het();

    % simulate
    sim{cond}.tf = 30;
    sim{cond}.opt = odeset('reltol',1.e-6,'abstol',1.e-9);
    sim{cond} = sim{cond}.simulate_model;

end


%% GET relevant time frame (start 1h before disturbance, end 4h after)
for cond=1:size(rescales,2)
    % find first point in the frame
    for i=1:size(sim{cond}.t,1)
        if(sim{cond}.t(i)>=16)
            first_pt=i;
            break
        end
    end
    
    % record relevant time points
    rel_t{cond}=sim{cond}.t(first_pt-1:end);
    
    % record states of the cell at these times
    rel_x{cond}=sim{cond}.x(first_pt-1:end,:);
end

%% CALCULATE growth rate and annihilator transcription regulation func.

for cond=1:size(rescales,2)
    ls{cond}=zeros(size(rel_t{cond}));
    F_antis{cond}=zeros(size(rel_t{cond}));   
    
    par=sim{cond}.parameters;
    for i=1:size(rel_t{cond})
        % STATE VECTOR TO SINGLE VARIABLES
        m_a = rel_x{cond}(i,1);
        m_r = rel_x{cond}(i,2);
        p_a = rel_x{cond}(i,3);
        R = rel_x{cond}(i,4);
        tc = rel_x{cond}(i,5);
        tu = rel_x{cond}(i,6);
        Bcm = rel_x{cond}(i,7);
        s = rel_x{cond}(i,8);
        h = rel_x{cond}(i,9);
        x_het=rel_x{cond}(i,10 : (9+2*sim{cond}.num_het) );
    
        % USEFUL PRE-CALCULATIONS
        % translation elongation rate
        e=sim{cond}.coll.e(par,tc);
    
        % ribosome inactivation rate due to chloramphenicol
        kcmh=par('kcm').*h;
    
        % ribosome dissociation constants
        k_a=sim{cond}.coll.k(e,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
        k_r=sim{cond}.coll.k(e,par('k+_r'),par('k-_r'),par('n_r'),kcmh);
        % heterologous genes
        k_het=ones(1,sim{cond}.num_het);
        if(sim{cond}.num_het>0)
            for j=1:sim{cond}.num_het
                k_het(j)=sim{cond}.coll.k(e,...
                sim{cond}.parameters(['k+_',sim{cond}.het.names{j}]),...
                sim{cond}.parameters(['k-_',sim{cond}.het.names{j}]),...
                sim{cond}.parameters(['n_',sim{cond}.het.names{j}]),...
                kcmh);
            end
        end
    
        T=tc./tu; % ratio of charged to uncharged tRNAs
        D=1+(m_a./k_a+m_r./k_r+sum(x_het(1:sim{cond}.num_het)./k_het))./...
            (1-par('phi_q')); % denominator in ribosome competition csim.optalculations
        B=R.*(1-1./D); % actively translating ribosomes - INCLUDING Q
    
        % growth rate
        l=sim{cond}.coll.l(par,e,B);
    
        % RECORD VALUES
        ls{cond}(i)=l;
        F_antis{cond}(i)=sim{cond}.het.regulation(sim{cond}.het.names{2},rel_x{cond}(i,:),0); % external input can be 0 as m_anti does not depend on it
    end
end

%% Figure 9a - anihilator transcription regulation

Fig9a = figure('Position',[0 0 316 280]);
set(Fig9a, 'defaultAxesFontSize', 9)
set(Fig9a, 'defaultLineLineWidth', 1.25)

hold on

% plot predictions for different cases
for cond=1:size(rescales,2)
    cond_colour=[0, 0.4470, 0.7410, 1-0.6*cond/size(rescales,2)];
    plot(rel_t{cond}-21,F_antis{cond},'Color',cond_colour)
end
plot([-5 5],[0.5 0.5],'k:') % plot ideal value

xlim([-5 5])
ylim([0.45 0.51])

xticks(-5:2.5:5)
grid on
grid minor

xlabel('t, time since disturbance [h]')
ylabel({'F_{anti}, reg. of annihlator transcription'})

hold off

%% Figure 9c - growth rate
Fig9b = figure('Position',[0 0 316 280]);
set(Fig9b, 'defaultAxesFontSize', 9)
set(Fig9b, 'defaultLineLineWidth', 1.25)

hold on

% plot predictions for different cases
for cond=1:size(rescales,2)
    cond_colour=[0, 0.4470, 0.7410, 1-0.6*cond/size(rescales,2)];
    plot(rel_t{cond}-21,ls{cond},'Color',cond_colour)
end
plot([-5 5],[1.252 1.252],'k:') % plot ideal value

xlim([-5 5])
ylim([1.2 1.6])

xticks(-5:2.5:5)
grid on
grid minor

xlabel('t, time since disturbance [h]')
ylabel({'\lambda, growth rate [1/h]'})

hold off