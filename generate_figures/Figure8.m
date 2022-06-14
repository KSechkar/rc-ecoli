%% Figure8.m
% Generate Figure 8 - ideal vs leaky controller

addpath(genpath('..'))

clear
close all 

%% SET UP and RUN the simulators for both cases

% 2 different setups - ideal AIF and realistic scenario
sim={advanced_simulator_ideal_aif,advanced_simulator};

for cond=1:2
    % disturbance signal parameters
    sim{cond}=sim{cond}.load_heterologous_and_external('integral_dirdist','pulse_inducer');
    sim{cond}.ext.input_func_parameters('inducer_base_level')=1;
    sim{cond}.ext.input_func_parameters('pulse_value_prop')=0;
    sim{cond}.ext.input_func_parameters('pulse_start_time')=0;
    sim{cond}.ext.input_func_parameters('pulse_duration')=30;
    sim{cond}.het.parameters('a_dist')=300;
    
    % integral controller parameters
    sim{cond}.het.parameters('K_dna-sens')=4000;
    sim{cond}.het.parameters('eta_dna-sens')=1;
    sim{cond}.het.parameters('kb_anti')=300;
    sim{cond}.het.parameters('a_sens')=15;
    sim{cond}.het.parameters('a_anti')=2000;
    sim{cond}.het.parameters('a_act')=1000;
       
    % push amended parameter values
    sim{cond}=sim{cond}.push_het();

    % simulate
    sim{cond}.tf = 60;
    sim{cond}.opt = odeset('reltol',1.e-6,'abstol',1.e-9);
    sim{cond} = sim{cond}.simulate_model;

end


%% GET relevant time frame (start 1h before disturbance, end 4h after)
for cond=1:2
    % find first point in the frame
    for i=1:size(sim{cond}.t,1)
        if(sim{cond}.t(i)>=15)
            first_pt=i;
            break
        end
    end
    
    % record relevant time points
    rel_t{cond}=sim{cond}.t(first_pt-1:end);
    
    % record states of the cell at these times
    rel_x{cond}=sim{cond}.x(first_pt-1:end,:);
end

%% Figure 8a,8c: evolution of heterologous mRNA levels with leaky and improved controller

colours=[[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.3010 0.7450 0.9330]; [0.6350 0.0780 0.1840]; [1 0 1]];% colours for plots

for cond=1:2
    if(cond==1)
        Fig8a = figure('Position',[0 0 316 280]);
        set(Fig8a, 'defaultAxesFontSize', 9)
        set(Fig8a, 'defaultLineLineWidth', 1.25)
    else
        Fig8c = figure('Position',[0 0 316 280]);
        set(Fig8c, 'defaultAxesFontSize', 9)
        set(Fig8c, 'defaultLineLineWidth', 1.25)
    end
    
    hold on
    for i=1:sim{cond}.num_het
        plot(rel_t{cond}-30,rel_x{cond}(:,9+i),'Color',colours(rem(i-1,5)+1,:));
    end
    
    xlim([-15 15])
    
    xticks(-15:5:15)
    grid on
    grid minor
    
    xlabel('t, time since disturbance [h]')
    ylabel('m_i, mRNA concentration [nM]')
    legend('m_{sens}','m_{anti}','m_{act}','m_{dist}', 'Location','northwest')
    if(cond==1)
        ylim([-500 20000])
        legend('m_{sens}','m_{anti}','m_{act}','m_{dist}', 'Location','west')
    else
        ylim([-500 5000])
        legend('m_{sens}','m_{anti}','m_{act}','m_{dist}', 'Location','northwest')
    end
    hold off
end


%% CALCULATE growth rate and annihilator transcription regulation func.

for cond=1:2
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

%% Figure 8b,d - anihilator transcription regulation

for cond=1:2
    if(cond==1)
        Fig8b = figure('Position',[0 0 316 280]);
        set(Fig8b, 'defaultAxesFontSize', 9)
        set(Fig8b, 'defaultLineLineWidth', 1.25)
    else
        Fig8d = figure('Position',[0 0 316 280]);
        set(Fig8d, 'defaultAxesFontSize', 9)
        set(Fig8d, 'defaultLineLineWidth', 1.25)
    end
    hold on
    
    % plot predictions
    plot(rel_t{cond}-30,F_antis{cond},'Color',[0, 0.4470, 0.7410])

    plot([-15 15],[0.5 0.5],'k:') % plot ideal value
    
    ylim([0.45 0.51])
    xlim([-15 15])
    
    xticks(-15:5:15)
    grid on
    grid minor
    
    xlabel('t, time since disturbance [h]')
    ylabel({'F_{anti}, reg. of annihlator transcription'})
    
    hold off
end