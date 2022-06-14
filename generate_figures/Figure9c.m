%% CLEAR

addpath(genpath('..'))

clear
close all

%% SET UP the simulator
sim=advanced_simulator; % initialise simulator
sim=sim.load_heterologous_and_external('integral_dirdist','constant_inducer');

% error tolerances during integration
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9);

sim.tf=10;

% params for finding steady state
Delta=0.01;
Max_iter=75;

%% DEFINE perturbation to be tested

sim.het.parameters('a_dist')=300;
sim=sim.push_het();
inducer_jump=1;



%% STATE default parameter values
defpar=sim.het.parameters;

% integral controller parameters
defpar('K_dna-sens')=4000;
defpar('eta_dna-sens')=1;
defpar('kb_anti')=300;
defpar('a_sens')=15;
defpar('a_anti')=2000;
defpar('a_act')=1000;

% wite these concentrations into the simulator
for k=keys(defpar)
    sim.het.parameters(k{1})=defpar(k{1});
end
sim=sim.push_het();

% HENCE CALCULATE more informative parameters
infopar=containers.Map('KeyType', 'char','ValueType', 'double');
infopar('u')=(defpar('c_act').*defpar('a_act'))./(defpar('c_anti').*defpar('a_anti')); % input, i.e. ratio of sense to antisense transcription
infopar('kappa')=(defpar('c_anti').*defpar('a_anti')); % integral gain, i.e. antisense transcription rate
infopar('theta')=defpar('kb_anti'); % binding rate of sense and antisense
infopar('zeta')=defpar('c_sens').*defpar('a_sens'); % transcription rate of sensor dna
infopar('K')=defpar('K_dna-sens');

%% DEFINE how we will vary different parameters

% each heatmap involves kappa and theta varied over a wide range
kappa_vals=infopar('kappa')*logspace(-2,2,7);
theta_vals=infopar('theta')*logspace(-2,2,7);
% kappa_vals=[infopar('kappa')];
% theta_vals=[infopar('theta')];

% heatmaps made for K and zeta varied over a wide range
% zeta_vals=infopar('zeta')*logspace(-1,1,3);
% K_vals=infopar('K')*logspace(-1,1,3);
zeta_vals=[infopar('zeta')];
K_vals=[infopar('K')];

%% RUN with varied parameters

%initialise results arrays
d0_errors={};
step_changes={};
step_mitigations={};
l0_saved={};
for i_zeta=1:size(zeta_vals,2)
    d0_errors{i_zeta}={};
    step_changes{i_zeta}={};
    step_mitigations{i_zeta}={};
    l0_saved{i_zeta}={};
end

% run
for i_zeta=1:size(zeta_vals,2)
    for i_K=1:size(K_vals,2)
        disp(['Simulating for i_zeta=,i_K= ',num2str([i_zeta i_K])])
        % reset parameters to default
        sim.het.parameters=defpar;
        
        % change parameters according to K and zeta values
        sim.het.parameters('c_sens')=1;
        sim.het.parameters('a_sens')=zeta_vals(i_zeta); % this way, product equals zeta for sure
        sim.het.parameters('K_dna-sens')=K_vals(i_K);

        sim=sim.push_het();

        % run
        out=test(sim,kappa_vals,theta_vals,inducer_jump,infopar('u'),Delta,Max_iter);
        
        % record results
        d0_errors{i_zeta}{i_K}=out{1};
        step_changes{i_zeta}{i_K}=out{2};
        step_mitigations{i_zeta}{i_K}=out{3};
        l0_saved{i_zeta}{i_K}=out{4};
    end
end
%% SAVE RESULTS
% results.d0_errors=d0_errors;
% results.step_changes=step_changes;
% results.step_mitigations=step_mitigations;
% results.l0_saved=l0_saved;
% save('integral_results/test_integral.mat','results')

%% LOAD RESULTS
% close all
% clear results d0_errors step_changes step_mitigations l0_saved
% loaded_data=load('integral_results/test_integral_wth.mat');
% results=loaded_data.results;
% 
% d0_errors=results.d0_errors;
% step_changes=results.step_changes;
% step_mitigations=results.step_mitigations;
% l0_saved=results.l0_saved;

%% Figure 9c - PLOT d0_errors - steady-state errors of F_{sens}
Fig9c=figure('Position',[0 0 600 350]);
set(Fig9c, 'defaultAxesFontSize', 9)

for i_zeta=1:size(zeta_vals,2)
    for i_K=1:size(K_vals,2)
        subplot(size(zeta_vals,2),size(K_vals,2),(size(zeta_vals,2)-i_zeta).*size(zeta_vals,2)+i_K)
        
        hmap=heatmap(theta_vals,flip(kappa_vals),flip(d0_errors{i_zeta}{i_K},1)); % make heatmap

        caxis([0 15]) % colours axis
        
        % CHANGING AXIS LABELS
        custom_x_lables={};
        hmap.XDisplayLabels = string(round(theta_vals,2));
        hmap.YDisplayLabels = string(flip(round(kappa_vals/1000,2)));

        title('Error between u and F_{anti} before disturbance [%]');
        ylabel('\kappa, annihilator synthesis rate [nM/h]')
        xlabel('\theta, actuator-annihilator binding rate constant [1/(nM \cdot h)]')

        
    end
end

%% DEFINE the METRIC FUNCTION
% sum-of-squares function; requires post-treatment
function out=test(sim,kappa_vals,theta_vals,inducer_jump,u,Delta,Max_iter)
    % initilaise output arrays
    d0_err=zeros(size(kappa_vals,2),size(theta_vals,2)); % SS error for p_sens signal with no disturbance (%)
    step_chan=zeros(size(kappa_vals,2),size(theta_vals,2)); % change in growth rate upon disturbance (%)
    step_mit=zeros(size(kappa_vals,2),size(theta_vals,2)); % reduction of growth rate change due to AIF (%)
    l0_saved=zeros(size(kappa_vals,2),size(theta_vals,2)); % SS growth rate with no disturbance (%) - for reference

    for i_kappa=1:size(kappa_vals,2)
        for i_theta=1:size(theta_vals,2)
            progress=((i_kappa-1).*size(kappa_vals,2)+i_theta-1)./(size(kappa_vals,2).*size(theta_vals,2)).*100;
            disp(['Progress: ',num2str(progress),' %'])

            % change parameters according to kappa and theta values
            sim.het.parameters('c_anti')=1;
            sim.het.parameters('a_anti')=kappa_vals(i_kappa); % this way, product equals kappa for sure
            sim.het.parameters('c_act')=1;
            sim.het.parameters('a_act')=kappa_vals(i_kappa).*u; % this way, product equals kappa for sure
            sim.het.parameters('kb_anti')=theta_vals(i_theta);

            % enable integral feedback
            sim.het.parameters('freeze_fbk')=0;
            sim=sim.push_het;

            % evaluate steady state without induction
            sim.ext.input_func_parameters('inducer_level')=0;
            sstss=get_steady(sim,Delta,Max_iter); % evaluate steady state value
            ss0=sstss{1};
            [F0,l0]=get_fl(ss0,sim);
            
            % for next stage, set initial conditions as required
            d0_init=containers.Map('KeyType', 'char','ValueType', 'double');
            d0_init('m_a')=ss0(1);
            d0_init('m_r')=ss0(2);
            d0_init('p_a')=ss0(3);
            d0_init('R')=ss0(4);
            d0_init('tc')=ss0(5);
            d0_init('tu')=ss0(6);
            d0_init('Bcm')=ss0(7);
            d0_init('s')=ss0(8);
            d0_init('h')=ss0(9);
            x_het=ss0(10 : (9+2*sim.num_het) );
            for j=1:sim.num_het
                % mRNA
                d0_init(['m_',sim.het.names{j}])=x_het(j);
                % protein
                d0_init(['p_',sim.het.names{j}])=x_het(sim.num_het+j);
            end
            
            % evaluate steady state now with jump in concentration
            sim.ext.input_func_parameters('inducer_level')=inducer_jump;
            sstss=get_steady(sim,Delta,Max_iter); % evaluate steady state value
            ss_jump=sstss{1};
            [F_jump,l_jump]=get_fl(ss_jump,sim);
            
            % evaluate steady state now with jump in concentration, but the feedback is frozen
            sim.ext.input_func_parameters('inducer_level')=inducer_jump;
            sim.het.parameters('freeze_fbk')=1;
            sim.het.parameters('fixed_F')=F0;
            sstss=get_steady(sim,Delta,Max_iter); % evaluate steady state value
            ss_jump=sstss{1};
            [F_freeze,l_freeze]=get_fl(ss_jump,sim);
            
            % CALCULATE OUTPUTS
            d0_err(i_kappa,i_theta)=abs((F0-u)./u).*100; % error in sensor signal (%)
            %disp([F0 u]);
            step_chan(i_kappa,i_theta)=abs((l_jump-l0)./l0).*100; % error in sensor signal (%)
            step_mit(i_kappa,i_theta)=((l0-l_freeze)-(l0-l_jump))./(l0-l_freeze).*100;
            l0_saved(i_kappa,i_theta)=l0;

            disp([theta_vals(i_theta), d0_err(i_kappa,i_theta)])
        end

        out{1}=d0_err;
        out{2}=step_chan;
        out{3}=step_mit;
        out{4}=l0_saved;
    end

    
end

function [F_anti,l]=get_fl(ss,sim)
    % get growth rate and ribosome mass fraction
    par=sim.parameters;
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

    e=sim.coll.e(par,tc); % translation elongation rate
    kcmh=par('kcm').*h; % ribosome inactivation rate due to chloramphenicol

    % ribosome dissociation constants
    k_a=sim.coll.k(e,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
    k_r=sim.coll.k(e,par('k+_r'),par('k-_r'),par('n_r'),kcmh);
    k_het=ones(sim.num_het,1);
    % heterologous genes
    if(sim.num_het>0)
        for j=1:sim.num_het
            k_het(j)=sim.coll.k(e,...
            sim.parameters(['k+_',sim.het.names{j}]),...
            sim.parameters(['k-_',sim.het.names{j}]),...
            sim.parameters(['n_',sim.het.names{j}]),...
            kcmh);
        end
    end

    D=1+(m_a./k_a+m_r./k_r+sum(ss_het(1:sim.num_het)./k_het))./...
        (1-par('phi_q')); % denominator in ribosome competition csim.optalculations
    B=R.*(1-1./D); % actively translating ribosomes - INCLUDING Q

    % growth rate
    l=sim.coll.l(par,e,B);

    F_anti=sim.het.regulation('anti',ss,sim.ext.input(ss,0)); % t=0 passed to ext.input as it's time invariant anyhow
end