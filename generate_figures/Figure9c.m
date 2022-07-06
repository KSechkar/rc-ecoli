%% Figure9c.m
% Explore parameter space for the actuator and annihilator RNA 
% synthesis and binding rates to combat leakiness, i.e. reduce the
% steady-state error in the value of the annhilitaor's transcription
% regulation function

%% CLEAR

addpath(genpath('..'))

clear
close all

%% SET UP the simulator
sim=cell_simulator; % initialise simulator
sim=sim.load_heterologous_and_external('aif_controller','constant_inducer');
sim.het.parameters('a_dist')=0; % we do not consider disturbance here
sim=sim.push_het(); % push heerologous genes' parameters

% error tolerances during integration
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9);

sim.tf=10;

% params for finding steady state
Delta=0.01;
Max_iter=75;

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

%% INITIALISE the array where the steady state errors will be recorded
d0_errors=zeros(size(kappa_vals,2),size(theta_vals,2))

%% RUN with different parameter combinations

for i_kappa=1:size(kappa_vals,2)
    for i_theta=1:size(theta_vals,2)
        progress=((i_kappa-1).*size(kappa_vals,2)+i_theta-1)./(size(kappa_vals,2).*size(theta_vals,2)).*100;
        disp(['Progress: ',num2str(progress),' %'])

        % change parameters according to kappa and theta values
        sim.het.parameters('c_anti')=1;
        sim.het.parameters('a_anti')=kappa_vals(i_kappa); % this way, product equals kappa for sure
        sim.het.parameters('c_act')=1;
        sim.het.parameters('a_act')=kappa_vals(i_kappa).*infopar('u'); % this way, product equals kappa for sure
        sim.het.parameters('kb_anti')=theta_vals(i_theta);

        % enable integral feedback
        sim.het.parameters('freeze_fbk')=0;
        sim=sim.push_het;

        % evaluate steady state without induction
        sim.ext.input_func_parameters('inducer_level')=0;
        ss0=get_steady(sim,Delta,Max_iter); % evaluate steady state value
        [F0,l0]=get_fl(ss0,sim);
        
        % calculate error in sensor signal (%)
        d0_errors(i_kappa,i_theta)=abs((F0-infopar('u'))./infopar('u')).*100;
    end
end


%% Figure 9c - PLOT d0_errors - steady-state errors of F_{sens}
Fig9c=figure('Position',[0 0 600 350]);
set(Fig9c, 'defaultAxesFontSize', 9)

hmap=heatmap(theta_vals,flip(kappa_vals),flip(d0_errors,1)); % make heatmap

caxis([0 15]) % colours axis

% CHANGING AXIS LABELS
custom_x_lables={};
hmap.XDisplayLabels = string(round(theta_vals,4));
hmap.YDisplayLabels = string(flip(round(kappa_vals/1000,2)));

title('Error between u and F_{anti} before disturbance [%]');
ylabel('\kappa, annihilator synthesis rate [nM/h]')
xlabel('\theta, actuator-annihilator binding rate constant [1/(nM \cdot h)]')

%% FUNCTION for getting the growth rate and the annihilator gene transcription regulation function
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

    e=sim.form.e(par,tc); % translation elongation rate
    kcmh=par('kcm').*h; % ribosome inactivation rate due to chloramphenicol

    % ribosome dissociation constants
    k_a=sim.form.k(e,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
    k_r=sim.form.k(e,par('k+_r'),par('k-_r'),par('n_r'),kcmh);
    k_het=ones(sim.num_het,1);
    % heterologous genes
    if(sim.num_het>0)
        for j=1:sim.num_het
            k_het(j)=sim.form.k(e,...
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
    l=sim.form.l(par,e,B);

    F_anti=sim.het.regulation('anti',ss,sim.ext.input(ss,0)); % t=0 passed to ext.input as it's time invariant anyhow
end