%% dram_fit.m
% fit parameters to experimental data using DRAM

%% CLEAR parameters, add paths of all files

addpath(genpath('..'))

clear
close all

timer=tic();

%% DEFINE starting parameter values

% All parameters are constrained to be positive. Initial values from Weisse et al. 2015 and Chure et al. 2022
params = {
    {'a_a', 55800,  0} % metabolic gene transcription rate
    {'a_r', 55800, 0} % max. ribosomal gene transcription rate
    {'nu_max', 6000,  0} % max. tRNA aminoacylatio rate
    {'K_t', 80000, 0} % MM constants for translation elongation and tRNA charging rates
    {'kcm', 0.3594/1000, 0} % chloramphenicol binding rate constant
    };

%% LOAD Experimental data

% read the experimental dataset (eq2 strain of Scott 2010)
dataset = readmatrix('data/scott2010_chure2022_notext.csv');

% nutrient qualities are equally log-spaced points
nutr_quals=logspace(log10(0.08),log10(0.5),6);

% get inputs: nutrient quality and h; get outputs: l and rib mass frac
data.xdata=[]; % initialise inputs array
data.ydata=[]; % intialise outputs array
for i = 1:size(dataset,1)
    if(dataset(i,1)>0.5)
        % inputs
        nutr_qual = nutr_quals(fix((i-1)/5)+1); % records start from worst nutrient quality
        h = dataset(i,4)*1000; % all h values for same nutr quality same go one after another. Convert to nM from uM!
        data.xdata=[data.xdata; [nutr_qual,h]];
    
        % outputs
        l = dataset(i,1); % growth rate (1/h)
        phi_r = dataset(i,3); % ribosome mass fraction
        data.ydata=[data.ydata; [l,phi_r]];
    end
end


%% SET UP the simulator

sim=cell_simulator; % initialise simulator

% parameters for getting steady state
sim.tf = 10; % single integraton step timeframe
Delta = 0.1; % threshold that determines if we're in steady state
Max_iter = 75; % maximum no. iterations (checking if SS reached over first 750 h)

sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9); % more lenient integration tolerances for speed

model.ssfun = @(theta,data) rc_ecoli_sos(theta,data,sim,Delta,Max_iter);

%% DEFINE prior distirbution parameters for errors
model.S20 = [0.0319 0.0142];
model.N0  = [7 8];

%% BURN-IN

options.nsimu = 1000;
[results, ~, ~]= mcmcrun(model,data,params,options);

%% MCMC SAMPLING

options.nsimu = 10000;
[results, chain, s2chain] = mcmcrun(model,data,params,options, results);


%% RECORD results in a .mat file

mcmc_outcome.results=results;
mcmc_outcome.chain=chain;
mcmc_outcome.s2chain=s2chain;
% save results
save('outcomes/mcmc_outcome.mat','mcmc_outcome')

%% DISPLAY RUNTIME
disp(toc(timer))