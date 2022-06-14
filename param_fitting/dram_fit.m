%% CLEAR

addpath(genpath('..'))

clear
close all

timer=tic();
%% LOAD EXPERIMENTAL DATA
% read the experimental dataset (eq2 strain of Scott 2010)
dataset = readmatrix('data/scott2010-eq2-notext.csv');

% [1] => nutrient qualities are equally log-spaced points
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

%% PLOT EXPERIMENTAL DATA
% Fdata = figure('Position',[0 0 540 480]);
% hold on
% plot(data.ydata(:,1),data.ydata(:,2),'o','Color','r')
% plot([0.5 0.5],[0 0.3],'--','Color','k')
% title('Experimental data being fitted');
% ylabel('Ribosome mass fraction');
% xlabel('Growth rate, 1/h')
% xlim([0 2])
% ylim([0 0.3])
% hold off

%% DEFINE SUM OF SQUARES FUNCTION
sim=advanced_simulator; % initialise simulator

% parameters for getting steady state
sim.tf = 10; % single integraton step timeframe
Delta = 0.1; % threshold that determines if we're in steady state
Max_iter = 75; % maximum no. iterations (checking if SS reached over first 750 h)
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9);

model.ssfun = @(theta,data) advanced_sos(theta,data,sim,Delta,Max_iter);

%% DEFINE PARAMETERS
% All parameters are constrained to be positive. Initial values from Weisse
% 2015 and Chure 2022
params = {
    {'a_a', 4.14*60,  0}
    {'a_r', 930*60, 0}
    {'nu_max', 6000,  0}
    {'psi_max', 1080000,   0}
    {'K_t', 80000, 0}
    {'tau', 1, 0}
    {'kcm', 0.3594/1000, 0}
    };

%% DEFINE PRIOR FOR ERROR TERM
% The default prior distribution is sigma2 ~ invchisq(S20,N0), the inverse chi
% squared distribution (Gelman et al., 1995).
% Scale and degrees of freedom estiated from measurement errors by Scott et
% al., 2010
model.S20 = [0.0319 0.0142];
model.N0  = [7 8];

%%
% First generate an initial chain.
options.nsimu = 1000;
[results, ~, ~]= mcmcrun(model,data,params,options);
%%
% Then re-run starting from the results of the previous run.
options.nsimu = 10000;
[results, chain, s2chain] = mcmcrun(model,data,params,options, results);

mcmc_outcome.results=results;
mcmc_outcome.chain=chain;
mcmc_outcome.s2chain=s2chain;
% save results
save('outcomes/mcmc_outcome.mat','mcmc_outcome')

disp(toc(timer))