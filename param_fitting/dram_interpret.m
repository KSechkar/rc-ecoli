%% dram_interpret.m
% load and interpret rsults of parameter fitting with DRAM

%% CLEAR parameters, add paths of all files

addpath(genpath('..'))

clear
close all

%% DEFINE starting parameter values (to compare with the fit)

params = {
    {'a_a', 55800,  0} % metabolic gene transcription rate
    {'a_r', 55800, 0} % max. ribosomal gene transcription rate
    {'nu_max', 6000,  0} % max. tRNA aminoacylatio rate
    {'K_t', 80000, 0} % MM constants for translation elongation and tRNA charging rates
    {'kcm', 0.3594/1000, 0} % chloramphenicol binding rate constant
    };

% record original params into a single vector
theta_origin=zeros([size(params,1) 1]);
for i=1:size(theta_origin,1)
    theta_origin(i) = params{i}{2};
end

%% LOAD Experimental data (to compare with the fit)

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

%% LOAD ftting outcome

% first 5000 points
load('outcomes/hpz2_outcome_2705.mat');
results = mcmc_outcome.results;
chain = mcmc_outcome.chain;
s2chain = mcmc_outcome.chain;

%% PLOT the chains

% pairwise parameter value scatter diagram
figure(2); clf
mcmcplot(chain,[],results,'pairs');

% obtained probability density functions
figure(3); clf
mcmcplot(chain,[],results,'denspanel',2); %

%% DISPLAY mcmc chain statistics
% Function |chainstats| calculates mean ans std from the chain and
% estimates the Monte Carlo error of the estimates. Number |tau| is
% the integrated autocorrelation time and |geweke| is a simple test
% for a null hypothesis that the chain has converged.
chainstats(chain,results)

%% ESTIMATE modes of probaility distributions to get fitted parameter values

pd = fitdist(chain(:,1),'Kernel');
points_in_pdf=500; % resolution of pdf approximation
modes=zeros([size(chain,2) 1]); % initialise array of mode values
pdfs=zeros([size(chain,2) points_in_pdf]);
for i = 1:size(chain,2)
    points = linspace(0.8*min(chain(:,i)),1.2*max(chain(:,i)),points_in_pdf);
    [pdf,points]=density(chain(:,i),points,2); % estimate pdf on discrete points on x-axis
    [~,mode_index]=max(pdf); % get position of mode
    modes(i)=points(mode_index); % get mode

    % save pdf
    for j=1:size(points,2)
        pdfs(i,j)=pdf(j);
    end
end

%% GET model predictions with fitted parameters
ymodel=rc_ecoli_modelfun(modes,data.xdata,sim,Delta,Max_iter);

disp(['SOS=',num2str(sum((ymodel-data.ydata).^2))]) % print resultant sum of squared errors

% get points for original/strarting parameter values (optional)
% ymodel_origin=rc_ecoli_modelfun(theta_origin,data.xdata,sim,Delta,Max_iter);

%% PLOT AND COMPARE

Fdata = figure('Position',[0 0 385 280]);
set(Fdata, 'defaultAxesFontSize', 9)
hold on
colours=['c','g','b','r','m'];
colourind=1;
last_nutr_qual=data.xdata(1,1);
for i=1:size(data.xdata,1)
    if(data.xdata(i,1)~=last_nutr_qual)
        colourind=colourind+1;
        last_nutr_qual=data.xdata(i,1);
    end
    plot(data.ydata(i,1),data.ydata(i,2),'o','Color',colours(colourind),'LineWidth',1) % real data
    plot(ymodel(i,1),ymodel(i,2),'+','Color',colours(colourind),'MarkerSize',8,'LineWidth',1.25) % model predictions
    %plot(ymodel_origin(i,1),ymodel_origin(i,2),'x','Color',colours(colourind),'MarkerSize',8,'LineWidth',1.25) % original model predictions (optional)

end
ylabel('Ribosome mass fraction \phi_r');
xlabel('Growth rate \lambda, 1/h')
xlim([0.25 1.75])
hold off

%% ADD lines for 1ST GROWTH LAW FITS

% group data by nutrient quality
xs_1={[],[],[]};
ys_1={[],[],[]};
nutrind=1;
chlorind=1;
last_nutr_qual=data.xdata(1,1);
for i=1:size(data.xdata,1)
    if(data.xdata(i,1)~=last_nutr_qual)
        nutrind=nutrind+1;
        chlorind=1;
        last_nutr_qual=data.xdata(i,1);
    end
    xs_1{chlorind}(end+1)=ymodel(i,1);
    ys_1{chlorind}(end+1)=ymodel(i,2);
    chlorind=chlorind+1;
end

% make linear fits
fit_coeffs=zeros([nutrind 2]);
for nutrind=1:size(xs_1,2)
    if(size(xs_1{nutrind},2)>=2)
        linfit=polyfit(xs_1{nutrind},ys_1{nutrind},1);
        fit_coeffs(nutrind,1)=linfit(1);
        fit_coeffs(nutrind,2)=linfit(2);
    end
end

% plot
dashings={'--','-.',':'};
hold on
for chlorind=size(xs_1,2):(-1):1
    if(fit_coeffs(chlorind,1)~=0)
        xpoints=linspace(0,xs_1{chlorind}(end)*1.15,100); % points for which we plot the linear fit
        ypoints=polyval(fit_coeffs(chlorind,:),xpoints);
        plot(xpoints,ypoints,'Color','k','LineStyle',dashings{chlorind},'LineWidth',0.5)
    end
end
xlim([0 2])
ylim([0 0.45])
hold off

%% ADD lines for 2ND GROWTH LAW FITS

% group data by nutrient quality
xs_2={[]};
ys_2={[]};
nutrind=1;
chlorind=1;
last_nutr_qual=data.xdata(1,1);
for i=1:size(data.xdata,1)
    if(data.xdata(i,1)~=last_nutr_qual)
        nutrind=nutrind+1;
        chlorind=1;
        last_nutr_qual=data.xdata(i,1);
        xs_2{nutrind}=[];
        ys_2{nutrind}=[];
    end
    xs_2{nutrind}(chlorind)=ymodel(i,1);
    ys_2{nutrind}(chlorind)=ymodel(i,2);
    chlorind=chlorind+1;
end

% make linear fits
fit_coeffs=zeros([nutrind 2]);
for nutrind=1:size(xs_2,2)
    if(size(xs_2{nutrind},2)>=2)
        linfit=polyfit(xs_2{nutrind},ys_2{nutrind},1);
        fit_coeffs(nutrind,1)=linfit(1);
        fit_coeffs(nutrind,2)=linfit(2);
    end
end

% plot
hold on
for nutrind=1:size(xs_2,2)
    if(fit_coeffs(nutrind,1)~=0)
        xpoints=linspace(0,xs_2{nutrind}(1)*1.1,100); % points for which we plot the linear fit
        ypoints=polyval(fit_coeffs(nutrind,:),xpoints);
        plot(xpoints,ypoints,'Color',colours(nutrind),'LineWidth',0.5)
    end
end
xlim([0 1.9])
ylim([0 0.45])
hold off

%% FIM AND EIGENVALUES
% estimate Fisher Info Matrix as inverse of the variance-covariance matrix for
% the chains
FIM=pinv(cov(chain)) % print FIM
eigenvals=eig(FIM);
eigenvals_norm=eigenvals./max(eigenvals);
eigenvals_log=log10(eigenvals_norm) % print normalised eigenvalue logarithms

%% PARAMETER SENSITIVITIES
% get the matrix C that diagonalises FIM
[C,DIAG]=eig(FIM);

% get parameter sensitivities
sens2=zeros([1 size(params,1)]);
for j=1:size(sens2,2)
    for i=1:size(sens2,2)
        sens2(j)=sens2(j)+eigenvals(i)*(C(i,j)^2);
    end
end
disp('Parameter sensitivities (for a_a, a_r, nu_max, K_t and kcm respectively) ')
sens_log=log10(sens2./(sum(sens2))) % print parameter sensitivities