%% CLEAR

addpath(genpath('..'))

clear
close all

%% ORIGINAL PARAMETER VALUES
params = {
    {'a_a', 55800,  0}
    {'a_r', 55800, 0}
    {'nu_max', 6000,  0}
    {'K_t', 80000, 0}
    {'kcm', 0.3594/1000, 0}
    };
% record original params into a single vector
theta_origin=zeros([size(params,1) 1]);
for i=1:size(theta_origin,1)
    theta_origin(i) = params{i}{2};
end

%% LOAD EXPERIMENTAL DATA (FOR COMPARISON WITH MODEL)
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

%% SET UP SIMULATOR
% params for finding steady state
Delta = 0.1; % threshold below which the changes in l and phi_r assumed negligible
Max_iter=25; % max. no iterations for finding steady state
sim=advanced_simulator; % initialise simulator
sim.tf=10;
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9);

%% DISPLAY CHAIN PLOTS
% Chain plots should reveal that the chain has converged and we can
% use the results for estimation and predictive inference.

% first 5000 points
load('outcomes/newcm_outcome_1306.mat');
results = mcmc_outcome.results;
chain = mcmc_outcome.chain;
s2chain = mcmc_outcome.chain;

% load('outcomes/azure_outcome_1805.mat');
% results1 = mcmc_outcome.results;
% chain1 = mcmc_outcome.chain;
% s2chain1 = mcmc_outcome.chain;
% chain=[chain;chain1];
% s2chain=[s2chain;s2chain1];

% discard unstable entries
%chain=chain(1:1500,:);


figure(2); clf
mcmcplot(chain,[],results,'pairs');
figure(3); clf
mcmcplot(chain,[],results,'denspanel',2); %

%% DISPLAY MCMC CHAIN STATISTICS
% Function |chainstats| calculates mean ans std from the chain and
% estimates the Monte Carlo error of the estimates. Number |tau| is
% the integrated autocorrelation time and |geweke| is a simple test
% for a null hypothesis that the chain has converged.
chainstats(chain,results)

%% ESTIMATE PROBABILITY DISTRIBUTION MODES
pd = fitdist(chain(:,1),'Kernel');
points_in_pdf=1000; % resolution of pdf approximation
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
means=mean(chain,1).';

% compare means and modes
% reduced_sos(modes,data,sim,Delta,Max_iter)
% reduced_sos(means,data,sim,Delta,Max_iter)

%% MODEL PREDICTIONS
% get points from the model
%modes=[248.4; 55800/250; 4000; 1080000; 80000; 1; 0.5/1000];
% modes=exp([ 14.0574      13.3057      8.74011      11.2213     -7.87053].');

% modes(1) = 55800; %11360; %194.1; % transcription rate (/h) [FITTED]
% modes(2) = 55800; %9077; %155.1; % transcription rate (/h) [FITTED]
% modes(5)=modes(5)*10;
% modes(1) = modes(1).*10000; %11360; %194.1; % transcription rate (/h) [FITTED]
% modes(2) = modes(2).*10000; %9077; %155.1; % transcription rate (/h) [FITTED]
% modes(2)=modes(1)/1.25;
% modes(4)=modes(4)/3;

ymodel=advanced_modelfun(modes,data.xdata,sim,Delta,Max_iter);

sum((ymodel-data.ydata).^2)

% get point from ORIGINAL model estimate
ymodel_origin=advanced_modelfun(theta_origin,data.xdata,sim,Delta,Max_iter);

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
    plot(ymodel(i,1),ymodel(i,2),'+','Color',colours(colourind),'MarkerSize',8,'LineWidth',1.25) % model prediction
    %plot(ymodel_origin(i,1),ymodel_origin(i,2),'x','Color',colours(colourind),'MarkerSize',8,'LineWidth',1.25) % original model predictions

end
ylabel('Ribosome mass fraction \phi_r');
xlabel('Growth rate \lambda, 1/h')
xlim([0.25 1.75])
% ylim([0 0.45])
hold off
% legend('RDM+Glucose, \sigma=0.5','RDM+Glycerol, \sigma=0.347','cAA+Glucose, \sigma=0.240',...
%     'cAA+Glycerol, \sigma=0.167','M63+Glucose, \sigma=0.115','FontSize',9)

%% 1ST GROWTH LAW FITS
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
%legend('4 uM chlorampenicol','2 uM chlorampenicol','0 uM chlorampenicol','FontSize',9)

%% 2ND GROWTH LAW FITS
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
% estimate Fisher Info Matrix as inverse of covar
FIM=pinv(cov(chain))
eigenvals=eig(FIM);
eigenvals_norm=eigenvals./max(eigenvals);
eigenvals_log=log10(eigenvals_norm)

%% PARAMETER SENSITIVITY
% get the matrix C that diagonalises FIM
[C,DIAG]=eig(FIM);

% get parameter sensitivities
sens2=zeros([1 size(params,1)]);
for j=1:size(sens2,2)
    for i=1:size(sens2,2)
        sens2(j)=sens2(j)+eigenvals(i)*(C(i,j)^2);
    end
end
sens_log=log10(sens2./(sum(sens2)))