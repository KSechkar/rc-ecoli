%% AppFigure2.m
% Generate figures displaying fits/predictions obtained using a more
% accurate model of ribosome inactivation by chloramphenicol

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

%% SET UP SIMULATOR
% params for finding steady state
Delta = 0.1; % threshold below which the changes in l and phi_r assumed negligible
Max_iter=25; % max. no iterations for finding steady state
sim=cell_simulator_corrected_inactivation; % initialise simulator - NOTE: NOT CELL_SIMULATOR
sim.tf=10;
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9);


%% AppFigure2a - display results of the fitting the corrected model to experimental data

fitted_vals=[ 263.35; ... % transcription rate (/h)
    461.08; % transcription rate (/h)
    8559.9; % max metabolic rate (/h) 
    62319; % translation elongation and tRNA charging rate Hill constant (nM)
    4.0723/100000; % chloramphenical binding rate constant (/h/nM)  
    ];

ymodel=rc_ecoli_modelfun(fitted_vals,data.xdata,sim,Delta,Max_iter);

disp(['SOS=',num2str(sum((ymodel-data.ydata).^2))]) % print resultant sum of squared errors

% PLOT AND COMPARE
AppFig2a = figure('Position',[0 0 385 280]);
set(AppFig2a, 'defaultAxesFontSize', 9)
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

end
ylabel('\phi_r, ribosomal mass fraction');
xlabel('\lambda, growth rate [1/h]')

% 1ST GROWTH LAW FITS
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
xticks([0:0.5:1.5,1.9])
yticks([0:0.1:0.4,0.45])
grid on
%legend('4 uM chlorampenicol','2 uM chlorampenicol','0 uM chlorampenicol','FontSize',9)

% 2ND GROWTH LAW FITS
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



%% AppFigure2b - model predictions obtained with parameters describing the 

fitted_vals=[ 2.706*10^5; ... % transcription rate (/h)
    2.189*10^5; % transcription rate (/h)
    4333.1; % max metabolic rate (/h) 
    6475.7; % translation elongation and tRNA charging rate Hill constant (nM)
    4.2020/10000; % chloramphenical binding rate constant (/h/nM)  
    ];

ymodel=rc_ecoli_modelfun(fitted_vals,data.xdata,sim,Delta,Max_iter);

disp(['SOS=',num2str(sum((ymodel-data.ydata).^2))]) % print resultant sum of squared errors

% PLOT AND COMPARE
AppFig2b = figure('Position',[0 0 385 280]);
set(AppFig2b, 'defaultAxesFontSize', 9)
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

end
ylabel('\phi_r, ribosomal mass fraction');
xlabel('\lambda, growth rate [1/h]')
xlim([0.25 1.75])
% ylim([0 0.45])
hold off
% legend('RDM+Glucose, \sigma=0.5','RDM+Glycerol, \sigma=0.347','cAA+Glucose, \sigma=0.240',...
%     'cAA+Glycerol, \sigma=0.167','M63+Glucose, \sigma=0.115','FontSize',9)

% 1ST GROWTH LAW FITS
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
xticks([0:0.5:1.5,1.9])
yticks([0:0.1:0.4,0.45])
grid on
%legend('4 uM chlorampenicol','2 uM chlorampenicol','0 uM chlorampenicol','FontSize',9)

% 2ND GROWTH LAW FITS
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