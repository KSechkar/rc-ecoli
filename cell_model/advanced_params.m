function params=advanced_params(~)
    params = containers.Map('KeyType', 'char', ...
                'ValueType', 'double');
    
    % GENERAL PARAMETERS
    params('M') = 11.9*10^8; % cell mass (aa) - taken for 1 div/h for order of magnitude [3]

    % GENE EXPRESSION PARAMETERS
    % housekeeping genes
    params('phi_q')=0.59; %0.55; % constant housekeeping protein fraction [4]
    
    % metabolic/aminoacylating genes
    params('c_a') = 1; % copy no. (nM) [1]*
    params('b_a') = 6; % mRNA decay rate (/h) [1]
    params('k+_a') = 60; % ribosome binding rate (/h/nM) [1]
    params('k-_a') = 60; % ribosome unbinding rate (/h) [1]
    params('n_a') = 300; % protein length (aa) [1]

    % ribosome gene
    params('c_r') = 1; % copy no. (nM) [1]*
    params('b_r') = 6; % mRNA decay rate (/h) [1]
    params('k+_r') = 60; % ribosome binding rate (/h/nM) [1]
    params('k-_r') = 60; % ribosome unbinding rate (/h) [1]
    params('n_r') = 7459; % protein length (aa) [1]

    % ACTIVATION & RATE FUNCTION PARAMETERS
    params('e_max')=20*3600; % max elongation rate (aa/h) [2]
    params('K_nus')=0.5*10^6; % substrate-wise Hill constant (nM), Monod const for glucose [2]
    params('psi_max')= 1080000; % max synthesis rate (aa/h) [2]   
    params('tau')= 1; % ppGpp sensitivity (ribosome transc. and tRNA synth. Hill const) [2]

    % optionally -switch to fixed ribosome transcription rate
    params('is_fixed_F_r') = 0; % 1 if fixed
    params('fixed_F_r') = 0.001; % fixed F_r value from 0 to 1

    %CHLORAMPHENICOL CONCENTRATION
    params('h') = 0; % chloramphenicol concentration (nM)

    % PARAMETERS TO BE FIT
    params('a_a') = 2.706*10^5; %11360; %194.1; % transcription rate (/h) [FITTED]
    params('a_r') = 2.189*10^5; %9077; %155.1; % transcription rate (/h) [FITTED]
    params('nu_max')= 4333.1; %3662.2; % max metabolic rate (/h) [FITTED]
    params('K_nut')= 6475.7; %12793; % tRNA charging rate Hill constant (nM) [FITTED]
    params('K_e')= 6475.7; % translation elongation rate Hill constant (nM) [FITTED]
    params('kcm')= 4.2020/10000; %0.0003148; %0.3594/1000; % chloramphenical binding rate constant (/h/nM) [FITTED]    
end

% SOURCES
% [1] - Weiße AY et al. Mechanistic links between cellular trade-offs, gene expression, and growth.
% [2] - Chure G et al. An Optimal Regulation of Fluxes Dictates Microbial Growth In and Out of Steady-State
% [3] - Bremer H et al. Modulation of Chemical Composition and Other Parameters of the Cell at Different Exponential Growth Rates
% [4] - Hui et al. Quantitative proteomic analysis reveals a simple strategy of global resource allocation in bacteria