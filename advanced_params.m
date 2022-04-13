function params=advanced_params(~)
    params = containers.Map('KeyType', 'char', ...
                'ValueType', 'double');
    
    % GENERAL PARAMETERS
    params('M') = 11.9*10^8; % cell mass (aa) - taken for 1 div/h for order of magnitude [3]
    params('s') = 0.5*10^6; % substrate concentration (nM) [WILL VARY]

    % GENE EXPRESSION PARAMETERS
    % housekeeping genes
    params('phi_q')=0.55; % constant housekeeping protein fraction [2]
    
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

    % optionally -switch to fixed ribosome transcription rate
    params('is_fixed_F_r') = 0; % 1 if fixed
    params('fixed_F_r') = 0.001; % fixed F_r value from 0 to 1

    % mostly, instead of Hill function for nutrient we just use 'quality coeff'
    params('using_nutr_qual') = 1; % 1 if using nutrient quality 
    params('nutr_qual') = 0.5;

    %CHLORAMPHENICOL CONCENTRATION
    params('h') = 0; % chloramphenicol concentration (nM)

    % PARAMETERS TO BE FIT
    params('a_a') = 194.6469; %181.4148; %4.14*60; % transcription rate (/h) [TO BE FIT]*
    params('a_r') = 155.1436; %118.8015; %930*60; % transcription rate (/h) [TO BE FIT]*
    params('nu_max')= 3662.2; %3690.7; %20.*params('n_a'); % max metabolic rate (/h) [TO BE FIT]
    params('psi_max')= 1200900; %1189500; %300*3600; % max synthesis rate (aa/h) [TO BE FIT]   
    params('K_t')= 12793; %6000; %1568.7; %80000; % elongation rate and tRNA charging rate Hill constant (nM) [TO BE FIT]
    params('tau')= 0.6450; % 0.6427; %1; % ppGpp sensitivity (ribosome transc. and tRNA synth. Hill const) [TO BE FIT]
    params('kcm')= 3.1323/10000; %0.0003148; %0.3594/1000; % chloramphenical binding rate constant (/h/nM) [TO BE FIT]    
end

% SOURCES
% [1] - Weiße AY et al. Mechanistic links between cellular trade-offs, gene expression, and growth.
% [2] - Chure G et al. An Optimal Regulation of Fluxes Dictates Microbial Growth In and Out of Steady-State
% [3] - Bremer H et al. Modulation of Chemical Composition and Other Parameters of the Cell at Different Exponential Growth Rates

% NOTES
% * Weisse et al. give c_x*a_x as a single variable (max transcription rate). So we just set c_x=1 and a_x to this value