%% cell_params.m
% values of all parameters decsribing the host cell

%%

function params=cell_params(~)
    params = containers.Map('KeyType', 'char', ...
                'ValueType', 'double');
    
    % GENERAL PARAMETERS
    params('M') = 11.9*10^8; % cell mass (aa) - taken for 1 div/h for order of magnitude [3]
    params('phi_q')=0.59; % constant housekeeping protein mass fraction [4]

    % GENE EXPRESSION PARAMETERS  
    % metabolic/aminoacylating genes
    params('c_a') = 1; % copy no. (nM) - convention
    params('b_a') = 6; % mRNA decay rate (/h) [1]
    params('k+_a') = 60; % ribosome binding rate (/h/nM) [1]
    params('k-_a') = 60; % ribosome unbinding rate (/h) [1]
    params('n_a') = 300; % protein length (aa) [1]

    % ribosome gene
    params('c_r') = 1; % copy no. (nM) - convention
    params('b_r') = 6; % mRNA decay rate (/h) [1]
    params('k+_r') = 60; % ribosome binding rate (/h/nM) [1]
    params('k-_r') = 60; % ribosome unbinding rate (/h) [1]
    params('n_r') = 7459; % protein length (aa) [1]

    % ACTIVATION & RATE FUNCTION PARAMETERS
    params('e_max')=20*3600; % max elongation rate (aa/h) [2]
    params('psi_max')= 1080000; % max synthesis rate (aa/h) [2]   
    params('tau')= 1; % ppGpp sensitivity (ribosome transc. and tRNA synth. Hill const) [2]

    % optionally - used to switch to fixed ribosome transcription rate
    params('is_fixed_F_r') = 0; % 1 if fixed
    params('fixed_F_r') = 0.001; % fixed F_r value from 0 to 1

    % FITTED PARAMETERS
    params('a_a') = 2.706*10^5; % metabolic gene transcription rate (/h) 
    params('a_r') = 2.189*10^5; % ribosomal gene transcription rate (/h) 
    params('nu_max')= 4333.1; % max tRNA amioacylation rate (/h)
    params('K_nut')= 6475.7; % tRNA charging rate Michaelis-Menten constant (nM) 
    params('K_e')= 6475.7; % translation elongation rate Michaelis-Menten constant (nM) 
    params('kcm')= 4.2020/10000; % chloramphenical binding rate constant (/h/nM)   
end

%% REFERENCES:
% [1] - Weiße AY et al. 2015 Mechanistic links between cellular trade-offs, gene expression, and growth
% [2] - Chure G et al. 2022 An Optimal Regulation of Fluxes Dictates Microbial Growth In and Out of Steady-State
% [3] - Bremer H et al. 2008 Modulation of Chemical Composition and Other Parameters of the Cell at Different Exponential Growth Rates
% [4] - Hui et al. 2015 Quantitative proteomic analysis reveals a simple strategy of global resource allocation in bacteria