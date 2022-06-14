function init_conds=advanced_init_conds(par)
    init_conds = containers.Map('KeyType', 'char', ...
            'ValueType', 'double');
    
    % mRNAs
    init_conds('m_a')=0;
    init_conds('m_r')=0;
    % proteins
    init_conds('p_a')=par('M').*(1-par('phi_q'))./(2.*par('n_a')); % [2] *
    init_conds('R')=par('M').*(1-par('phi_q'))./(2.*par('n_r')); % [2] *
    % tRNAs
    init_conds('tc')=2.*80000; % charged tRNAs [2] **
    init_conds('tu')=2.*80000; % free tRNAs [2] **
    % ribosomes inactivated by chloramphenicol
    init_conds('Bcm')=0;
    % nutrient quality and chloramphenicol
    init_conds('s')=0.5;
    init_conds('h')=0;
    
end

% SOURCES
% [1] - Weiße AY et al. Mechanistic links between cellular trade-offs, gene expression, and growth.
% [2] - Chure G et al. An Optimal Regulation of Fluxes Dictates Microbial Growth In and Out of Steady-State

% NOTES
% * In [2] (see advanced_params.m for article name), we start with 50/50 a/R allocation
% ** 3E-5 abundance units = 80 uM = 80000 nM