%% cell_init_conds.m
% default initial conditions for the host cell

%%
function init_conds=cell_init_conds(par)
    init_conds = containers.Map('KeyType', 'char', ...
            'ValueType', 'double');
    
    % mRNAs
    init_conds('m_a')=0; % metabolic
    init_conds('m_r')=0; % ribosomal
    
    % proteins
    init_conds('p_a')=par('M').*(1-par('phi_q'))./(2.*par('n_a')); % metabolic *
    init_conds('R')=par('M').*(1-par('phi_q'))./(2.*par('n_r')); % ribosomal *

    % tRNAs
    init_conds('tc')=2.*80000; % charged tRNAs **
    init_conds('tu')=2.*80000; % free tRNAs **

    % ribosomes inactivated by chloramphenicol
    init_conds('Bcm')=0;

    % nutrient quality s and chloramphenicol conc. h
    init_conds('s')=0.5;
    init_conds('h')=0; % no translation inhibition
    
end

%% NOTES
% * Start with 50/50 a/R allocation as a convention
% ** 3E-5 abundance units in Chure et al., 2022 = 80 uM = 80000 nM