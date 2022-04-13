%% Function for getting desired steady state values
function [l,e,phi_het]=get_lephihet(sim,ss)
    % PAREMETER MAP TO VARIABLE PAR
    par = sim.parameters;
    
    % STATE VECTOR TO SINGLE VARIABLES
    m_a = ss(1);
    m_r = ss(2);
    p_a = ss(3);
    R = ss(4);
    tc = ss(5);
    tu = ss(6);
    Bcm = ss(7);

    % USEFUL PRE-CALCULATIONS
    % translation elongation rate
    e=sim.coll.e(par,tc);

    % ribosome inactivation rate due to chloramphenicol
    kcmh=par('kcm').*par('h');

    % ribosome dissociation constants
    k_a=sim.coll.k(e,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
    k_r=sim.coll.k(e,par('k+_r'),par('k-_r'),par('n_r'),kcmh);

    % consider the contribution of heterologous genes to resource competition if there are any
    mk_het=0;
    if(sim.het.num_genes>0)
        ss_het=ss(8:end); % current state of heterologous genes
        sim.het = sim.het.find_current_ks(e,kcmh); % ribosome dissociation constants for heterolgous genes
        for i=1:sim.het.num_genes
            mk_het = mk_het + ss_het(i)./sim.het.current_ks(i); % get sum of mx/kx for heterologous genes
        end
    end

    T=tc./tu; % ratio of charged to uncharged tRNAs
    D=1+(m_a./k_a+m_r./k_r+mk_het)./(1-par('phi_q')); % denominator in ribosome competition cobj.optalculations
    B=R.*(1-1./D); % actively translating ribosomes - INCLUDING Q

    % tRNA charging and synthesis
    nu=sim.coll.nu(par,tu);
    psi=sim.coll.psi(par,T);

    % growth rate
    l=sim.coll.l(par,e,B);

    % heterologous protein mass fraction
    phi_het=0;
    for i=1:sim.het.num_genes
        phi_het=phi_het+ss_het(sim.het.num_genes+i).*sim.het.parameters(['n_',sim.het.names{i}])./sim.parameters('M');
    end
end