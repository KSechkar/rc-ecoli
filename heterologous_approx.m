%% heterologous_approx.m
% collection of functions for approximating the effects of heterologous
% gene expression

classdef heterologous_approx
   
    methods (Access = public)
        % setady state mRNA levels
        function m=ss_m(obj, ...
                plasmid_conc, ... % concentration of plasmid with heterologous genes
                F, ... % gene of interest (GOI) ss regulation function value
                c, ... % GOI concentration
                a, ... % GOI transcription rate
                b, ... % GOI mRNA decay rate
                ss0, ... % steady state without heterologous gene expression
                e0, ... % steady state translation elongation rate without heterologous gene expression
                sim) % simulator (required to get parameters and regulatory functions for all genes)

            % get growth rate
            l=obj.ss_l(plasmid_conc,ss0,e0,sim);

            %  get mRNA level
            m=F.*c.*a./(l+b);
        end
        
        % steady state protein mass fraction
        function phi=ss_phi(obj, ...
                plasmid_conc, ... % concentration of plasmid with heterologous genes
                F, ... % gene of interest (GOI) ss regulation function value
                c, ... % GOI concentration
                a, ... % GOI transcription rate
                kplus, ... % GOI mRNA-ribosome binding rate
                kminus, ... % GOI mRNA-ribosome unbinding rate
                n, ...  % GOI length in aa
                ss0, ... % steady state without heterologous gene expression
                e0, ... % steady state translation elongation rate without heterologous gene expression
                sim) % simulator (required to get parameters and regulatory functions for all genes)
            
            % store model parameters in arrays
            par_main=sim.parameters;
            par_het=sim.het.parameters;

            % GET the denominator for phi_i expression
            denominator=obj.denominator(plasmid_conc,ss0,e0,par_main,par_het,sim);

            % get the phi value
            k = sim.coll.k(e0,kplus,kminus,n,par_main('kcm').*par_main('h'));
            phi = F.*c.*a./k./denominator;
        end

        % steady state protein mass fraction OF ALL HETEROLOGOUS GENES
        function phi_het=ss_phi_het(obj, ...
                plasmid_conc, ... % concentration of plasmid with heterologous genes
                ss0, ... % steady state without heterologous gene expression
                e0, ... % steady state translation elongation rate without heterologous gene expression
                sim) % simulator (required to get parameters and regulatory functions for all genes)
            
            % store model parameters in arrays
            par_main=sim.parameters;
            par_het=sim.het.parameters;

            % GET the denominator for phi_i expression
            denominator=obj.denominator(plasmid_conc,ss0,e0,par_main,par_het,sim);

            % get numerator (sum of all heterologous genes' contributions)
            numerator=0; % initialise numerator
            for i=1:sim.het.num_genes
                k_i=sim.het.coll_main.k(e0,par_het(['k+_',sim.het.names{i}]),...
                    par_het(['k-_',sim.het.names{i}]),par_het(['n_',sim.het.names{i}]),par_main('kcm').*par_main('h'));
                F_i=sim.het.regulation(sim.het.names{i}); % right now, regulation assumed constant!
                numerator=numerator+F_i.*plasmid_conc.*par_het(['a_',sim.het.names{i}])./k_i;
            end

            % calculate the phi_het value
            phi_het=numerator./denominator;
        end

        % steady state growth rate
        function l=ss_l(obj, ...
                plasmid_conc, ... % concentration of plasmid with heterologous genes
                ss0, ... % steady state without heterologous gene expression
                e0, ... % steady state translation elongation rate without heterologous gene expression
                sim) % simulator (required to get parameters and regulatory functions for all genes)
            
            % store model parameters in arrays
            par_main=sim.parameters;
            par_het=sim.het.parameters;

            % get the l value
            F_r=sim.coll.F_r(par_main,ss0(5)./ss0(6));
            phi_r=obj.ss_phi(plasmid_conc,F_r,par_main('c_r'),par_main('a_r'),...
                par_main('k+_r'),par_main('k-_r'),par_main('n_r'),ss0,e0,sim);
            l=e0./par_main('n_r').*phi_r;
        end

        % steady state growth rate using FULL expression for l
        % defined through a quadratic equation A2*l^2+A1*l+A0=0
        function l=ss_l_full(obj, ...
                plasmid_conc, ... % concentration of plasmid with heterologous genes
                ss0, ... % steady state without heterologous gene expression
                e0, ... % steady state translation elongation rate without heterologous gene expression
                sim) % simulator (required to get parameters and regulatory functions for all genes)
            
            % store model parameters in arrays
            par_main=sim.parameters;
            par_het=sim.het.parameters;
            kcmh=par_main('kcm').*par_main('h');

            % get A2
            A2=1;

            % get A1
            A1=obj.denominator(plasmid_conc,ss0,e0,par_main,par_het,sim)+par_main('b_a');

            % get A0
            k_r=sim.coll.k(e0,par_main('k+_r'),par_main('k-_r'),par_main('n_r'),kcmh);
            F_r=sim.coll.F_r(par_main,ss0(5)./ss0(6));
            A0=e0./par_main('n_r').*(F_r.*par_main('c_r').*par_main('a_r')./k_r);
            
            % get l
            l=(-A1+sqrt(A1.^2+4*A0.*A2))./(2*A2);
        end
        
        % denominator used in the expressions for phi and l
        function denom=denominator(obj,plasmid_conc,ss0,e0,par_main,par_het,sim)
            denom = 0; % initialise denominator
            kcmh=par_main('kcm').*par_main('h');

            % add contribution of metabolic genes
            k_a=sim.coll.k(e0,par_main('k+_a'),par_main('k-_a'),par_main('n_a'),kcmh);
            denom=denom+par_main('c_a').*par_main('a_a')./k_a;

            % add contribution of ribosomal gene
            k_r=sim.coll.k(e0,par_main('k+_r'),par_main('k-_r'),par_main('n_r'),kcmh);
            F_r=sim.coll.F_r(par_main,ss0(5)./ss0(6));
            denom=denom+F_r.*par_main('c_r').*par_main('a_r')./k_r;

            % add contributions of heterologous genes
            for i=1:sim.het.num_genes
                k_i=sim.het.coll_main.k(e0,par_het(['k+_',sim.het.names{i}]),...
                    par_het(['k-_',sim.het.names{i}]),par_het(['n_',sim.het.names{i}]),kcmh);
                F_i=sim.het.regulation(sim.het.names{i}); % right now, regulation assumed constant!
                denom=denom+F_i.*plasmid_conc.*par_het(['a_',sim.het.names{i}])./k_i;
            end

            % rescale to account for housekeeping genes
            denom=denom./(1-sim.parameters('phi_q'));
        end
    end

end