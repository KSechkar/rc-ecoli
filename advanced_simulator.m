%% simulator.m
%%% OCTOBER 4, 2021

classdef advanced_simulator
    
    properties (SetAccess = public)
        x0;
        t;
        x;

        opt = odeset('RelTol', 1e-12, 'AbsTol', 1e-16);
        tf = 100;
        init_conditions;
        parameters;
        coll=advanced_collec; % collection of rate and activation functions

        het; % defining heterologous genes
    end
    
    methods (Access = public)
        
        function obj = advanced_simulator(tf) % Constructor
            if nargin == 1
                obj.tf = tf;
            end
            obj.parameters=advanced_params(); % set defualt parameters
            obj.init_conditions=advanced_init_conds(obj.parameters); % set init conditions

            obj.het=heterologous_genes(); % set up heterologous genes (default: one extra gene)
        end
        
        function obj = simulate_model(obj)
            obj = obj.set_x0; % set initial condition for main cell genes
            % set initial condition for heterologous genes
            obj.het = obj.het.set_x0();
            [t_full, x_full] = ode15s(@obj.ss_model, [0, obj.tf], [obj.x0; obj.het.x0], obj.opt);
            obj.x=x_full(:,1:7);
            obj.t=t_full;
            obj.het.x=x_full(:,8:end);
            obj.het.t=t_full;
        end
        
        function plot_simulation(obj)
            % Protein massess (aa)
            str_var = {'q','a','R','Bcm','het'};
            amino_q=obj.parameters('M')*ones(size(obj.t))*obj.parameters('phi_q');
            amino_a=obj.x(:,3)*obj.parameters('n_a');
            amino_r=obj.x(:,4)*obj.parameters('n_r');
            amino_Bcm=obj.x(:,7)*obj.parameters('n_r');
            amino_het=zeros(size(obj.t)); % initialise masses of heterologous proteins
            if(obj.het.num_genes~=0)
                for i=1:obj.het.num_genes
                    amino_het=sum((amino_het+obj.het.x(:,(obj.het.num_genes+i))...
                        *obj.het.parameters(['n_',obj.het.names{i}])),2);
                end
            end
            Faa = figure('Position',[0 0 540 480]);
            set(Faa, 'defaultLineLineWidth', 2)
            set(Faa, 'defaultAxesFontSize', 16)
            hold all
            patch([obj.t; flip(obj.t)],[amino_q+amino_a+amino_r+amino_Bcm+amino_het; flip(amino_a+amino_r+amino_Bcm+amino_het)], [0.75, 0.75, 0.75])
            patch([obj.t; flip(obj.t)],[amino_a+amino_r+amino_Bcm+amino_het; flip(amino_r+amino_Bcm)], [0.9290, 0.6940, 0.1250])
            patch([obj.t; flip(obj.t)],[amino_r+amino_Bcm+amino_het; flip(amino_Bcm+amino_het)], [0.4940, 0.1840, 0.5560])
            patch([obj.t; flip(obj.t)],[amino_Bcm+amino_het; amino_het], [0.4660 0.6740 0.1880])
            patch([obj.t; flip(obj.t)],[amino_het; zeros(size(obj.t))], [0 0.4470 0.7410])

            xlabel('Time (h)');
            ylabel('Protein mass (aa)');
            legend(str_var{[1,2,3,4,5]});
            title('Prot. mass (aa) evolution');
            % ylim([0 obj.parameters('M')])
            hold off
        end
        
    end
    
    methods (Access = protected)
        
        function obj = set_x0(obj)
            obj.x0 = [
                      % mRNAs;
                      obj.init_conditions('m_a');
                      obj.init_conditions('m_r');
                      % proteins
                      obj.init_conditions('p_a');
                      obj.init_conditions('R'); % non-inactivated ribosomes
                      % tRNAs
                      obj.init_conditions('tc'); % charged
                      obj.init_conditions('tu'); % free
                      % ribosomes inactivated by chloramphenicol
                      obj.init_conditions('Bcm');
                      ];
        end
        
        function dxdt = ss_model(obj, t, x)
            % PAREMETER MAP TO VARIABLE PAR
            par = obj.parameters;
            
            % STATE VECTOR TO SINGLE VARIABLES
            m_a = x(1);
            m_r = x(2);
            p_a = x(3);
            R = x(4);
            tc = x(5);
            tu = x(6);
            Bcm = x(7);

            % USEFUL PRE-CALCULATIONS
            % translation elongation rate
            e=obj.coll.e(par,tc);

            % ribosome inactivation rate due to chloramphenicol
            kcmh=par('kcm').*par('h');

            % ribosome dissociation constants
            k_a=obj.coll.k(e,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
            k_r=obj.coll.k(e,par('k+_r'),par('k-_r'),par('n_r'),kcmh);

            % consider the contribution of heterologous genes to resource competition if there are any
            mk_het=0;
            if(obj.het.num_genes>0)
                x_het=x(8:end); % current state of heterologous genes
                obj.het = obj.het.find_current_ks(e,kcmh); % ribosome dissociation constants for heterolgous genes
                for i=1:obj.het.num_genes
                    mk_het = mk_het + x_het(i)./obj.het.current_ks(i); % get sum of mx/kx for heterologous genes
                end
            end

            T=tc./tu; % ratio of charged to uncharged tRNAs
            D=1+(m_a./k_a+m_r./k_r+mk_het)./(1-par('phi_q')); % denominator in ribosome competition cobj.optalculations
            B=R.*(1-1./D); % actively translating ribosomes - INCLUDING Q

            % tRNA charging and synthesis
            nu=obj.coll.nu(par,tu);
            psi=obj.coll.psi(par,T);

            % growth rate
            l=obj.coll.l(par,e,B);
            
            % DEFINE DX/DT
            dxdt = [
                    % mRNAs
                    par('c_a').*par('a_a')-(par('b_a')+l).*m_a;
                    obj.coll.F_r(par,T).*par('c_r').*par('a_r')-(par('b_r')+l).*m_r;
                    % ,metabolic protein a
                    (e./par('n_a')).*(m_a./k_a./D).*R-l.*p_a;
                    % ribosomes
                    (e./par('n_r')).*(m_r./k_r./D).*R-l.*R-kcmh.*B;
                    % tRNAs
                    nu.*p_a-l.*tc-e.*B;
                    psi-l.*tu-nu.*p_a+e.*B;
                    % ribosomes inactivated by chloramphenicol
                    kcmh.*B-l.*Bcm;
                    ];
            % concatenate dxdt for heterologous genes if there are any
            if(obj.het.num_genes>0)
                dxdt=[dxdt;...
                    obj.het.dxdt(t,x,par,e,kcmh,D,l)];
            end

            % TEST ONLY
%             eB=(e).*(m_r./k_r./D).*R+...
%                 (e).*(m_a./k_a./D).*R+...
%                 (e).*(x_het(1)./obj.het.current_ks(1)./D).*R+...
%                 (e).*(x_het(2)./obj.het.current_ks(1)./D).*R;
%             eB=eB./0.45;
%             [eB/par('M'),l]
        %(R*par('n_r')+p_a*par('n_a')+x(10)*obj.het.parameters('n_out')+x(11)*obj.het.parameters('n_ctrl'))./0.45
        end
        
    end
    
end