%% simulator.m
%%% OCTOBER 4, 2021

classdef advanced_simulator_corrected_inactivation
    
    properties (SetAccess = public)
        x0;
        t;
        x;

        opt = odeset('RelTol', 1e-12, 'AbsTol', 1e-16);
        tf = 100;
        init_conditions;
        parameters; % main cell parameters
        coll=advanced_collec; % collection of rate and activation functions
        
        num_het=0; % number of heterologous genes
        num_misc=0;
        het; % object describing heterologous genes

        num_ext=0; % number of external inputs
        ext; % object describing external inputs
    end
    
    methods (Access = public)
        
        function obj = advanced_simulator_corrected_inactivation(tf) % Constructor
            if nargin == 1
                obj.tf = tf;
            end
            % set default parameters and initial conditions of the native system
            obj=obj.set_default_parameters();
            obj=obj.set_default_init_conditions();
            
            % set up heterologous genes and external inputs (SELECT APPROPRIATE CLASS FILES)
            obj=obj.load_heterologous_and_external('no_het','no_ext'); % none by default

            % push parameters and initial conditions of het. system into main framework
            obj=obj.push_het();
        end

        function obj=set_default_parameters(obj)
            obj.parameters=advanced_params();
        end

        function obj=set_default_init_conditions(obj)
            obj.init_conditions=advanced_init_conds(obj.parameters);
        end

        function obj=push_het(obj)
            obj=obj.push_het_parameters();
            obj=obj.push_het_init_conditions();
        end

        function obj = load_heterologous_and_external(obj,het_sys,ext_sys)
            % LOAD HETEROLOGOUS GENES MODULE
            % access the relevant file
            addpath(genpath([pwd, filesep,'het_modules']));
            het_func=str2func(het_sys);

            % class file describing a synthetic gene system
            obj.het=het_func();

            % get number of heterologous genes
            obj.num_het=size(obj.het.names,2);

            % get number of miscellaneous species
            obj.num_misc=size(obj.het.misc_names,2);

            % LOAD EXTERNAL INPUTS MODULE
            % access the relevant file
            addpath(genpath([pwd, filesep,'ext_modules']));
            ext_func=str2func(ext_sys);

            % class file describing a synthetic gene system
            obj.ext=ext_func();

            % get number of heterologous genes
            obj.num_ext=size(obj.ext.names,2);

            % push parameters
            obj=obj.push_het();

            % CHECK COMPAITIBILITY
            % is het gene module compatible with external input module?
            if(size(obj.ext.compatible_hets,2)~=0)
                het_in_compatible=false;
                for i=1:size(obj.ext.compatible_hets,2)
                    if strcmp(obj.het.module_name,obj.ext.compatible_hets{i})
                        het_in_compatible=true;
                        break;
                    end
                end
            else
                het_in_compatible=true;
            end

            % does external input module allow the het gene module to work?
            if(size(obj.het.prerequisite_exts,2)~=0)
                ext_in_prerequisite=false;
                for i=1:size(obj.het.prerequisite_exts,2)
                    if strcmp(obj.ext.module_name,obj.het.prerequisite_exts{i})
                        ext_in_prerequisite=true;
                        break;
                    end
                end
            else
                ext_in_prerequisite=true;
            end

            % check compaitibility
            if ~(het_in_compatible && ext_in_prerequisite)
                disp('Incompatible modules! Expect errors')
            end
        end

        function obj = push_het_parameters(obj)
            % add heterologous genes if there are any
            if(obj.num_het>0)
                for key=keys(obj.het.parameters)
                    obj.parameters(key{1})=obj.het.parameters(key{1});
                end
            end
        end

        function obj = push_het_init_conditions(obj)
            % add heterologous genes if there are any
            if(obj.num_het>0)
                for key=keys(obj.het.init_conditions)
                    obj.init_conditions(key{1})=obj.het.init_conditions(key{1});
                end
            end
        end
        
        function obj = simulate_model(obj)
            obj = obj.set_x0; % set initial condition
            [obj.t, obj.x] = ode15s(@obj.ss_model, [0, obj.tf], [obj.x0], obj.opt);
        end
        
        function obj = set_x0(obj)
            % NATIVE GENES
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
                      % (constant) nutrient quality and chloramphenicol
                      obj.init_conditions('s');
                      obj.init_conditions('h');
                      ];

            % ...ADD HETEROLOGOUS GENES IF THERE ARE ANY
            if(obj.num_het>0)
                x0_het=zeros(2*obj.num_het,1); % initialise
                % calculate
                for i=1:obj.num_het
                    % mRNA
                    x0_het(i)=obj.init_conditions(['m_',obj.het.names{i}]);
                    % protein
                    x0_het(i+obj.num_het)=obj.init_conditions(['p_',obj.het.names{i}]);
                end
                obj.x0=[obj.x0;x0_het]; % concantenate
            end
            % ...ADD MISCELLANEOUS SPECIES IF THERE ARE ANY
            if(obj.num_misc>0)
                x0_misc=zeros(obj.num_misc,1); % initialise
                % calculate
                for i=1:obj.num_misc
                    x0_misc(i)=obj.init_conditions(obj.het.misc_names{i});
                end
                obj.x0=[obj.x0;x0_misc]; % concantenate
            end
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
            s = x(8);
            h = x(9);
            x_het=x(10: (9+2*obj.num_het)); % heterologous genes 

            % USEFUL PRE-CALCULATIONS
            % translation elongation rate
            e=obj.coll.e(par,tc);

            % ribosome inactivation rate due to chloramphenicol
            kcmh=par('kcm').*h;

            % ribosome dissociation constants
            k_a=obj.coll.k(e,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
            k_r=obj.coll.k(e,par('k+_r'),par('k-_r'),par('n_r'),kcmh);
            % heterologous genes
            k_het=ones(obj.num_het,1); % initialise with default value 1
            if(obj.num_het>0)
                for i=1:obj.num_het
                    k_het(i)=obj.coll.k(e,...
                    obj.parameters(['k+_',obj.het.names{i}]),...
                    obj.parameters(['k-_',obj.het.names{i}]),...
                    obj.parameters(['n_',obj.het.names{i}]),...
                    kcmh);
                end
            end

            T=tc./tu; % ratio of charged to uncharged tRNAs
            D=1+(m_a./k_a+m_r./k_r+sum(x_het(1:obj.num_het)./k_het))./...
                (1-par('phi_q')); % denominator in ribosome competition cobj.optalculations
            B=R.*(1-1./D); % actively translating ribosomes - INCLUDING Q

            % tRNA charging and synthesis
            nu=obj.coll.nu(par,tu,s);
            psi=obj.coll.psi(par,T);

            % growth rate
            l=obj.coll.l(par,e,B);

            % GET EXTERNAL INPUT
            ext_inp=obj.ext.input(x,t);
            
            % DEFINE DX/DT...
            % ...NATIVE SYSTEM
            dxdt = [
                    % mRNAs
                    par('c_a').*par('a_a')-(par('b_a')+l).*m_a-kcmh.*(m_a./k_a./D).*R;
                    obj.coll.F_r(par,T).*par('c_r').*par('a_r')-(par('b_r')+l).*m_r-kcmh.*(m_r./k_r./D).*R;
                    % ,metabolic protein a
                    (e./par('n_a')).*(m_a./k_a./D).*R-l.*p_a;
                    % ribosomes
                    (e./par('n_r')).*(m_r./k_r./D).*R-l.*R-kcmh.*B;
                    % tRNAs
                    nu.*p_a-l.*tc-e.*B;
                    psi-l.*tu-nu.*p_a+e.*B;
                    % ribosomes inactivated by chloramphenicol
                    kcmh.*B-l.*Bcm;
                    % nutrient quality assumed constant
                    0;
                    % chloramphenicol concentration assumed constant
                    0;
                    ];

            % ...HETEROLOGOUS GENES
            if(obj.num_het>0)
                dxdt_het=zeros(2*obj.num_het,1); % initialise
                % calculate
                for i=1:obj.num_het
                    dxdt_het(i)=obj.het.regulation(obj.het.names{i},x,ext_inp)...
                        .*par(['c_',obj.het.names{i}]).*par(['a_',obj.het.names{i}])...
                        -(par(['b_',obj.het.names{i}])+l).*x_het(i)...
                        +obj.het.extra_m_term(obj.het.names{i},x,ext_inp);%...
                        %-kcmh.*(x_het(i)./k_het(i)./D).*R;

                    % protein
                    dxdt_het(i+obj.num_het)=(e./par(['n_',obj.het.names{i}])).*(x_het(i)./k_het(i)./D).*R...
                        -l.*x_het(i+obj.num_het)+...
                        obj.het.extra_p_term(obj.het.names{i},x,ext_inp);
                end
                dxdt=[dxdt;dxdt_het]; % concantenate
            end
            % ...MISCELLANEOUS SPECIES
            if(obj.num_misc>0)
                dxdt=[dxdt;obj.het.misc_ode(t,x,ext_inp,l)];
            end

        end

        function plot_simulation(obj,species,plot_type)
            if strcmp(species,'native')
                if strcmp(plot_type,'mass fractions')
                    % Protein massess (aa)
                    str_var = {'q','a','R','Bcm','het'};
                    amino_q=obj.parameters('M')*ones(size(obj.t))*obj.parameters('phi_q');
                    amino_a=obj.x(:,3)*obj.parameters('n_a');
                    amino_r=obj.x(:,4)*obj.parameters('n_r');
                    amino_Bcm=obj.x(:,7)*obj.parameters('n_r');
                    amino_het=zeros(size(obj.t)); % initialise masses of heterologous proteins
                    if(obj.num_het~=0)
                        x_het=obj.x(:,10: (9+2*obj.num_het) );
                        for i=1:obj.num_het
                            amino_het=sum((amino_het+x_het(:,(obj.num_het+i))...
                                *obj.parameters(['n_',obj.het.names{i}])),2);
                        end
                    end
                    Faa = figure('Position',[0 0 540 480]);
                    set(Faa, 'defaultLineLineWidth', 2)
                    set(Faa, 'defaultAxesFontSize', 16)
                    hold all
                    patch([obj.t; flip(obj.t)],[amino_q+amino_a+amino_r+amino_Bcm+amino_het; flip(amino_a+amino_r+amino_Bcm+amino_het)], [0.75, 0.75, 0.75])
                    patch([obj.t; flip(obj.t)],[amino_a+amino_r+amino_Bcm+amino_het; flip(amino_r+amino_Bcm)], [0.9290, 0.6940, 0.1250])
                    patch([obj.t; flip(obj.t)],[amino_r+amino_Bcm+amino_het; flip(amino_Bcm+amino_het)], [0.4940, 0.1840, 0.5560])
                    patch([obj.t; flip(obj.t)],[amino_Bcm+amino_het; flip(amino_het)], [0.4660 0.6740 0.1880])
                    patch([obj.t; flip(obj.t)],[amino_het; zeros(size(obj.t))], [0 0.4470 0.7410])
        
                    xlabel('Time (h)');
                    ylabel('Protein mass (aa)');
                    legend(str_var{[1,2,3,4,5]});
                    title('Prot. mass (aa) evolution');
                    % ylim([0 obj.parameters('M')])
                    hold off

                elseif strcmp(plot_type,'concentrations')
                    Fm = figure('Position',[0 0 1000 800]);
                    set(Fm, 'defaultLineLineWidth', 2)
                    set(Fm, 'defaultAxesFontSize', 14)
                    
                    % mRNA levels
                    subplot(2,2,1)
                    hold on
                    plot(obj.t,obj.x(:,1),'Color',[0.9290, 0.6940, 0.1250]);
                    plot(obj.t,obj.x(:,2),'Color',[0.4940, 0.1840, 0.5560]);
                    xlabel('Time (h)');
                    ylabel('mRNA concentration (nM)');
                    legend('a','r');
                    title('mRNA levels');
                    hold off
                    
                    % protein levels
                    subplot(2,2,2)
                    hold on
                    plot(obj.t,obj.x(:,3),'Color',[0.9290, 0.6940, 0.1250]);
                    plot(obj.t,obj.x(:,4),'Color',[0.4940, 0.1840, 0.5560]);
                    xlabel('Time (h)');
                    ylabel('mRNA concentration (nM)');
                    legend('a','r');
                    title('protein levels');
                    hold off
                    
                    % tRNA levels
                    subplot(2,2,3)
                    plot(obj.t,obj.x(:,5),obj.t,obj.x(:,6))
                    xlabel('Time (h)');
                    ylabel('tRNA concentration (nM)');
                    legend('charged','uncharged');
                    title('tRNA charging');
                end

            % HETEROLOGOUS PROTEINS
            elseif strcmp(species,'heterologous')
                colours=[[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.3010 0.7450 0.9330]; [0.6350 0.0780 0.1840]; [1 0 1]];% colours for plots
                if strcmp(plot_type,'mass fractions')
                    x_het=obj.x(:,10: (9+2*obj.num_het) );

                    Faa_het = figure('Position',[0 0 540 480]);
                    hold all
                    amino=zeros(size(obj.t,1),obj.num_het);

                    % record masses evolution
                    for i=1:size(obj.het.names,2)
                        amino(:,i)=obj.parameters(['n_',obj.het.names{i}])*x_het(:,obj.num_het+i);
                    end
                    % make plots
                    for i=1:obj.num_het
                        if(i~=obj.num_het)
                            patch([obj.t; flip(obj.t)],[sum(amino(:,i:end),2);...
                                flip(sum(amino(:,(i+1):end),2))], colours(rem(i-1,5)+1,:))
                        else
                            patch([obj.t; flip(obj.t)],[amino(:,end); zeros(size(obj.t))], colours(rem(i-1,5)+1,:))
                        end
                    end
                    % format plot
                    xlabel('Time (h)');
                    ylabel('Protein mass (aa)');
                    legend(obj.het.names);
                    title('Prot. mass (aa) evolution');
                    hold off
    
                elseif strcmp(plot_type,'concentrations') % plotting CONCENTRATIONS (aa)
                    Fc_het = figure('Position',[0 0 900 400]);
                    x_het=obj.x(:,10: (9+2*obj.num_het) );
                    
                    % protein levels
                    subplot(1,2,1)
                    hold on
                    for i=1:obj.num_het
                        plot(obj.t,x_het(:,(i+obj.num_het)),'Color',colours(rem(i-1,5)+1,:));
                    end
                    xlabel('Time (h)');
                    ylabel('protein concentration (nM)');
                    legend(obj.het.names);
                    title('Protein levels');
                    hold off
    
                    % mRNA levels
                    subplot(1,2,2)
                    hold on
                    for i=1:obj.num_het
                        plot(obj.t,x_het(:,i),'Color',colours(rem(i-1,5)+1,:));
                    end
                    xlabel('Time (h)');
                    ylabel('mRNA concentration (nM)');
                    legend(obj.het.names);
                    title('mRNA levels');
                    hold off

                elseif strcmp(plot_type,'regulation') % plotting REGULATION FUNCTION VALUES
                    F_reg = figure('Position',[0 0 500 400]);
                    % get external inputs
                    ext_inps=zeros(size(obj.t,1),1); % initialise
                    for i=1:size(obj.t,1)
                        ext_inps(i)=obj.ext.input(obj.x(i),obj.t(i));
                    end

                    % calculate regulation function values
                    F_het = zeros([size(obj.t,1) obj.num_het]);
                    for i=1:size(obj.t,1)
                        for j=1:obj.num_het
                            F_het(i,j) = obj.het.regulation(obj.het.names(j),obj.x(i,:),ext_inps(i));
                        end
                    end
                    
                    % plot
                    hold on
                    for i=1:obj.num_het
                        plot(obj.t,F_het(:,i),'Color',colours(rem(i-1,5)+1,:))
                    end
                    xlabel('Time (h)');
                    ylabel('F_x');
                    legend(obj.het.names);
                    title('Transcription regulation');
                    ylim([0 1])
                    hold off
                    
                end
             elseif strcmp(species,'miscellaneous')
                Fmisc = figure('Position',[0 0 500 400]);

                hold on
                for i=1:obj.num_misc
                    plot(obj.t,obj.x(:,9+obj.num_het*2+i))
                end
                xlabel('Time (h)');
                ylabel('Concentration, nM');
                legend(obj.het.misc_names);
                title('Miscellaneous species');
                hold off
            elseif strcmp(species,'external')
                Fext = figure('Position',[0 0 500 400]);

                if(obj.num_ext>0)
                    ext_inps=zeros(size(obj.t,1),1); % initialise
                    for i=1:size(obj.t,1)
                        ext_inps(i)=obj.ext.input(obj.x(i),obj.t(i));
                    end
                end
                
                plot(obj.t,ext_inps,'k');
                xlabel('Time (h)');
                ylabel('External input');
                legend(obj.ext.names);
                title('External input over time');
                hold off
            end
        end
    
        % return total heterologous mRNA transcription (conc.*transcription rate) for a given state of the system
        function tht = total_heterologous_transcription(obj, x, t)
            tht=zeros(1,size(x,1));
            if(obj.num_het>0)
                ext_inp=obj.ext.input(x,t);
                for i=1:obj.num_het
                    F_i=obj.het.regulation(obj.het.names(i),x,ext_inp);
                    tht=tht+F_i.*obj.parameters(['a_',obj.het.names{i}]).*obj.parameters(['c_',obj.het.names{i}]);
                end
            end
        end
    end
end