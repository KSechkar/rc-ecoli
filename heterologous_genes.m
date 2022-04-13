% Class describing all heterologous genes in the cell
classdef heterologous_genes
    
    properties (SetAccess = public)
        x0;
        x;
        t;

        coll_main=advanced_collec; % heterologous protein espression may also involve functions form the basic collection
        
        init_conditions; % intitial conditions for heterologous genes
        parameters; % parameters for heterologous gene expression

        num_genes; % number of heterologous genes
        names; % names of heterologous gene
        current_ks; % current ribosome binding affinities
    end
    
    methods (Access = public)
        % Constructor (EDIT TO SPECIFY GENES)
        function obj = heterologous_genes(obj)
            % ENTER GENE NAMES HERE
            obj.names = {'xtra'}; % default: one extra gene

            % set parameters and initial conditions
            obj.num_genes=size(obj.names,2);
            obj = obj.set_default_parameters(); % set defualt parameters
            obj = obj.set_default_init_conditions(); % set init conditions
            obj.current_ks=zeros(obj.num_genes,1); % initialise ribosome dissociation constants
        end
        
        % reload the heterologous gene simulator for a given list of names
        function obj = reset_for_names(obj,names)
            obj.names = names;

            % set parameters and initial conditions
            obj.num_genes=size(obj.names,2);
            obj = obj.set_default_parameters(); % set defualt parameters
            obj = obj.set_default_init_conditions(); % set init conditions
            obj.current_ks=zeros(obj.num_genes,1); % initialise ribosome dissociation constants
        end

        % Setting default parameter values (EDIT TO SPECIFY GENES)
        function obj = set_default_parameters(obj)
            obj.parameters = containers.Map('KeyType', 'char', ...
                'ValueType', 'double');
            for i=1:obj.num_genes
                % metabolic/aminoacylating genes
                obj.parameters(['c_',obj.names{i}]) = 1; % copy no. (nM) default: [1]*
                obj.parameters(['a_',obj.names{i}]) = 10; % max transcription rate (nM)
                obj.parameters(['b_',obj.names{i}]) = 6; % mRNA decay rate (/h) default: [1]
                obj.parameters(['k+_',obj.names{i}]) = 60; % ribosome binding rate (/h/nM) default: [1]
                obj.parameters(['k-_',obj.names{i}]) = 60; % ribosome unbinding rate (/h) default: [1]
                obj.parameters(['n_',obj.names{i}]) = 300; % protein length (aa) default: [1]
            end
            
            % ENTER YOUR PARAMETERS HERE IF NEEDED (inc. extras, such as Hill constants fo regulation)
        end
        
        % Setting default initial conditions (EDIT TO SPECIFY GENES)
        function obj = set_default_init_conditions(obj)
            obj.init_conditions = containers.Map('KeyType', 'char', ...
            'ValueType', 'double');
            for i=1:obj.num_genes
                obj.init_conditions(['m_',obj.names{i}])=0;
                obj.init_conditions(['p_',obj.names{i}])=0;
            end

            % RE-ENTER YOUR INITIAL CONDITIONS HERE IF NEEDED
        end
        
        % Regulatory functions for proteins (EDIT TO SPECIFY GENES)
        function F = regulation(obj,gene_name,varargin)
            % EDIT ELSEIFs WITH GENE NAMES (you may pass extra arguments using varargin)
            if strcmp(gene_name,'out')
                F=1;
            elseif strcmp(gene_name,'ctrl')
                F=1;
            else
                F=1; % constitutive expression by default
            end
        end
        
        function obj = find_current_ks(obj,e,kcmh)
            for i=1:obj.num_genes
                obj.current_ks(i)=obj.coll_main.k(e,...
                    obj.parameters(['k+_',obj.names{i}]),...
                    obj.parameters(['k-_',obj.names{i}]),...
                    obj.parameters(['n_',obj.names{i}]),...
                    kcmh);
            end
        end

        % Plotting simulation outcomes
        function plot_simulation(obj,which)
            colours=[[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.3010 0.7450 0.9330]; [0.6350 0.0780 0.1840]; [1 0 1]];% colours for plots
            
            if strcmp(which,'masses') % plotting PROTEIN MASSES (aa)
                Faa_het = figure('Position',[0 0 540 480]);
                hold all
                amino=zeros(size(obj.x(1)));
                % record masses evolution
                for i=size(obj.names,2)
                    amino(i)=obj.x(:,obj.num_genes+i)*obj.parameters(['n_',obj.names{i}]);
                end
                % make plots
                for i=1:obj.num_genes
                    if(i~=obj.num_genes)
                        patch([obj.t; flip(obj.t)],[sum(obj.x(:,(obj.num_genes+1):(obj.num_genes+i)),2);...
                            flip(sum(obj.x(:,(obj.num_genes+1):(obj.num_genes+i-1)),2))], colours(rem(i,5)))
                    else
                        patch([obj.t; flip(obj.t)],[sum(obj.x(:,2*obj.num_genes)); zeros(size(obj.t))], colours(rem(i,5)))
                    end
                end
                % format plot
                xlabel('Time (h)');
                ylabel('Protein mass (aa)');
                legend(obj.names);
                title('Prot. mass (aa) evolution');
                hold off

            elseif strcmp(which,'concentrations') % plotting CONCENTRATIONS (aa)
                Fc_het = figure('Position',[0 0 1000 400]);
                
                % protein levels
                subplot(1,2,1)
                hold on
                for i=1:obj.num_genes
                    plot(sim.t,sim.x(:,(i+obj.num_genes)),'Color',colours(rem(i,5)));
                end
                xlabel('Time (h)');
                ylabel('protein concentration (nM)');
                legend(obj.names);
                title('Protein levels');
                hold off

                % mRNA levels
                subplot(1,2,2)
                hold on
                for i=1:obj.num_genes
                    plot(sim.t,sim.x(:,i),'Color',colours(rem(i,5)));
                end
                xlabel('Time (h)');
                ylabel('mRNA concentration (nM)');
                legend(obj.names);
                title('mRNA levels');
                hold off
            end
        end
        
        % setting x0
        function obj = set_x0(obj)
            obj.x0 = zeros(2*obj.num_genes,1); % initialise
            % mRNAs
            for i=1:obj.num_genes
                obj.x0(i)=obj.init_conditions(['m_',obj.names{i}]);
            end
            % proteins
            for i=1:obj.num_genes
                obj.x0(obj.num_genes+i)=obj.init_conditions(['p_',obj.names{i}]);
            end
        end

        % returning x0
        function x0_het = get_x0(obj)
            x0_het=obj.x0;
        end
        
        % return dxdt for heterologous proteins
        function dxdt = dxdt(obj, ~,...
                x_main,... % cell's native proteins/mRNAs/tRNAs
                par_main,... %cell-wide parameters
                e,... % translation elongation rate
                kcmh,... % chloramphenicol conc. and binding rate
                D,... % resource competition denominator
                l... % growth rate
                )
            % GET HETEROLOGOUS PROTEIN AND MRNA LEVELS
            x_het=x_main(8:end);

            % PAREMETERS MAPPED TO VARIABLE PAR
            par_het = obj.parameters; % parameters specific for heterologous protein expression
            
            % DEFINE DX/DT
            dxdt = zeros(2*obj.num_genes,1); %initialise
            
            for i=1:obj.num_genes
                % mRNA
                dxdt(i)=obj.regulation(obj.names{i}).*...
                    par_het(['c_',obj.names{i}]).*par_het(['a_',obj.names{i}])...
                    -(par_het(['b_',obj.names{i}])+l).*x_het(i);
                % protein
                dxdt(i+obj.num_genes)=(e./par_het(['n_',obj.names{i}])).*(x_het(i)./obj.current_ks(i)./D).*x_main(4)...
                    -l.*x_het(i+obj.num_genes);
            end
        end
        
    end
    
end