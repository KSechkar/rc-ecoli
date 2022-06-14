classdef induction
    % describe genes and their parameters
    properties (SetAccess = public)
        module_name;
        names;
        misc_names;
        parameters;
        init_conditions;

        prerequisite_exts;
    end

    % regulatory functions
    methods (Access = public)
        function obj =induction(obj)
            obj.module_name='induction';

            % -------------------------------------------------------------
            % SPECIFY GENE INFO--------------------------------------------

            % SPECIFY PREREQUISITE EXTERNAL INPUT
            % the external input must be one of the listed
            % if no external inputs required, leave empty
            obj.prerequisite_exts={'constant_inducer','pulse_inducer'};

            % SPECIFY GENE NAMES HERE
            obj.names={'nir',... % repressor in the negative-inducible system
                'out'... % output protein
                };

            obj.misc_names={};

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------


            % default parameters
            obj.parameters=containers.Map('KeyType', 'char', ...
                    'ValueType', 'double');
             for i=1:size(obj.names,2)
                obj.parameters(['c_',obj.names{i}]) = 1; % copy no. (nM) default: [1]*
                obj.parameters(['a_',obj.names{i}]) = 10; % max transcription rate (nM)
                obj.parameters(['b_',obj.names{i}]) = 6; % mRNA decay rate (/h) default: [1]
                obj.parameters(['k+_',obj.names{i}]) = 60; % ribosome binding rate (/h/nM) default: [1]
                obj.parameters(['k-_',obj.names{i}]) = 60; % ribosome unbinding rate (/h) default: [1]
                obj.parameters(['n_',obj.names{i}]) = 300; % protein length (aa) default: [1]
             end

            % default initial conditions
            obj.init_conditions=containers.Map('KeyType', 'char', ...
                    'ValueType', 'double');
             for i=1:size(obj.names,2)
                obj.init_conditions(['m_',obj.names{i}]) = 0;
                obj.init_conditions(['p_',obj.names{i}]) = 0;
             end
             for i=1:size(obj.misc_names,2)
                obj.init_conditions(obj.misc_names{i}) = 0;
             end


            % -------------------------------------------------------------
            % SPECIFY GENE INFO--------------------------------------------

            % SPECIFY NON-DEFAULT PARAMETERS HERE
            % Hill function for Input-p_nir binding
            obj.parameters('K_niri')=1000;
            obj.parameters('eta_niri')=1;

            % Hill function for p_nir-promoter binding
            obj.parameters('K_ind')=800;
            obj.parameters('eta_ind')=2;
            
            % max transcription rates
            obj.parameters('a_nir')=5;
            obj.parameters('a_out')=50;

            % initial condition (high p_nir levels)
            obj.init_conditions('p_nir')=0;

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------
        end
        
        % regulation functions for heterologous mRNA transcription
        % sensing state of the system
        function F = regulation(obj,gene_name,x,ext_inp)

            % -------------------------------------------------------------
            % SPECIFY REGULATORY FUNCTIONS---------------------------------

            if strcmp(gene_name,'out')
                het_par=obj.parameters;
                x_het=x( 10 : 9+size(obj.names,2)*2 ); % get heterologous gene info
                p_nir=x_het(3); % repressor protein in induction system
                I=ext_inp(1); % inducer concentration
                
                % Hill function for Input-p_nir binding
                Hill_niri=het_par('K_niri').^het_par('eta_niri')./ ...
                    (het_par('K_niri').^het_par('eta_niri') + I.^het_par('eta_niri'));
                
                % Hill function for p_nir-promoter binding
                F=het_par('K_ind').^het_par('eta_ind')./ ...
                    (het_par('K_ind').^het_par('eta_ind') + (p_nir.*Hill_niri).^het_par('eta_ind'));

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------
            else
                F=1; % constitutive expression by default
            end
        end

        % any extra terms in mRNA ODEs apart from synthesis & growth dilution
        function extra_term=extra_m_term(obj,gene_name,x,ext_inp)

            % -------------------------------------------------------------
            % SPECIFY TERMS---------------------------------

            if strcmp(gene_name,'out')
                extra_term=0;

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------

            else
                extra_term=0; % none by default
            end
        end

        % any extra terms in protein ODEs apart from synthesis & growth dilution
        function extra_term=extra_p_term(obj,gene_name,x,ext_inp)

            % -------------------------------------------------------------
            % SPECIFY REGULATORY FUNCTIONS---------------------------------

            if strcmp(gene_name,'out')
                extra_term=0;

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------
            
            else
                extra_term=0; % none by default
            end
        end
        
        % ODEs for miscellaneous species in the system
        function dxdt=misc_ode(obj,t,x,ext_inp,l)
            dxdt=[];
        end
    end
end
