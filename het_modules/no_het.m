classdef no_het
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
        function obj = no_het(obj)
            obj.module_name='no_het';

            % -------------------------------------------------------------
            % SPECIFY GENE INFO--------------------------------------------

            % SPECIFY PREREQUISITE EXTERNAL INPUT
            % the external input must be one of the listed
            % if no external inputs required, leave empty
            obj.prerequisite_exts={};

            % SPECIFY GENE NAMES HERE
            obj.names={};

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

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------
        end
        
        % regulation functions for heterologous mRNA transcription
        % sensing state of the system
        function F = regulation(obj,gene_name,x,ext_inp)

            % -------------------------------------------------------------
            % SPECIFY REGULATORY FUNCTIONS---------------------------------

            if strcmp(gene_name,'xtra')
                F=1;

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

            if strcmp(gene_name,'xtra')
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

            if strcmp(gene_name,'xtra')
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