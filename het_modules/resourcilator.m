classdef resourcilator
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
        function obj =resourcilator(obj)
            obj.module_name='resourcilator';

            % -------------------------------------------------------------
            % SPECIFY GENE INFO--------------------------------------------

            % SPECIFY PREREQUISITE EXTERNAL INPUT
            % the external input must be one of the listed
            % if no external inputs required, leave empty
            obj.prerequisite_exts={'no_ext'};

            % SPECIFY GENE NAMES HERE
            obj.names={'rep',... % protein repressing miRNA expression
                'trig',... % RNA that triggers the toehold repressor on competitor's mRNA
                'comp' % gene competing with repressor for resources
                };
            obj.misc_names={'toe'}; % repressed toehold complex with competitor's mRNA

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
            % non-default transcription rates
            obj.parameters('a_rep')=0.5;
            obj.parameters('a_trig')=75;
            obj.parameters('a_comp')=50;

            % trig is a small mRNA - NOT TRANSLATED!
            % achieve that with a zero binding rate
            obj.parameters('k+_trig')=0;

            % Hill function for repression of trig
            obj.parameters('K_trig')=3000;
            obj.parameters('eta_trig')=2;

            % Toehold binding rate constant 
            obj.parameters('kb_toe')=10; % for the time being, arbitrary (??)

            % Degradation rate of bound toehold complex
            obj.parameters('b_toe')=6; % for the time being, same as all RNA species

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------
        end
        
        % regulation functions for heterologous mRNA transcription
        % sensing state of the system
        function F = regulation(obj,gene_name,x,ext_inp)

            % -------------------------------------------------------------
            % SPECIFY REGULATORY FUNCTIONS---------------------------------

            if strcmp(gene_name,'trig')
                het_par=obj.parameters;
                x_het=x( 10 : 9+size(obj.names,2)*2 ); % get heterologous gene info
                p_rep=x_het(4); % repressor protein in induction system         
                
                F=het_par('K_trig').^het_par('eta_trig')./...
                    (het_par('K_trig').^het_par('eta_trig') + p_rep.^het_par('eta_trig')); % Hill repression function

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
            het_par=obj.parameters;
            x_het=x( 10 : 9+size(obj.names,2)*2 ); % get heterologous gene info
            m_trig=x_het(2);
            m_comp=x_het(3);

            % extra term for m_comp - binding m_trig
            if strcmp(gene_name,'comp')
                extra_term=-het_par('kb_toe').*m_trig.*m_comp;

            % etra term for m_trig - binding m_comp
            elseif strcmp(gene_name,'trig')
                extra_term=-het_par('kb_toe').*m_trig.*m_comp;

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
            % m_trig and m_comp bound in a repressed toehold
            het_par=obj.parameters;
            x_het=x( 10 : 9+size(obj.names,2)*2 ); % get heterologous gene info
            m_trig=x_het(2);
            m_comp=x_het(3);

            misc=x(9+size(obj.names,2)*2+1:end);
            toe=misc(1);

            dxdt=het_par('kb_toe').*m_trig.*m_comp - (het_par('b_toe')+l).*toe;
        end
        
    end
end