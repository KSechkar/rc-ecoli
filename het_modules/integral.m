% WARNING: MIND THE *100 FACTORS FOR C!

classdef integral
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
        function obj = integral(obj)
            obj.module_name='integral';

            % -------------------------------------------------------------
            % SPECIFY GENE INFO--------------------------------------------

            % SPECIFY PREREQUISITE EXTERNAL INPUT
            % the external input must be one of the listed
            % if no external inputs required, leave empty
            obj.prerequisite_exts={'constant_inducer','pulse_inducer'};

            % SPECIFY GENE NAMES HERE
            obj.names={'sens', ... % burden sensor, activates m_toe exp
                'toe', ... % toehold RNA, annihilates m_act
                'act', ... % actuator, affects cell burden
                'dist_i',... % disturbance: inducer protein
                'dist_o'... % disturbance: output protein
                };

            obj.misc_names={'bound toehold'};

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------


            % default parameters
            obj.parameters=containers.Map('KeyType', 'char', ...
                    'ValueType', 'double');
             for i=1:size(obj.names,2)
                obj.parameters(['c_',obj.names{i}]) = 100; % copy no. (nM) default: [1]*
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
            % Hill function for p_sens-DNA binding (arbitrary)
            obj.parameters('K_dna-sens')=1000;
            obj.parameters('eta_dna-sens')=1;

            % m_toe-m_act binding
            obj.parameters('kb_toe')=10; % (what would a realistic value be?)
            obj.parameters('b_bound_toehold')=6; % toehold degradation rate

            % m_toe NOT transcribed
            obj.parameters('k+_toe')=0;

            % max transcription rates
            obj.parameters('a_sens')=1;
            obj.parameters('a_toe')=50;
            obj.parameters('a_act')=25;

            % possibility to 'freeze the feedback'
            obj.parameters('freeze_fbk')=0; % 1 if true, 0 if not
            obj.parameters('fixed_F')=0.5; % fixed F_act value in case feedback is frozen

            % --- Disturbance ---

            % Input-p_dist_i and I-p_dist_i-DNA binding (Qian et al. 2017)
            obj.parameters('K_I-dist_i')=1000;
            obj.parameters('eta_I-dist_i')=1;
            obj.parameters('K_dna-Idist_i')=6000;
            obj.parameters('eta_dna-Idist_i')=2;
            
            % max transcription rates (arbitrary)
            obj.parameters('a_dist_i')=5;
            obj.parameters('a_dist_o')=100;

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------
        end
        
        % regulation functions for heterologous mRNA transcription
        % sensing state of the system
        function F = regulation(obj,gene_name,x,ext_inp)

            % -------------------------------------------------------------
            % SPECIFY REGULATORY FUNCTIONS---------------------------------
            het_par=obj.parameters;
            x_het=x( 10 : 9+size(obj.names,2)*2 ); % get heterologous gene info
            
            % regulating m_toe expression
            if strcmp(gene_name,'toe')
                p_sens=x_het(6);

                % Hill repression function
                F = het_par('K_dna-sens').^het_par('eta_dna-sens')./ ...
                    (het_par('K_dna-sens').^het_par('eta_dna-sens') + p_sens.^het_par('eta_dna-sens'));

                F = F.*(1-het_par('freeze_fbk')) + het_par('fixed_F').*het_par('freeze_fbk');

            elseif (strcmp(gene_name,'dist_o'))
                p_dist_i=x_het(9); % LuxR protein
                I=ext_inp(1); % inducer concentration (AHL)
                
                % I-dist_i complex: given by Hill function for Input-p_dist_i binding
                Idist_i = p_dist_i .* ...
                    I.^het_par('eta_I-dist_i')./(het_par('K_I-dist_i').^het_par('eta_I-dist_i') + I.^het_par('eta_I-dist_i'));
                
                % Hill function for p_dist_i-promoter binding
                F = Idist_i.^het_par('eta_dna-Idist_i')./ ...
                    (het_par('K_dna-Idist_i').^het_par('eta_dna-Idist_i') + Idist_i.^het_par('eta_dna-Idist_i'));              

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

            % extra term for m_act-m_toe binding
            if strcmp(gene_name,'act')
                m_toe=x_het(2);
                m_act=x_het(3);
                extra_term=-het_par('kb_toe').*m_toe.*m_act;

            % etra term for m_trig - binding m_comp
            elseif strcmp(gene_name,'toe')
                m_toe=x_het(2);
                m_act=x_het(3);
                extra_term=-het_par('kb_toe').*m_toe.*m_act;

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
            het_par=obj.parameters;
            x_het=x( 10 : 9+size(obj.names,2)*2 ); % get heterologous gene info
            m_toe=x_het(2);
            m_act=x_het(3);

            misc=x(10+size(obj.names,2)*2:end);
            toe=misc(1);

            dxdt=het_par('kb_toe').*m_toe.*m_act - (het_par('b_bound_toehold')+l).*toe;
        end
    end
end