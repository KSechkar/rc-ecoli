classdef pulse_inducer
    % describe genes and their parameters
    properties (SetAccess = public)
        module_name;
        names;
        compatible_hets;
        input_func_parameters;
    end

    % regulatory functions
    methods (Access = public)
        function obj =pulse_inducer(obj)
            obj.module_name='pulse_inducer';
            % -------------------------------------------------------------
            % SPECIFY EXTERNAL INPUT INFO----------------------------------
            
            % SPECIFY INPUT NAMES HERE
            obj.names={'inducer'};

            % SPECIFY COMPATIBLE HETEROLOGOUS MODULES
            % (e.g. input depends on flu. readout => only modules with GFP compatible)
            % leave empty if it does not matter
            obj.compatible_hets={};


            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------


            % default parameters
            obj.input_func_parameters=containers.Map('KeyType', 'char', ...
                    'ValueType', 'double');

            % -------------------------------------------------------------
            % SPECIFY INPUT INFO--------------------------------------------

            % SPECIFY NON-DEFAULT PARAMETERS HERE
            obj.input_func_parameters('inducer_base_level')=100; % constant conc of inducer(nM)
            obj.input_func_parameters('pulse_value_prop')=1.25; % specify value during the pulse relative to baseline
            obj.input_func_parameters('pulse_start_time')=10; % specify start time of the pulse
            obj.input_func_parameters('pulse_duration')=1; % pulse duration (h)

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------
        end
        
        % regulation functions for heterologous mRNA transcription
        % sensing state of the system
        function inp = input(obj,x,t)

            % -------------------------------------------------------------
            % SPECIFY EXTERNAL INPUT FUNCTION------------------------------
            ifpar=obj.input_func_parameters;

            if (t>ifpar('pulse_start_time') && t<ifpar('pulse_start_time')+ifpar('pulse_duration'))
                inp=ifpar('inducer_base_level').*ifpar('pulse_value_prop');
            else
                inp=ifpar('inducer_base_level');
            end

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------
        end

    end
end