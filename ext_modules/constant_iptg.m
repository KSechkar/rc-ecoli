classdef constant_iptg
    % describe genes and their parameters
    properties (SetAccess = public)
        module_name;
        names;
        compatible_hets;
        input_func_parameters;
    end

    % regulatory functions
    methods (Access = public)
        function obj =constant_iptg(obj)
            obj.module_name='constant_iptg';
            % -------------------------------------------------------------
            % SPECIFY EXTERNAL INPUT INFO----------------------------------
            
            % SPECIFY INPUT NAMES HERE
            obj.names={'iptg'};

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
            % SPECIFY GENE INFO--------------------------------------------

            % SPECIFY NON-DEFAULT PARAMETERS HERE
            obj.input_func_parameters('iptg_level')=100; % constant conc of IPTG (nM)

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------
        end
        
        % regulation functions for heterologous mRNA transcription
        % sensing state of the system
        function inp = input(obj,x,t)

            % -------------------------------------------------------------
            % SPECIFY EXTERNAL INPUT FUNCTION------------------------------

            inp=obj.input_func_parameters('iptg_level');

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------
        end

    end
end