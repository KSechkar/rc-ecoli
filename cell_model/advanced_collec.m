%% collec.m
% collection of rate and activation functions

classdef advanced_collec
   
    methods (Access = public)
        function k=k(obj,epsilon,kplus,kminus,n,kcmh) % risbome dissociation constant
            k=(kminus+epsilon./n+kcmh)./kplus;
        end

        function l = l(obj,par,epsilon,B) % growth rate
            l=epsilon.*B./par('M');
        end

        function F_q = F_q(obj,par,p_q) % housekeeping gene activation
            F_q=(par('K_q').^par('h_q'))./(p_q.^par('h_q')+par('K_q').^par('h_q'));
            %F_q=0;
        end

        function F_r = F_r(obj,par,T) % ribosome gene activation
            F_r=T./(T+par('tau')).*(1-par('is_fixed_F_r'))+par('fixed_F_r').*par('is_fixed_F_r');
        end
        
        function e = e(obj,par,tc) % elongation rate (aa/h)
            e = par('e_max')*tc./(tc+par('K_e'));
        end

        function nu = nu(obj,par,tu,s) % tRNA charging rate
            nu = par('nu_max').*(tu./(tu+par('K_nut'))).*s; % OR use a Hill function of nutrient conc. instead of s
        end

        function psi = psi(obj,par,T) % tRNA synthesis rate
            psi = par('psi_max').*T./(T+par('tau'));
        end
        
    end

end
