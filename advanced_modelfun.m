function ymodel=advanced_modelfun(theta,xdata,sim,Delta,Max_iter)
    % reset parameters and initial condition
    sim.parameters=advanced_params();
    sim.parameters('using_nutr_qual')=1;
    sim.init_conditions=advanced_init_conds(sim.parameters);

    % change fitted parameters to current values
    sim.parameters('a_a') = theta(1); % metabolic prot. transcription rate (/h)
    sim.parameters('a_r') = theta(2); % ribosome transcription rate (/h)
    sim.parameters('nu_max') = theta(3); % max metabolic rate (/h)
    sim.parameters('psi_max') = theta(4); % max tRNA synthesis rate (aa/h)  
    sim.parameters('K_t') = theta(5); % elongation rate and tRNA charging rate Hill constant (nM)
    sim.parameters('tau') = theta(6); % ppGpp sensitivity (ribosome transc. and tRNA synth. Hill const)
    sim.parameters('kcm') = theta(7); % chloramphenical binding rate constant (/h/nM)
    
    % estimate steady state values for every input
    ymodel=zeros([size(xdata,1) 2]); % initialise array of avlues predicted by the model
    for i=1:size(xdata,1)
        % change parameter values to relevant inputs
        sim.parameters('nutr_qual')=xdata(i,1);
        sim.parameters('h')=xdata(i,2);
        [l, R, Bm]=get_steady(sim,Delta,Max_iter); % evaluate steady state value

        ymodel(i,1)=l; % record growth rate
        ymodel(i,2)=(R+Bm).*sim.parameters('n_r')./sim.parameters('M'); % record ribosome mass fraction
    end
end