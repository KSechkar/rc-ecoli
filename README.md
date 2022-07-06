# rc-ecoli
MATLAB code for modelling how synthetic gene circuits interact with the host cell's (_E. coli_) native genes and affect its growth rate, used in Kirill Sechkar's MEng Molecular Bioengineering Master's thesis 'Modelling and Simulation of Biomolecular Controllers for Indirect Regulation of Gene Expression via Resource Competition Coupling'. The folders are organised as follows:

## cell_model
Contains scripts implementing the model of an E.coli cell. These include:
- _cell_simulator.m_ - a Matlab class object enabling simulations of the host cell. The expression of synthetic circuits can be simulated by loading the 'heterologous gene' and 'external input' modules (see _het_modules_ and _ext_inputs_ folders). Note that the associated modules' parameters must be pushed into the main simulator's memory using the function _push_het_ every time they altered.
- _cell_params.m_ - provides default values of all parameters decsribing the host cell
- - _cell_params.m_ - provides default initial conditions for the state of the cell
- _cell_formulae.m_ - contains a collection of formulae for rate and activation functions used by the cell simulator
- _get_steady.m_ - an auxiliary function that allows to determine the cell's steady state in given conditions by running _cell_simulator_ simulations.

## het_modules
Contains scripts defining class objects that allow to model the expression of different heterologous gene circuits. These include:
- _no_het.m_ - no heterologous gene being expressed. Due to being 'empty', this script can be copied and filled in to describe a synthetic gene circuit of interest.
- _one_constit.m_ - one constitutive heterologous gene.
- _two_constit.m_ - two constitutive heterologous genes.
- _aif_controller.m_ - antithetic integral feedback controller for maintaining a constant extent of competition for ribosomes in the cell (see Results 2.3.2 in the thesis)

## ext_modules
Contains scripts defining class objects that describe external inputs administered to the cells, such as chemical inducer concentrations or light stimuli. These include:
- _no_ext.m_ - no external signal. Due to being 'empty', this script can be copied and filled in to describe an externa input of interest.
- _constant_inducer.m_ - constant concentration of a chemical inducer.
- _pulse_inducer.m_ - time-limited step pulse above or below the baseline concentration of a chemical inducer

## generate_figures
Running the corresponding scripts allows to reproduce the simulation results presented in the thesis. Files not named _Figure(...).m_ or _AppFigure(...).m_ (for Appendix figures) are auxiliary scripts required by some of the figure-generating programs.

## param_fitting
Contains scripts allowing to fit the model's parameters to experimental data obtained by Scott et al. [^1] and processed by Chure et al. [^2] using the Delayed Rejection Adaptive Markoc Chain Monte Carlo (DRAM) algorithm [^3]. Note that these scripts require the [MCMCSTAT](https://github.com/mjlaine/mcmcstat) package to work.
- _dram_fit.m_ - run the DRAM algorithm and record the outcome in a _.mat_ file
- _dram_interpret.m_ - read a _.mat_ file describing a run of the algorithm and allow to interpret its results by plotting the MCMC chains, posterior distributions and model predictions with fitted parameters, as well as calculating the corresponding Fisher Information Matrix and parameter sensitivities
- _outcomes_ folder holds _.mat_ files containing the MCMC chain.
    -- _best_mcmc_outcome.mat_ - best outcome (i.e. provides the closest fit), which was used to determine the model's parameter values
    -- _mcmc_outcome_1.mat_, _mcmc_outcome_2.mat_, _mcmc_outcome_3.mat_ - other runs of the DRAM algorithm that yielded worse results

## data
Experimental data used to fit the parameters and compare model predictions with real-life measurements, taken from Chure et al.'s publication (2). The text file _annotation.txt_ explains the meaning of each dataset present.

---

REFERENCES
[^1] Matthew Scott, Carl W. Gunderson, Eduard M. Mateescu, Zhongge Zhang, and Terence Hwa. Interdependence of cell growth and gene expression: Origins and consequences. Science, 330(6007):1099–1102, 2010.

[^2] Griffin Chure and Jonas Cremer. An optimal regulation of fluxes dictates microbial growth in and out of steady-state. biorXiv, 2022.

[^3] Heikki Haario, Marko Laine, Antonietta Mira, and Eero Saksman. DRAM: Efficient adaptive MCMC. Statistics and Computing, 16(4):339–354, 2006.
