function load_default_recon_options(recon_options = Dict())

    if (!@isdefined(recon_options)) 
        recon_options = Dict()
    end 
    BLAKJac_load_default_recon_options!(recon_options)
    return recon_options
end


function BLAKJac_load_default_recon_options!(recon_options)

    recon_options["kPattern"]  = 0           # function to generate (ky,kz)-pattern
    recon_options["rfName"]  = "from_optim" # name of the RF-shape function
    recon_options["rfFile"] = "RFoptimCFullRand_210530H" # name of the file to read the RF shape from
    recon_options["rflabel"] = ""           # For graph labels
    recon_options["nTR"]     = 1120         # number of simulated time points            
    recon_options["TR"]      = 0.01         # in seconds     
    recon_options["T1ref"]   = 0.67          # The reference T1 value around which linearity will be assumed
    recon_options["T2ref"]   = 0.076
    recon_options["startstate"] = 1.0       # Starting state of the z-magnetisation; +1 for no prepulse, -1 for inversion
    recon_options["halfrange_lt1"]=1.5      # The range over which T1-effect is calculated, spreads over exp(+-halfrange) around T1ref
    recon_options["halfrange_lt2"]=2.5      # The range over which T2-effect is calculated, spreads over exp(+-halfrange) around T2ref
    recon_options["maxstate"] = 40          # Length of history taken into account in simulating magnetisation from sequence
    recon_options["nky"]      = 224         # Number of different encoding values simulated; nTR/nky/nkz is number of samples per encoding    
    recon_options["nkz"]      = 1               
    recon_options["maxMeas"]  = 2000        # Maximum number of measurements per phase-encoding (ky) value
    recon_options["nPars"]    = 3           # 3 if (rho,T1,T2) 
                       #   !!!!!!!!!!!!!!!!!!! Note: nPars is disused (hard-coded) in DiagonalizedRun; the below is a somewhat clumsy surrogate 
    recon_options["handleB1"] = "no"        # "no", "sensitivity", "co-reconstruct"
    recon_options["considerCyclic"]=false   # If true, then the RF pattern (and only the RF pattern) is considered cyclic
    recon_options["useSymmetry"] = false     # Assume real-ness of rho, T1 and T2 and therefor symmetry in k-space
    recon_options["timeProbeset"]=[ 200 500 1000] # time points for which plots are generated
    recon_options["invregval"] = [200.0^(-2), 200.0^(-2), 200.0^(-2)] # The inverses of the regularization matrix
    recon_options["sigma_m"]  = 0.01        # Assumed inaccuracy level of magnetization estimates
    recon_options["fin_step"] = 0.7         # step size (log scale) to calculate higher orders via finite differences
    recon_options["Teighlor_probing_range"] = -4:4 # range of points from which Teighlor-expansion series are estimated
    recon_options["expansion_components"] = 1+2+3+4+5+6+7+8 # number of components (zeroth+linear+higherOrder); by pref, a triangular number
    recon_options["analysis_stepRange"] = -12:12 # Range for the maps to study the T1T2-dependency of (JʰJ)⁻¹
    recon_options["analysis_step"] = 0.2    # step*stepRange is the log(domain) around T1ref,T2ref that will be probed
    recon_options["useSurrogate"] = false   # Using the surrogate model tends to be beneficial if at least 24 (T1,T2)-combinations are to be evaluated
    recon_options["sigma_ref"] = 0.2        # Reference normalized noise level for information-content estimation. See logbook around 2021-06-03
    recon_options["sar_limit"] = 40^2/0.01  # SAR limit. Initially set to an RMS level of 40 degrees at TR=10ms

    r3 = sqrt(1/3.0)
    recon_options["T1T2Probeset"] = [(0.0,0.0), (0.0,1.0), (r3,0.5), (r3,-0.5), (0.0, -1.0), (-r3,-0.5), (-r3,0.5)] # deprecated !!!!
                                            # set of (T1,T2) combinations (in the log(T/Tref) space) for which n-factor and a-factor will be evaluated
    recon_options["T1T2set"]      = [(0.8183, 0.0509), (0.67, 0.2066), (1.2208, 0.1253), (0.2465, 0.0461), (2.2245, 0.3082), (0.3677, 0.0461), (0.3677, 0.1253)]
    recon_options["T1T2Plotaxis"] = (0.0,0.0)                                        
    recon_options["sizeSteps"] = [15 20 30 40 50]
                                            # In the first stage of the optimization procedure, the resolutions of the sequence that will be optimized upon

    recon_options["portion_size"] = 50      # In an optimization procedure, the size of the sequence-portion that will be optimized "in-context" 
    recon_options["portion_step"] = 30      # The portion-to-be-optimized is stepped over this number of TRs
                                            # Parameters for the NelderMead optimization algorithm
    recon_options["optpars"] = Optim.Options(time_limit = 3000.0, iterations = 100000, f_tol=1.0e-2, g_tol = 1.0e-3)
    recon_options["opt_method"] = NelderMead
    recon_options["opt_expand_type"] = "piecewise" # in ["piecewise", "spline", "nodes"]
    recon_options["opt_numberOfGaussians"] = 15 # for the DiOptimizeMultiGauss approach

    recon_options["opt_keep_positive"] = false;    # if true, then disallow negative RF values
    recon_options["optimization_set"] = ["RFreal"] # the set of variables to optimize on
    recon_options["opt_iterations_limit"] = 999 # cutting short the optimization-iteration process (while not yet confident on processing)
    recon_options["opt_initialize"] = "cRandom30" # type of initialization for optimizer, 
                                                  # in ["cRandom30", "ones", "ernst", "init_angle", "quadraticPhase30", "RF_shape"]
    recon_options["opt_complex"] = true      # complex optimization, each excitation having its own phase
    recon_options["opt_slow_phase"] = false # allow quadratic phase increment
    recon_options["opt_imposed_2nd_derivative_of_phase"] = [] # allows to provide a preselected 2nd derivative op the phase, in radiants
    recon_options["opt_focus"] = "mean"     # to focus on one of ["rho","T1","T2","mean","weighted","max"]
    recon_options["emphasize_low_freq"] = false # if switched on, then the optimization will steer more on reducing the noise factor on low spatial frequencies
    recon_options["opt_criterion"] = "information content" # "noise_level", "low-flip corrected noise", "information content", "information content penalized", "sar-limited noise"
    recon_options["opt_account_maxFlip"] = true # take into account that RF pulses need time
    recon_options["opt_emergeCriterion"] = 1000 # how seldomly do we want an in-between result?
    recon_options["opt_insertion_method"]="split" # {split, random} - for the multiGauss approach, the procedure to add one new Gaussian

    recon_options["plot"]     = true         # is plotting of intermediate results required?
    recon_options["plottypes"] = ["first"]   # collection of "first", "trajectories", "weights", "angles", "noisebars", "noiseSpectrum", "infocon"
    recon_options["T1T2testset"] = [(0.67,0.076), (3.2,0.019), (2.0,0.2), (0.1,0.076)]

    # The following are not really options but mis-used to store information over optimization-iterations 
    recon_options["optcount"]    = 0
    recon_options["stage_text"]  = ""


end