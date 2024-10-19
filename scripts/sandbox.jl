## 2024-09-04 copy from 2024-06-29 BLAKJacScript_B1.jl, to check on compatibility
include("setup.jl")
recon_options = Dict() # erase all existing settings
nsweeps = 6                                               
nky = 224                                                 
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)
rfFunction = rfDictionary["from_file"]
saved_H = Dict()

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.01
recon_options["startstate"] = -1 
recon_options["sigma_ref"] = 1.4 # See logbook 20220815
# recon_options["optpars"]   = Optim.Options(time_limit = 20000.0, iterations = 100000, f_tol=1.0e-10, g_tol = 1.0e-10)
recon_options["optpars"]   = Optim.Options(time_limit = 20000.0, iterations = 100000, f_tol=1.0e-5, g_tol = 1.0e-5)  
# recon_options["opt_method"] = SimulatedAnnealing   
# recon_options["optpars"]   = Optim.Options(time_limit = 3000.0, iterations = 20000)  

recon_options["opt_criterion"] = "noise_level" 
recon_options["account_SAR"]   = true     

recon_options["sar_limit"] = 40^2/0.01 
recon_options["emphasize_low_freq"] = true 
recon_options["handleB1"] = "sensitivity" # "no" # "sensitivity" #
recon_options["lambda_B1"] = 10.0 # 100.0 # was 10.0     
recon_options["opt_initialize"] = "cRandom30" 
# recon_options["opt_initialize"] = "ernst" 
recon_options["opt_focus"] = "max"      
recon_options["opt_complex"] = false      
recon_options["opt_account_maxFlip"] = false
recon_options["opt_keep_positive"] = false                           
recon_options["opt_slow_phase"] =  true # false # true                         
recon_options["considerCyclic"] = false  
recon_options["opt_emergeCriterion"] = 2000 # 500 # 2000
ph = [] 
# ph = zeros(nTR); ph .= 2.0;
recon_options["opt_imposed_2nd_derivative_of_phase"] = ph
recon_options["opt_iterations_limit"] = 1
recon_options["sizeSteps"] = [6]  
recon_options["B1metric"] = "multi_point_values" # "multi_point"         
nRealizations = 2     

# In the original script: fn_base = "20240703X"
fn_base = "20240911A"
for i in 1:nRealizations; goodseed = i
    stageText = ""
    portionRange = 0:0
    fn = "/home/mfuderer/Documents/Julia/Capture/$fn_base($i).jld2"
    RFdeg = BLAKJac.BLAKJac_optimize(trajectorySet, recon_options, goodseed);
    FileIO.save(fn,"RFdeg",RFdeg)
end 


## 2024-10-16 comparisons
include("setup.jl")
recon_options = Dict() # erase all existing settings
nsweeps = 6                                               
nky = 224                                                 
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)
rfFunction = rfDictionary["from_file"]
saved_H = Dict()

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.01
recon_options["startstate"] = -1 
recon_options["sigma_ref"] = 1.4 # See logbook 20220815

recon_options["emphasize_low_freq"] = true 
recon_options["handleB1"] = "sensitivity" # "no" # "sensitivity" #                      
recon_options["considerCyclic"] = false  

recon_options["B1metric"] = "multi_point_values" # "multi_point"         

@show recon_options["T1T2set"]
include("BLAKJac_B1_figures.jl")
testCase = "20240703X"; Nr=10
map_B1_sensitivities(testCase,5)
@show recon_options["T1T2set"]






