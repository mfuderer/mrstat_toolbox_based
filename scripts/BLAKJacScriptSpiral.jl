## Intended for generation of cartesian-spiral 3D sequences as run by Oscar
#  Script code initially copied from BLAKJacScript_B1_1.jl, with main differences: 
#  no B1 optimization, no use of symmetry, no enhancement around the persumed ky=0   
#      (the latter is because it will be BLAKJac-ed as being a regular 2D sequence, 
#       but the factual sequence has ky=kz=0 at the beginning or end of each 'k-space'  )

## 2024-10-19 First try of the above 
include("setup.jl")
recon_options = Dict() # erase all existing settings
nsweeps = 5                               # not 6                                       
nky = 112                                 # not 224         
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)
rfFunction = rfDictionary["from_file"]
saved_H = Dict()

recon_options["useSymmetry"] = false      # differs from most applications     
recon_options["TR"]      = 0.008
recon_options["startstate"] = 1 # -1      # no inversion, contiguous 3D
recon_options["sigma_ref"] = 1.4 # See logbook 20220815
# recon_options["optpars"]   = Optim.Options(time_limit = 20000.0, iterations = 100000, f_tol=1.0e-10, g_tol = 1.0e-10)
recon_options["optpars"]   = Optim.Options(time_limit = 20000.0, iterations = 100000, f_tol=1.0e-5, g_tol = 1.0e-5)  
# recon_options["opt_method"] = SimulatedAnnealing   
# recon_options["optpars"]   = Optim.Options(time_limit = 3000.0, iterations = 20000)  

recon_options["opt_criterion"] = "noise_level" 
recon_options["account_SAR"]   = true     

recon_options["sar_limit"] = 40^2/0.01 
recon_options["emphasize_low_freq"] = false      # (see initial comments) 
recon_options["handleB1"] = "no"                 #               
recon_options["lambda_CSF"] = 5.0 # 0.0       
recon_options["opt_initialize"] = "cRandom30" 
# recon_options["opt_initialize"] = "ernst" 
recon_options["opt_focus"] = "max"      
recon_options["opt_complex"] = false      
recon_options["opt_account_maxFlip"] = false
recon_options["opt_keep_positive"] = false                           
recon_options["opt_slow_phase"] =  true # false # true                         
recon_options["considerCyclic"] = true # false  
recon_options["opt_emergeCriterion"] = 2000 # 500 # 2000
ph = [] 
# ph = zeros(nTR); ph .= 2.0;
recon_options["opt_imposed_2nd_derivative_of_phase"] = ph
recon_options["opt_iterations_limit"] = 1
recon_options["sizeSteps"] = [5]  
recon_options["B1metric"] = "multi_point_values" # "multi_point"         
nRealizations = 5      

fn_base = "20241019B"
for i in 1:nRealizations; goodseed = i
#for (i,goodseed) in enumerate([5])
    stageText = ""
    portionRange = 0:0
    fn = "/home/mfuderer/Documents/Julia/Capture/$fn_base($i).jld2"
    RFdeg = BLAKJac.BLAKJac_optimize(trajectorySet, recon_options, goodseed);
    FileIO.save(fn,"RFdeg",RFdeg)
end 

##

@show recon_options["T1T2set"]
include("BLAKJac_B1_figures.jl")
testCase = "20241019B"; Nr=5
map_B1_sensitivities(testCase,5)
@show recon_options["T1T2set"]

## 2024-10-19 text file output
include("setup.jl")
recon_options = Dict() # erase all existing settings
              
caselist     = [("20241019B(5)", "Cspiral560_noInv", 11, 112)]

for case in caselist
    nsweeps = 5
    (rfFile,name,tag,nky) = case
    nTR = round(Int64, nsweeps*nky)
    ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
    kz = ones(nsweeps*nky)
    trajectorySet = BLAKJac.TrajectorySet(ky,kz)
    BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)

    recon_options["rfName"]  = "from_file"
    recon_options["rfFile"]  = rfFile; 
    recon_options["rfFunction"] = rfDictionary[recon_options["rfName"]]
    rfFunction = recon_options["rfFunction"]
    RFdeg = rfFunction(recon_options["nTR"], recon_options["nky"])

    mkdir("/home/mfuderer/mrstat/data/txtRFfiles/"*name)
    afn = "/home/mfuderer/mrstat/data/txtRFfiles/"*name*"/rf_angle_list$tag.txt"
    pfn = "/home/mfuderer/mrstat/data/txtRFfiles/"*name*"/rf_phase_list$tag.txt"
    open(afn,"w") do io
        for f in RFdeg
            s=@sprintf("%7.3f \n",abs(f))
            write(io, s)
        end
    end
    open(pfn,"w") do io
        for f in RFdeg
            s=@sprintf("%8.3f \n", rad2deg(angle(f)))
            write(io, s)
        end
    end
end

