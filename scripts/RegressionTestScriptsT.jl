## 2024-10-30 Intended for a series of regression-test scripts, intended to compare the output of MRSTAT-embedded
## BLAKJac with the output of the split-off BLAKJac repository (this instance)

## 2024-10-24 Start simple
include("setup.jl")
recon_options = Dict() # erase all existing settings
nsweeps = 6                                               
nky = 64 # keep it short, for processing time                                                 
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)
rfFunction = rfDictionary["from_file"]
saved_H = Dict()

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.008
recon_options["startstate"] = 1 # -1 
recon_options["sigma_ref"] = 1.4 # See logbook 20220815
recon_options["maxstate"] = 128
recon_options["optpars"]   = Optim.Options(time_limit = 200.0, iterations = 100000, f_tol=1.0e-5, g_tol = 1.0e-5)  ## time_limit!!!
recon_options["opt_criterion"] = "noise_level" 
recon_options["account_SAR"]   = false     
recon_options["emphasize_low_freq"] = true 
recon_options["handleB1"] = "no"       
recon_options["opt_initialize"] = "cRandom30" 
recon_options["opt_focus"] = "max"      
recon_options["opt_complex"] = false      
recon_options["opt_account_maxFlip"] = false
recon_options["opt_keep_positive"] = false                           
recon_options["opt_slow_phase"] =  false # true                         
recon_options["considerCyclic"] = false # true # false  
recon_options["opt_emergeCriterion"] = 500 # 2000
ph = [] 
# ph = zeros(nTR); ph .= 2.0;
recon_options["opt_imposed_2nd_derivative_of_phase"] = ph
recon_options["opt_iterations_limit"] = 1
recon_options["sizeSteps"] = [6]  
nRealizations = 2      

fn_base = "20241101G"
for i in 1:nRealizations; goodseed = i
#for (i,goodseed) in enumerate([5])
    stageText = ""
    portionRange = 0:0
    fn = "/home/mfuderer/Documents/Julia/Capture/$fn_base($i).jld2"
    RFdeg = BLAKJac.BLAKJac_optimize(trajectorySet, recon_options, goodseed);
    FileIO.save(fn,"RFdeg",RFdeg)
    noises, ItotAll, b1f = BLAKJac.BLAKJac_analysis!(cpu, RFdeg.+0im, trajectorySet, recon_options, saved_H)
    @show noises;
end 


## 2024-10-30 cross-compare  
include("setup.jl")
recon_options = Dict() # erase all existing settings
nsweeps = 6                                               
nky = 64 # keep it short, for processing time                                                 
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)
rfFunction = rfDictionary["from_file"]
saved_H = Dict()

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.008
recon_options["startstate"] = 1 # -1 
recon_options["sigma_ref"] = 1.4 # See logbook 20220815
recon_options["account_SAR"]   = false     
recon_options["emphasize_low_freq"] = true 
recon_options["handleB1"] = "no"                            
recon_options["considerCyclic"] = false # true # false 

nRealizations = 2      

fn_base = "20241030D"
for i in 1:nRealizations
    fn = "$fn_base($i)"
    recon_options["rfFile"]  = fn
    RFdeg = rfFunction(recon_options["nTR"], recon_options["nky"])
    noises, ItotAll, b1f = BLAKJac.BLAKJac_analysis!(cpu, RFdeg, trajectorySet, recon_options, saved_H)
    @show noises;
end 

## (dependent!)
fn_base = "20241101G"
for i in 1:nRealizations
    fn = "$fn_base($i)"
    recon_options["rfFile"]  = fn
    RFdeg = rfFunction(recon_options["nTR"], recon_options["nky"])
    noises, ItotAll, b1f = BLAKJac.BLAKJac_analysis!(cpu, RFdeg, trajectorySet, recon_options, saved_H)
    @show noises;
end 

## 2024-10-30 examine sensitivity of numerical error
aaa = [ 300 -50 100; -300 50 -140; 100 -90 100; -100 90 -80; 
    -220 140 -140; 200 -130 150; -120 100 -80; 120 -100 80; 
     230 -120 70; -240 120 -80; 150 -90 80; -140 80 -70]
perturb=zeros(12,3)
#perturb[11,3]=10^-4
perturb=10^-4*randn(12,3)
@show perturb

h1 = (aaa'*aaa)^-1 
h2 = ((aaa+perturb)'*(aaa+perturb))^-1
@show h1
@show (h2-h1) 


## 2024-11-01 Same as 2024-10-24 'Start simple' but with phase      
include("setup.jl")
recon_options = Dict() # erase all existing settings
nsweeps = 6                                               
nky = 64 # keep it short, for processing time                                                 
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)
rfFunction = rfDictionary["from_file"]
saved_H = Dict()

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.008
recon_options["startstate"] = 1 # -1 
recon_options["sigma_ref"] = 1.4 # See logbook 20220815
recon_options["maxstate"] = 128
recon_options["optpars"]   = Optim.Options(time_limit = 200.0, iterations = 100000, f_tol=1.0e-5, g_tol = 1.0e-5)  ## time_limit!!!
recon_options["opt_criterion"] = "noise_level" 
recon_options["account_SAR"]   = false     
recon_options["emphasize_low_freq"] = true 
recon_options["handleB1"] = "no"       
recon_options["opt_initialize"] = "cRandom30" 
recon_options["opt_focus"] = "max"      
recon_options["opt_complex"] = false      
recon_options["opt_account_maxFlip"] = false
recon_options["opt_keep_positive"] = false                           
recon_options["opt_slow_phase"] =   false # true                         
recon_options["considerCyclic"] = false # true # false  
recon_options["opt_emergeCriterion"] = 500 # 2000
ph = [] 
# ph = zeros(nTR); ph .= 2.0;
recon_options["opt_imposed_2nd_derivative_of_phase"] = ph
recon_options["opt_iterations_limit"] = 1
recon_options["sizeSteps"] = [6]  
nRealizations = 4      

fn_base = "20241101M"
for i in 1:nRealizations; goodseed = i
#for (i,goodseed) in enumerate([5])
    stageText = ""
    portionRange = 0:0
    fn = "/home/mfuderer/Documents/Julia/Capture/$fn_base($i).jld2"
    RFdeg = BLAKJac.BLAKJac_optimize(trajectorySet, recon_options, goodseed);
    FileIO.save(fn,"RFdeg",RFdeg)
    noises, ItotAll, b1f = BLAKJac.BLAKJac_analysis!(cpu, RFdeg.+0im, trajectorySet, recon_options, saved_H)
    @show noises;
end 


## 2024-11-01 Analysis 
include("setup.jl")
recon_options = Dict() # erase all existing settings
nsweeps = 6                                               
nky = 64 # keep it short, for processing time                                                 
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)
rfFunction = rfDictionary["from_file"]
saved_H = Dict()

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.008
recon_options["startstate"] = 1 # -1 
recon_options["sigma_ref"] = 1.4 # See logbook 20220815
recon_options["account_SAR"]   = false     
recon_options["emphasize_low_freq"] = true 
recon_options["handleB1"] = "no"                            
recon_options["considerCyclic"] = false # true # false 
recon_options["maxstate"] = 64

nRealizations = 4      

fn_base = "20241101M"
for i in 1:nRealizations
    fn = "$fn_base($i)"
    recon_options["rfFile"]  = fn
    RFdeg = rfFunction(recon_options["nTR"], recon_options["nky"])
    noises, ItotAll, b1f = BLAKJac.BLAKJac_analysis!(cpu, RFdeg, trajectorySet, recon_options, saved_H)
    @show noises;
end 

## 2024-11-02 (dependent!)
fn_base = "20241101J"
nRealizations = 4
(fig,ax)=(subplots(Int(ceil(nRealizations/4)),4,figsize=(12,3)))
for i in 1:nRealizations
    fn = "$fn_base($i)"
    recon_options["rfFile"]  = fn
    RFdeg = rfFunction(recon_options["nTR"], recon_options["nky"])
    ax[i].plot(abs.(RFdeg))
    anglesdd = zeros(length(RFdeg))
    for i in 1:length(RFdeg)-2
        anglesdd[i] = (rad2deg(angle(conj(RFdeg[i])*RFdeg[i+1]*RFdeg[i+1]*conj(RFdeg[i+2])))+270.0) % 180.0 -90.0
    end
    ax[i].plot(anglesdd)
    noises, ItotAll, b1f = BLAKJac.BLAKJac_analysis!(cpu, RFdeg, trajectorySet, recon_options, saved_H) 
    @show noises, mean(noises), ItotAll, b1f
end











