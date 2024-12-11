# Function that simulates a 1D profile, reconstructs using polynomial-recon and compares/plots the results
function WrecSimulRun(rng = MersenneTwister(2))

# simulate object 
    ;
    nrSegments = recon_options["objectSegs"]
    spread = recon_options["simulationSpread"]
    rhoVar = recon_options["RhoVariation"]
    segsT1 = recon_options["simulationT1center"] .* exp.(spread .* randn(rng, Float64,nrSegments))
    segsT2 = recon_options["simulationT2center"] .* exp.(spread .* randn(rng, Float64,nrSegments))
    segsRho = recon_options["maxRho"] .* ((1-rhoVar) .+ rhoVar .*rand(rng, Float64,nrSegments))
    segsMask = ones(nrSegments)
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!hack #presumed to be temporary 2021-07-07
    #segsMask[1:100] .= 0.0; segsMask[121:155] .=0.0; segsMask[173:224] .=0.0
    #segsMask[1:100] .= 0.0; segsMask[111:155] .=0.0; segsMask[166:224] .=0.0
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!hack
    segsRho .*= segsMask

    f = recon_options["simulationT1frequency"]
    if (f>0.0)
        #  overwrite to test on 'critical'frequencies
        segsT1 = [recon_options["simulationT1center"]*(1.0+spread*(sin(f*i*2*pi/nrSegments))) for i in 1:nrSegments]
        segsT2 = [recon_options["simulationT2center"] for i in 1:nrSegments]
        segsRho = [recon_options["maxRho"] for i in 1:nrSegments]  
    end
    
    hifri = false #presumed to be temporary (patches around 2021-07-01)
    if hifri
        segsT2 = [recon_options["simulationT2center"]+0.004*(-1)^i for i in 1:nrSegments]
    end

@show mean(segsT1), mean(segsT2)

    repfac = Int64(round(recon_options["objSize"]/nrSegments))

    simRho = repeat(segsRho, inner=repfac)
    simT1  = repeat(segsT1,  inner=repfac)
    simT2  = repeat(segsT2,  inner=repfac)
    simMask =repeat(segsMask,inner=repfac)


    close("all")

# Make and plot simulated sequence
    nTR = recon_options["nTR"]
    rfFunction = recon_options["rfFunction"]
    RFdeg = rfFunction(recon_options["nTR"], recon_options["nky"])
    RFdegsim = RFdeg .* recon_options["simulatedB1errfactor"]
    RFrad = RFdegsim * (π/180) # convert to radians

    figure();
    plot(RFrad);
    xlabel("nTR");
    ylabel("Flip angle (rad)");
    title("Simulated sequence");

    TR = recon_options["TR"]
    TE = TR/2.01
    TI = (recon_options["startstate"]==1) ? 20.0 : 0.01;
    nky = recon_options["nky"]
    ky = 1:nTR
    # per 2021-03-15, found this line of code unexplainable (but same quirk was in the recon code, so it went well)    
    # ky = (1 .+ ky .% nky )   # replaced by ...   
    ky = (1 .+ (ky.-1) .% nky )     
    plot(ky./nky)

# Assemble SPGR sequence simulator
    nPars = recon_options["nPars"]
    sliceprofile = ones(ComplexF64, nTR, 1) # ideal slice profile
    #sliceprofile = SVector( ComplexF64(1.0))  # ideal slice profile

    spgrSim = BlochSimulators.FISP2D(RFdegsim, sliceprofile,TR,TE,Val(recon_options["maxstate"]),TI)
    spgr    = BlochSimulators.FISP2D(RFdeg,    sliceprofile,TR,TE,Val(recon_options["maxstate"]),TI)
    #surrogate = BlochSimulators.PolynomialSurrogate(spgrSim, recon_options)
    #surrogate = BlochSimulators.PolynomialSurrogate(recon_options)         # hacked for debugging

# Calculate magnetization history for all segments

    smapSim = zeros(nTR, nrSegments);
    for i =1:nrSegments
        T1 = segsT1[i]
        T2 = segsT2[i]
        parameters = [BlochSimulators.T₁T₂(T1,T2)]

        s = simulate(cpu,spgrSim,parameters) 
        @assert (!recon_options["surrogateSimulation"]) "option not implemented here"
        # if recon_options["surrogateSimulation"]
        #     s = simulate(cpu,surrogate,parameters)
        # else 
        #     s = simulate(cpu,spgrSim,parameters) 
        # end

        smapSim[:,i] = -imag.(s)
    end

# Per timepoint: reassemble object, IFT and pick the relevant phase-encode

    simData = Array{ComplexF64,1}(undef, nTR)     # will hold the simulated measurements
    simData.= 0

    for i=1:nTR
        # generate image corresponding to this weighting and transform
        segsSim  = segsRho.* smapSim[i,:]
        simImage = repeat(segsSim, inner=repfac)
        simKspace = ifftshift(ifft(ifftshift(complex.(simImage))))

        # pick the datapoint simulated
        thisky = ky[i]
        simData[i] = simKspace[thisky]
    end

# add noise
    SNR = recon_options["simulatedSNR"]
    appliedNoiseFac = 200.0*0.038/SNR
    noiseData = randn(rng,ComplexF64,nTR)
    # temporary tweak 2021-07-12
            # centerky = nky/2.0
            # centeredky = ky.-centerky
            # tweakFilter = ((0.1 .* centeredky).^2 ./(1.0 .+(0.1 .* centeredky).^2))
            # noiseData .*= tweakFilter;
    simData += appliedNoiseFac./sqrt(nky).*noiseData
# apply model error
    nSweeps = Int(floor(nTR/nky))
    errorFactor = exp.(recon_options["sampleBias"].*randn(rng,nSweeps))
    errorFactor = repeat(errorFactor, inner=nky)
    simData[1:nky*nSweeps] .*= errorFactor 
    figure(); plot(errorFactor)

# end of simulation

# intermezzo: calculate full Jacobian given the simualted object
    calculateFullJacobi = false
    if (calculateFullJacobi)
        parameters = BlochSimulators.T₁T₂.(segsT1,segsT2) |> vec
        nY = recon_options["objSize"]
        nYh = floor(nY//2)
        nUnknowns = nY*nPars
        y = -nYh:nY-nYh-1
        row = repeat(y, outer=nPars)

        fullJ = zeros(ComplexF64,nTR,nUnknowns)
        fullJdf= zeros(ComplexF64,nTR,nUnknowns)
        tmpJ  = zeros(ComplexF64, nrSegments, nPars)
        derivs = simulate_derivatives(cpu,surrogate,parameters) 

        # 2022-01-22
        # from here on, a bunch of hard-coding for generating figures for an IEEE-TMI publication

        for i=1:nTR
            kyim2pi = im*2*pi*ky[i]/nky
            tmpJ[:,1] = derivs.m[i,:].*200
            tmpJ[:,2] = segsRho.*derivs.∂T₁[i,:].*recon_options["T1ref"]
            tmpJ[:,3] = segsRho.*derivs.∂T₂[i,:].*recon_options["T2ref"]
            #tmpJ[:,1] = derivs.m[i,:]
            #tmpJ[:,2] = segsRho.*derivs.∂T₁[i,:]
            #tmpJ[:,3] = segsRho.*derivs.∂T₂[i,:]
            fullJ[i,:] = repeat(vec(tmpJ), inner=repfac).*exp.(kyim2pi.*row)
            tmpJrep = repeat(tmpJ, inner=(repfac,1)).*repeat(exp.(kyim2pi.*y),inner=(1,3))
            tmpJrep = fft(fftshift(tmpJrep,1),1)
            fullJdf[i,:] = vec(tmpJrep)
        end

        fullJdfR = reshape(permutedims(reshape(fullJdf,224,5,224,3),[2 1 4 3]),5*224,3*224)

        figure(); imshow(abs.(fullJ), cmap="gray")
        figure(); imshow(abs.(fullJdf), cmap="gray")
        figure(); imshow(abs.(fullJdfR), cmap="gray")
        figure(); imshow(abs.(fullJdfR[1:40,1:24]), cmap="gray")

        # 2022-01-22 dirty hacks end here

        fullH = fullJ'*fullJ
        dHinv = diag(fullH^-1)
        dHinv = reshape(dHinv,nky,nPars)
        figure(); plot(abs.(dHinv[:,1])); title("full Jacobi on Rho")
        figure(); plot(abs.(dHinv[:,2])); title("full Jacobi on T1")
        figure(); plot(abs.(dHinv[:,3])); title("full Jacobi on T2")
    end

# reconstruct
    rhomap = zeros(size(simRho))
    t1map = zeros(size(simRho))
    t2map = zeros(size(simRho))

    if (recon_options["recon_type"]=="polynomial")
        @assert false "polynomial solver not implemented here"
        #close("all")
        T1ref = recon_options["T1ref"]
        T2ref = recon_options["T2ref"]
        
        outmap, normSquared, history = wrecReconIterTargetless(spgr, simData)
        hhh = 20.0 .* history .+61 
        # outmap = wrecReconIterTargetless2(spgr, simData)  # "2" referring to logbook writeup 2021-09-14
        rhomap = abs.(outmap[:,1])
        t1map = abs.(T1ref.*exp.(outmap[:,2]./rhomap))
        t2map = abs.(T2ref.*exp.(outmap[:,3]./rhomap));
    elseif (recon_options["recon_type"]=="reflective_solver")
        qmaps = wrecReconInvokeReflective(spgr, simData)
        rhomap = qmaps.ρ .* nky      # the factor is for a presumed difference of definition of FFT-direction
        t1map = qmaps.T₁ 
        t2map = qmaps.T₂ 
    else
        ro = recon_options["recon_type"]
        @assert false "unknown option $ro"
    end

# Analyze surrogate mismatch for artefatct-level analysis
    analyzeSurrogateMismatch = false
    if (analyzeSurrogateMismatch)
        surrogateError = zeros(nTR, nrSegments);
        for i =1:nrSegments
            T1 = segsT1[i]
            T2 = segsT2[i]
            parameters = [BlochSimulators.T₁T₂(T1,T2)]
            ssurr = simulate(cpu,surrogate,parameters)
            sspgr = simulate(cpu,spgr,parameters) 
            surrogateError[:,i] = -imag.(ssurr.-sspgr)
        end
        rmsSurrogateError = norm(surrogateError)/norm(smapSim)
        @show rmsSurrogateError
    end

# plot results 

    errRho = [simMask[i]>0 ? (abs(rhomap[i]) - abs(simRho[i]))/abs(simRho[i]) : 0.0 for i in 1:nky]
    errT1  = simMask.*(t1map  .- simT1 )./simT1
    errT2  = simMask.*(t2map  .- simT2 )./simT2
    absErrRho = simMask.*(abs.(rhomap) .- abs.(simRho))
    absErrT1  = simMask.*(t1map  .- simT1 )
    absErrT2  = simMask.*(t2map  .- simT2 )
    errMRho = (mean(rhomap) - mean(simRho))/mean(simRho)
    errMT1  = (mean(t1map) - mean(simT1) )/mean(simT1)
    errMT2  = (mean(t2map)  - mean(simT2) )/mean(simT2)
    figure();title("relative errors")
    plot(abs.(errRho),label="Rho")
    plot(abs.(errT1),label="T1")
    plot(abs.(errT2),label="T2")
    legend()
    @show minimum(abs.(errT2)), maximum(abs.(errT2))

    rmserrT1 = sqrt(mean((abs.(errT1).^2)))
    rmserrT2 = sqrt(mean((abs.(errT2).^2)))
    rmserrRho = sqrt(mean((abs.(errRho).^2)))
    rmsAbsErrT1 = 1000.0 * sqrt(mean((abs.(absErrT1).^2)))
    rmsAbsErrT2 = 1000.0 * sqrt(mean((abs.(absErrT2).^2)))
    rmsAbsErrRho = sqrt(mean((abs.(absErrRho).^2)))
    fom = (rmserrRho+rmserrT1+rmserrT2)/3.0
    @show rmserrRho, rmserrT1, rmserrT2, fom;
    @show abs(errMRho), abs(errMT1), abs(errMT2)

    (fig,(ax1,ax2,ax3))=(subplots(1,3,figsize=(9,2)))
    subplots_adjust(top=0.7)
    
    ax1.plot(abs.(rhomap),label="reconstructed rho")
    ax1.plot(simRho,label="input rho")
    ax1.set_title("Rho map");

    ax2.plot(abs.(t1map),label="reconstructed T1")
    ax2.plot(simT1,label="input T1")
    ax2.set_title("T1 map");

    ax3.plot(abs.(t2map),label="reconstructed T2")
    ax3.plot(simT2,label="input T2")
    ax3.set_title("T2 map");

    modelTxt = recon_options["surrogateSimulation"] ? "Surrotgate model" : "Bloch model"
    legendText = @sprintf("%s, SNR %.0f, %s ----> relative error %4.1f%% (%6.3f, %8.6f, %9.7f)", 
                            recon_options["rfName"], SNR, modelTxt, 100*fom, rmsAbsErrRho, rmsAbsErrT1, rmsAbsErrT2)
    suptitle(legendText)

# some belated results of the "intermezzo"
    if (calculateFullJacobi)
        #for i=1:3; @show sqrt(mean(dHinv[:,i])*nky); end
        fullNoise = [sqrt(mean(abs.(dHinv[:,i]))*nky) for i in 1:3]
            # The factor 2 is to account for being affected by just a single complex component of the noise

        # generate BLAKJac results for comparison
        recon_options["T1T2testset"] = [(segsT1[i],segsT2[i]) for i in 1:nrSegments]  
        H,HN = BLAKJac(RFdeg); 
        varρ = mean(abs.(H[:,:,1,1]))
        varT₁s = [mean(abs.(H[:,i,2,2]))/(segsRho[i])^2 for i in 1:nrSegments]
        varT₂s = [mean(abs.(H[:,i,3,3]))/(segsRho[i])^2 for i in 1:nrSegments]
        varT₁ = mean(varT₁s)
        varT₂ = mean(varT₂s)
        bljNoise = [sqrt(varρ), sqrt(varT₁), sqrt(varT₂)]

        # temporary line for plotting BLAKJac-on-average-T1T2 result
        # @show appliedNoiseFac

        showFullNoise = appliedNoiseFac.*sqrt.(dHinv[:,2].*nky./2.0)
        showBljNoise = repeat(appliedNoiseFac.*sqrt.(varT₁s./2.0), inner=repfac)
            # In above two statements, the factor 2 is to account for being affected by just a single complex component of the noise

        figure(); 
        plot(showFullNoise,label="full-Jacobi estimated deviation on T1")
        plot(showBljNoise,label="BLAKJac estimated deviation on T1")
        plot(0.002.*showBljNoise./showFullNoise, label="Ratio * 0.002")
        # temporary line for plotting BLAKJac-on-average-T1T2 result
        # showBljConst = [0.125398*0.076/sqrt(2.0) for i in 1:nky]; plot(showBljConst,label="BLAKJac assuming average T1,T2")
        legend()

        @show fullNoise
        @show bljNoise
        discrep = log.(bljNoise./fullNoise)
        @show sqrt(mean(discrep.*discrep))
    end

    return fom, t2map, errT2, t1map, errT1;
end


include("../numerical_phantom/load_data_phantom.jl");
include("../numerical_phantom/RF_Shapes.jl");

simData = rand(ComplexF64, 128);

function plotfn_line(x; figtitle="test")

    figure()
    subplot(1,3,1)
        plot(x.T₁)
    subplot(1,3,2)
        plot(x.T₂)
    subplot(1,3,3)
        plot(x.ρˣ)
        plot(x.ρʸ)
end

function wrecReconInvokeReflective(spgr, simData::Array{ComplexF64})

        # Select hardware resource
    resource = has_cuda_gpu() ? CUDALibs() : CPU1()

    # Numerical phantom setup
    recon_options = load_default_recon_options();

    recon_options["description"] = "";
    recon_options["recon_folder"] = "tmp";
    recon_options["recon_subfolder"] = "tmp";

    recon_options["numphantom"] = true;
    recon_options["numphantom_rf_shape"] = "insync_varpeaks";

    recon_options["numphantom_type"] = "line" # "line", "brainweb", "checkerboard", "shepp_logan", "tubes", "flat"
    recon_options["numphantom_size"] = (1,256) # (224,224)
    recon_options["numphantom_sequence"] = "Spoiled" # "Spoiled", "Balanced"
    recon_options["numphantom_trajectory"] = "Cartesian" # "Cartesian", "Radial"
    recon_options["numphantom_noise"] = false

    recon_options["slice_profile_correction"] = "shinnar_leroux" # small_tip_angle or shinnar_leroux
    recon_options["slice_thickness_multiplier"] = 3 # compute slice profiles from 0 to multiplier * nominal slice thickness

    recon_options["slice"] = 1;
    recon_options["coils"] = 1; # should be a range
    recon_options["lsqr_its"] = 10;
    recon_options["trf_max_iter_steihaug"] = 30;
    recon_options["trf_max_iter"] = 3;

    sequence, coordinates, coilmaps, trajectory, mask, phantom = load_data_phantom(recon_options);

    p = [T₁T₂ρˣρʸ(rand(), rand(), 1.0, 0.0) for i = 1:16]
    phantom = repeat(p, inner = 16) |> StructArray


    # output = mrstat_recon(sequence, trajectory, phantom, mask, coordinates, coilmaps, recon_options, CUDALibs(), plotfn_line);

    @assert recon_options["numphantom"]

    vx,vy = recon_options["numphantom_size"]
    parameters = collect(vec(phantom))[mask]

    # fucntion to adapt things to the currently used hardware resource
    adpt(x) = MRSTAT.Reconstruction.adapt_to(resource, x)

    ρ = complex.(phantom.ρˣ[mask], phantom.ρʸ[mask])
    ρ = eltype(adpt(sequence).RF_train).(ρ)
    echos = MRSTAT.BlochSimulators.simulate(resource, adpt(sequence), adpt(parameters))

    @info "Simulate data"
    @time raw_data = Mv(resource, echos, adpt(parameters), adpt(coordinates), adpt(coilmaps), adpt(trajectory), ρ);

    if vx > 1
        mask = reshape(mask,vx,vy);
    end

    B₁map = complex.(ones(size(mask))) #|> vec

    @warn "Proton density will not be initialized with LSQR for numerical phantoms"
    recon_options["initialize_pd"] = false
    recon_options["generate_mask"] = false

    qplot(phantom, figtitle = "Ground truth")

    if recon_options["numphantom_noise"]
        # add noise to simulated data
        max_signal = maximum([abs.(x[1]) for x in raw_data]);
        noise = max_signal * 2 * 0.0002 * randn(ComplexF64, length(raw_data)) # just some factor to get ok-ish noise level

        # noise = deserialize("/home/imaging_mrstat/ovanderheide/mrstat/noise")
        # noise = deserialize("noise_cartesian")
        # noise = load("tmp/numphantom/noisy_slr/mrstat_output.jld2", "mrstat_output").recon_options["numphantom_noise_vector"]

        raw_data = [raw_data[i] .+ noise[i] for i in eachindex(raw_data)]
        recon_options["numphantom_noise_std"] = std(noise)
        recon_options["numphantom_noise_vector"] = noise

        @show 100 * (std(noise) / max_signal)
        # serialize("noise_cartesian", noise)
    end

    output = mrstat_recon(sequence, trajectory, raw_data, mask, coordinates, coilmaps, B₁map, recon_options, resource, qplot, recon_options["recon_subfolder"])


    return output.qmaps[end] # return final iteration
end


## 2023-06-23 Test script by Oscar
# test diffusion code
function thisIsAtest()
include("../numerical_phantom/load_data_phantom.jl");
include("../numerical_phantom/RF_Shapes.jl");

# Select hardware resource
resource = has_cuda_gpu() ? CUDALibs() : CPU1()

# Numerical phantom setup
recon_options = load_default_recon_options();

recon_options["numphantom"] = true;
recon_options["numphantom_rf_shape"] = "insync_varpeaks";

recon_options["numphantom_type"] = "line" # "line", "brainweb", "checkerboard", "shepp_logan", "tubes", "flat"
recon_options["numphantom_size"] = (1,256) # (224,224)
recon_options["numphantom_sequence"] = "Spoiled" # "Spoiled", "Balanced"
recon_options["numphantom_trajectory"] = "Cartesian" # "Cartesian", "Radial"
recon_options["numphantom_noise"] = false

recon_options["slice_profile_correction"] = "shinnar_leroux" # small_tip_angle or shinnar_leroux
recon_options["slice_thickness_multiplier"] = 3 # compute slice profiles from 0 to multiplier * nominal slice thickness


sequence, coordinates, coilmaps, trajectory, mask, phantom = load_data_phantom(recon_options);
parameters = collect(vec(phantom));

parameters = [T₁T₂(rand(), rand()) for i = 1:length(coordinates)]

parameters_with_diffusion = [T₁T₂D(p.T₁, p.T₂, 0.1) for p in parameters]

simulate(CPU1(), sequence, parameters)
simulate(CPU1(), sequence, parameters_with_diffusion)

sequence_d = MRSTAT.Reconstruction.adapt_to(CUDALibs(), sequence)
parameters_d = MRSTAT.Reconstruction.adapt_to(CUDALibs(), parameters)
parameters_with_diffusion_d = MRSTAT.Reconstruction.adapt_to(CUDALibs(), parameters_with_diffusion)

@time simulate(CUDALibs(), sequence_d, parameters_d)
@time simulate(CUDALibs(), sequence_d, parameters_with_diffusion_d)
end


