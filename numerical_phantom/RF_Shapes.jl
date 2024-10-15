# Collection of various RF-shape definitions

function RF_linear(nTR,ny)
    RF = complex.(LinRange(1,90,nTR))
    return RF
end

function RF_linear1070(nTR,ny)
    RF = complex.(LinRange(10,70,nTR))
    return RF
end

function RF_two_segments(nTR,ny)
    RF=complex.(zeros(nTR))
    half = Int32(nTR / 2)
    RF[1:half]    .=complex(90)
    RF[half+1:nTR].=complex(30)
    return RF
end

function RF_partly_smoothly_oddeven(nTR,ny)
    RF=complex.(zeros(nTR))
    half = Int32(nTR / 2)
    # first half has a cos2 shape
    # second half alternates between small and large, ending at 0 / 180

    for i=1:half
        RF[i]      = 90*complex(1+cospi(2*i/half))/2
        RF[i+half] = 90*complex(1+(2*(i%2)-1)*(1-cospi(2*i/half))/2)
    end
    return RF
end

function RF_partly_hard_oddeven(nTR,ny)
    RF=complex.(zeros(nTR))
    half = Int32(nTR / 2)
    # first half is constant at 30
    # second half changes abruptly to 0 / 60 alternation

    for i=1:half
        RF[i]      = 30
        RF[i+half] = 30*complex(2*(i%2))
    end
    return RF
end

function RF_fast_rollercoaster(nTR,ny)
    RF=complex.(zeros(nTR))
    # period in number of TRs
    period = 50
    # flip angle at peak
    amplitude = 80
    base=10

    time = LinRange(1,nTR,nTR)
    @. RF = complex(base+amplitude*(1+cospi(2*time/period))/2)
end

function RF_growing_rollercoaster(nTR,ny)
    RF=complex.(zeros(nTR))
    # period in number of TRs
    period = 20
    # flip angle at peak
    amplitude = 140
    base=10

    time = LinRange(1,nTR,nTR)
    @. RF = complex((time/nTR)*(base+amplitude*(1+cospi(2*time/period))/2))
end

function RF_slow_grow_rollercoaster(nTR,ny)
    RF=complex.(zeros(nTR))
    # period in number of TRs
    period = 22
    # flip angle at peak
    amplitude = 50
    base = 15

    time = LinRange(1,nTR,nTR)
    @. RF = complex((time/nTR)*(base+amplitude*(1+cospi(2*time/period))/2))
    return RF
end

function RF_hiflip_withdip(nTR,ny)
    # only for understanding behavior - what is the proportion of direct signal?
    RF=complex.(zeros(nTR))
    @. RF = 70
    RF[500] = 0
    return RF
end

function RF_flatflatwobble(nTR,ny)
    RF=complex.(zeros(nTR))
    # period in number of TRs
    period = 750/26
    # flip angles, from base to base+2*amplitude
    amplitude = 30
    base = 15

    time = LinRange(1,nTR,nTR)
    @. RF = complex(base+amplitude*(1+cospi(2*time/period)))
    for j=1:260 RF[j]=17.0+0im end
    return RF
end

function RF_peakflatwobble(nTR,ny)
    RF=complex.(zeros(nTR))
    # period in number of TRs
    period = 30.2
    # flip angles, from base to base+2*amplitude
    amplitude = 30
    base = 15

    time = LinRange(1,nTR,nTR)
    @. RF = complex(base+amplitude*(1+cospi(2*time/period)))
    for j=1:260 RF[j]=17.0+0im end
    RF[1]=70.0+0im
    RF[2]=30.0+0im
    RF[3]=22.0+0im
    return RF
end

function RF_insync(nTR,ny)
    RF=complex.(zeros(nTR))
    # period in number of TRs
    period = ny
    # flip angles, from base to base+2*amplitude
    amplitude = 35
    base = 5

    time = LinRange(1,nTR,nTR)
    @. RF = complex(base+amplitude*(1-cospi(2*time/period))+0.1*time/nTR)
    return RF
end

function RF_outofsync(nTR, ny)
    RF=complex.(zeros(nTR))
    # period in number of TRs
    period = Int(1.25 * ny)
    # flip angles, from base to base+2*amplitude
    amplitude = 35
    base = 0

    time = LinRange(1,nTR,nTR)
    @. RF = complex(base+amplitude*(1-cospi(2*time/period))+0.1*time/nTR)
    return RF
end

function RF_10_30_60(nTR,ny)
    # somehow reasoned to be beneficial: flip angle variation between 10 and 30 degrees,
    # with a period of approx. 600ms which is out of sync
    RF=complex.(zeros(nTR))
    # period in number of TRs
    period = nTR/17 # the denominator not being a multiple of 4  makes it out of sync
    # flip angles, from base to base+2*amplitude
    amplitude = 10
    base = 10

    time = LinRange(1,nTR,nTR)
    @. RF = complex(base+amplitude*(1+cospi(2*time/period)))
    return RF
end

function RF_random_10_30(nTR,ny)
    RF=complex.(rand(Float64,nTR))
    RF .*= 20
    RF .+= 10
    return RF
end

function RF_low_flip_grow(nTR,ny)
    RF=complex.(zeros(nTR))
    time = LinRange(1,nTR,nTR)
    # period in number of TRs
    period = nTR/17 # the denominator not being a multiple of 4  makes it out of sync
    # flip angles, from base to base+2*amplitude
    amplitude = 40.0 .- 40 .* (time./nTR)
    base = 4.0 .* (time./nTR)

    @. RF = complex(base+amplitude*(1+cospi(2*time/period)))
    return RF
end


function RF_out_var(nTR,ny)
    RF=complex.(zeros(nTR))
    # period in number of TRs
    period = nTR / 3.0
    # flip angles, from base to base+2*amplitude
    amplitude = 22
    base = 6

    time = LinRange(1,nTR,nTR)
    @. RF = complex(base+amplitude*(1-cospi(2*time/period)))
    return RF
end

function RF_exp_best(nTR,ny)
    RF=complex.(zeros(nTR))
    # The number of periods is chosen to be odd (assuming nTR/nky is even) and lead to a period of approx 2.8s
    periods = floor(0.01*nTR / 2.8 / 2) * 2 -1
    # period in number of TRs
    period = nTR / periods
    # flip angles, from base to base+2*amplitude
    amplitude = 25
    base = 0

    time = LinRange(1,nTR,nTR)
    @. RF = complex(base+amplitude*(1-cospi(2*time/period))+0.1*time/nTR)
    return RF
end

function RF_exp_best_4_over_5(nTR,ny)
    RF=complex.(zeros(nTR))
    @assert(nTR==5*ny) # this function is only meaningful for a 5-measurement "actual" setup
    periods = 4
    # period in number of TRs
    period = nTR / periods
    # flip angles, from base to base+2*amplitude
    amplitude = 25
    base = 0

    time = LinRange(1,nTR,nTR)
    @. RF = complex(base+amplitude*(1-cospi(2*time/period))+0.1*time/nTR)
    return RF
end

function RF_equal_height_lobes(nTR,ny,lobes=5)
    RF=complex.(zeros(nTR))
    # period in number of TRs
    period = nTR / lobes
    # flip angles, from base to base+2*amplitude
    amplitude = 25
    base = 0

    time = LinRange(1,nTR,nTR)
    @. RF = complex(base+amplitude*(1-cospi(2*time/period))+0.1*time/nTR)
    return RF
end


function RF_actual(nTR,ny)
    dir = pwd(); @show dir
    a = readdlm("numerical_phantom/rf_insync_angle_list.txt")
    p = readdlm("numerical_phantom//rf_insync_phase_list.txt") .*pi ./180.0
    RFdeg = a .* exp.(im .*p)
    RF_train = vec(complex.(RFdeg))

    if length(RF_train) < nTR
        repfac = nTR÷length(RF_train) + 1
        RF_train = repeat(RF_train,outer=repfac)
    end
    if length(RF_train) > nTR
        RF_train = RF_train[1:nTR]
    end
    return RF_train
end

function RF_insync_varpeaks(nTR, ny)
        sl = 1 # sweeplength
        nk = 8
        RF_train = sin.(LinRange(0,π,sl*ny)).^2
        # RF_train = LinRange(0,90, nr)
        RF_train = repeat(RF_train,1,nk÷sl)
        Random.seed!(1234);
        RF_train = RF_train .* (80 * rand(1,nk÷sl)) .+ 15

        # make sure first sweep starts from zero otherwise aliasing artefacts?!
        RF_train[1:end÷2,1] .= sin.(LinRange(0,π/2,sl*ny÷2)).^2 .* RF_train[end÷2,1]
        RF_train[1,1] = RF_train[2,1] / 2 # start with alpha over 2 to prevent oscillations
        RF_train = vec(complex.(RF_train))
    return RF_train
end

function RF_almostflat(nTR, ny)
    RF_train = LinRange(80,81,nTR)
    RF_train = vec(complex.(RF_train))
    return RF_train
end

function RF_from_optim(nTR, ny)
    RF=complex.(zeros(nTR))
    @assert(nTR==5*ny) # this function is only meaningful for a 5-measurement "actual" setup
    vars = FileIO.load("/home/mfuderer/Documents/Julia/Capture/RFoptimRes10novFullItNM1.jld2")
    #vars = FileIO.load("/home/mfuderer/Documents/Julia/Capture/RFoptimRes11nov20NMc1.jld2")
    RFshort = vars["RFshort"]
    itp = interpolate(RFshort,BSpline(Cubic(Natural(OnGrid()))))
    RFdegR = itp.(LinRange(1,length(RFshort),nTR))
    RFdegRcrop = [((v>180) ? v-360 : v) for v in RFdegR]
    RF_train = vec(complex.(RFdegRcrop))
    return RF_train
end

function RF_from_optim7pt(nTR, ny)
    RF=complex.(zeros(nTR))
    @assert(nTR==5*ny) # this function is only meaningful for a 5-measurement "actual" setup
    vars = FileIO.load("/home/mfuderer/Documents/Julia/Capture/RFoptimRes0302.jld2")
    RFshort = vars["RFdegR"]
    itp = interpolate(RFshort,BSpline(Cubic(Natural(OnGrid()))))
    RFdegR = itp.(LinRange(1,length(RFshort),nTR))
    RFdegRcrop = [((v>180) ? v-360 : v) for v in RFdegR]
    RF_train = vec(complex.(RFdegRcrop))
    return RF_train
end

function RF_from_optim7pt_edited(nTR, ny) # manual editing to remove inserted sign swaps and inserted 180-deg pulses
    RF=complex.(zeros(nTR))
    @assert(nTR==5*ny) # this function is only meaningful for a 5-measurement "actual" setup
    vars = FileIO.load("/home/mfuderer/Documents/Julia/Capture/RFoptimRes0302.jld2")
    RFshort = vars["RFdegR"]
    itp = interpolate(RFshort,BSpline(Cubic(Natural(OnGrid()))))
    RFdegR = itp.(LinRange(1,length(RFshort),nTR))
    RFdegRcrop = [((v>180) ? v-360 : v) for v in RFdegR]

    #RFdegRcrop[189:200] .= 13.0      # temporary dive into negatives removed
    #RFdegRcrop[317:320] .= 105.0     # remove two 180 pulses + surrounding
    #RFdegRcrop[331:332] .= 80.0      # remove two 180 pulses
    #RFdegRcrop[346:349] .= 50.0      # remove three 180 pulses
    #RFdegRcrop[443:448] .= 30.0      # temporary dive into negatives removed
    #RFdegRcrop[524:529] .= 8.0      # temporary dive into negatives removed
    #RFdegRcrop[693:716] .= 109.0      # jitter on a bulge removed
    #RFdegRcrop[671:719] .= 70.0      # alternatively, the whole bulge removed
    #RFdegRcrop[1:10] .= 0.0             # for reference, remove 10 TRs
    #RFdegRcrop[332:332] .= 80.0         # of a doublet, remove one

    RF_train = vec(complex.(RFdegRcrop))
    return RF_train
end

function RF_full_7_cpx(nTR, ny)
    RF=complex.(zeros(nTR))
    @assert(nTR==5*ny) # this function is only meaningful for a 5-measurement "actual" setup
    #vars = FileIO.load("/home/mfuderer/Documents/Julia/Capture/RFoptimCFullRand_210529.jld2")
    vars = FileIO.load("/home/mfuderer/Documents/Julia/Capture/RFoptimCFullRand_210530H.jld2")
    RFdeg = vars["RFdeg"]
    RF_train = vec(complex.(RFdeg))
    return RF_train
end

function RF_from_file(nTR, ny)
    RF=complex.(zeros(nTR))
    #@assert(nTR==5*ny) # this function is only meaningful for a 5-measurement "actual" setup
    fn = recon_options["rfFile"]
    fPath = string("/home/mfuderer/Documents/Julia/Capture/",fn,".jld2")
    vars = FileIO.load(fPath)
    RFdeg = vars["RFdeg"]
    RF_train = vec(complex.(RFdeg))
    return RF_train
end

function RF_from_txt_file(nTR,ny) 
    fn = recon_options["rfFile"];
    a = readdlm("/home/mfuderer/ownCloud - Fuderer, M. (Miha)@surfdrive.surf.nl/data/rf_angle_list"*fn*".txt")
    phaseFile = "/home/mfuderer/ownCloud - Fuderer, M. (Miha)@surfdrive.surf.nl/data/rf_phase_list"*fn*".txt"
    p = zeros(size(a))
    if isfile(phaseFile)
        p = readdlm(phaseFile) .*pi ./180.0
    end
    RFdeg = a .* exp.(im .*p)
    RF_train = vec(complex.(RFdeg))

    if length(RF_train) < nTR
        repfac = nTR÷length(RF_train) + 1
        RF_train = repeat(RF_train,outer=repfac)
    end
    if length(RF_train) > nTR
        RF_train = RF_train[1:nTR]
    end
    return RF_train
end

function RF_25pts_7_cpx(nTR, ny)
    RF=complex.(zeros(nTR))
    @assert(nTR==5*ny) # this function is only meaningful for a 5-measurement "actual" setup
    vars = FileIO.load("/home/mfuderer/Documents/Julia/Capture/RFoptimC25ptRand_210530.jld2")
    RFdeg = vars["RFdeg"]
    RF_train = vec(complex.(RFdeg))
    return RF_train
end

function RF_fromOptimCpx15(nTR, ny)
    RF=complex.(zeros(nTR))
    @assert(nTR==5*ny) # this function is only meaningful for a 5-measurement "actual" setup
    vars = FileIO.load("/home/mfuderer/Documents/Julia/Capture/zzzzzzzzzzzzzzzzzzz.jld2")
    RFin = vars["RFdeg"]
    itp = interpolate(RFin,BSpline(Cubic(Natural(OnGrid()))))
    RFdeg = itp.(LinRange(1,length(RFin),nTR))
    RF_train = vec(complex.(RFdeg))
    return RF_train
end

function RF_fromOptimOn20(nTR, ny)
    RF=complex.(zeros(nTR))
    @assert(nTR==5*ny) # this function is only meaningful for a 5-measurement "actual" setup
    vars = FileIO.load("/home/mfuderer/Documents/Julia/Capture/RFoptimRes10nov20NMstep1.jld2")
    RFshort = vars["RFshort"]
    itp = interpolate(RFshort,BSpline(Cubic(Natural(OnGrid()))))
    RFdegR = itp.(LinRange(1,length(RFshort),nTR))
    RFdegRcrop = [((v>180) ? v-360 : v) for v in RFdegR]
    RF_train = vec(complex.(RFdegRcrop))
    return RF_train
end

function RF_from_optim_intervened(nTR, ny)
    RF=complex.(zeros(nTR))
    @assert(nTR==5*ny) # this function is only meaningful for a 5-measurement "actual" setup
    vars = FileIO.load("/home/mfuderer/Documents/Julia/Capture/RFoptimRes10novFullItNM1.jld2")
    RFshort = vars["RFshort"]
    itp = interpolate(RFshort,BSpline(Cubic(Natural(OnGrid()))))
    RFdegR = itp.(LinRange(1,length(RFshort),nTR))
    RFdegRcrop = [((v>180) ? v-360 : v) for v in RFdegR]

    RFdegRcrop = abs.(RFdegRcrop)     # This is the intervention

    RF_train = vec(complex.(RFdegRcrop))
    return RF_train
end

function RF_flat(nTR, ny) # only for validation purposes
    RF=complex.(zeros(nTR))
    level = 30
    @. RF = complex(level)
    return RF
end

function RF_random_from_20(nTR, ny)
    # pattern interpolated from 20 real random-level points
    seed = recon_options["random_RF_seed"];
    rng = MersenneTwister(seed)
    level = 35.0;
    variation = 35.0;
    RFbase = level.+variation.*randn(rng, Float64,20)
    RFshort = zeros(nTR)
    itp = interpolate(RFbase[:], BSpline(Cubic(Natural(OnGrid()))))
    RFshort = itp.(LinRange(1,size(RFbase,1),nTR))
    return(complex.(RFshort))
end

rfDictionary = Dict{String, Function}();

rfDictionary["linear"]      = RF_linear
rfDictionary["linear1070"]  = RF_linear1070
rfDictionary["two_segments"]= RF_two_segments
rfDictionary["partly_smoothly_oddeven"]= RF_partly_smoothly_oddeven
rfDictionary["partly_hard_oddeven"]= RF_partly_hard_oddeven
rfDictionary["fast_rollercoaster"]= RF_fast_rollercoaster
rfDictionary["growing_rollercoaster"]= RF_growing_rollercoaster
rfDictionary["slow_grow_rollercoaster"]= RF_slow_grow_rollercoaster
rfDictionary["hiflip_withdip"]= RF_hiflip_withdip
rfDictionary["flatflatwobble"]= RF_flatflatwobble
rfDictionary["peakflatwobble"]= RF_peakflatwobble
rfDictionary["insync"]     = RF_insync
rfDictionary["outofsync"]  = RF_outofsync
rfDictionary["10_30_60"]  = RF_10_30_60
rfDictionary["random_10_30"]  = RF_random_10_30
rfDictionary["low_flip_grow"]  = RF_low_flip_grow
rfDictionary["out_var"]       = RF_out_var
rfDictionary["exp_best"]      = RF_exp_best
rfDictionary["exp_best_4_over_5"]= RF_exp_best_4_over_5
rfDictionary["actual"]        = RF_actual
rfDictionary["insync_varpeaks"] = RF_insync_varpeaks
rfDictionary["almostflat"]    = RF_almostflat
rfDictionary["from_optim"]    = RF_from_optim
rfDictionary["from_optim7pt"] = RF_from_optim7pt
rfDictionary["from_optim7pt_edited"] = RF_from_optim7pt_edited
rfDictionary["full_7_cpx"]    = RF_full_7_cpx
rfDictionary["25pts_7_cpx"]   = RF_25pts_7_cpx
rfDictionary["fromOptimOn20"] = RF_fromOptimOn20
rfDictionary["from_optim_intervened"] = RF_from_optim_intervened
rfDictionary["flat"]          = RF_flat
rfDictionary["random_from_20"]= RF_random_from_20
rfDictionary["from_file"]     = RF_from_file
rfDictionary["from_txt_file"] = RF_from_txt_file
rfDictionary["equal_height_lobes"]= RF_equal_height_lobes