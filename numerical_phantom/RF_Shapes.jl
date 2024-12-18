# Collection of various RF-shape definitions

function RF_linear(nTR,ny)
    RF = complex.(LinRange(1,90,nTR))
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

function RF_from_file(nTR, ny)
    RF=complex.(zeros(nTR))
    #@assert(nTR==5*ny) # this function is only meaningful for a 5-measurement "actual" setup
    fn = recon_options["rfFile"]
    fPath = string(recon_options["rfFolder"],fn,".jld2")
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

function RF_flat(nTR, ny) # only for validation purposes
    RF=complex.(zeros(nTR))
    level = 30
    @. RF = complex(level)
    return RF
end

rfDictionary = Dict{String, Function}();

rfDictionary["linear"]      = RF_linear
rfDictionary["insync"]     = RF_insync
rfDictionary["actual"]        = RF_actual
rfDictionary["flat"]          = RF_flat
rfDictionary["from_file"]     = RF_from_file
rfDictionary["from_txt_file"] = RF_from_txt_file
